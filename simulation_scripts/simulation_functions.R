## Simulation functions

# Load libraries
library(tidyverse)
library(mvtnorm)
library(viridis)
library(scales)
library(WeightIt)
library(mgcv)
library("GPCERF", lib.loc = "/n/home_fasse/mcork/R/ifxrstudio/RELEASE_3_16")
library(SuperLearner)
library(stats)
library(spdep)
library(sf)
library(ape)

# Load parallel packages (currently not using)
library(future)
library(furrr)
library(parallel)


# Function to simulate data for our simulation
simulate_data <- function(n = 1000, sigma = 1, exp_relationship = "linear") {
  x <- rnorm(n, sd = 1)
  y <- rnorm(n, sd = 1)
  v <- rnorm(n, sd = 1)
  dat <- tibble(x, y, v)
  d <- as.matrix(dist(dat[, c('x', 'y')]))
  cov_mat <- exp(-d^2 / sigma)
  
  # Check for negative eigenvalues
  eigenvalues <- eigen(cov_mat)$values
  if(any(eigenvalues < 0)) message("Negative eigenvalue detected in covariance matrix")
  
  # Spatial confounder s, sampled from multivariate normal distribution
  # for now multiply by 10 
  s <- mvtnorm::rmvnorm(1, sigma = 2 * cov_mat)

  # Calculate Moran's I statistic to see how spatially correlated it is
  location_dists_inv <- 1/d
  diag(location_dists_inv) <- 0
  Moran.I(as.vector(s), location_dists_inv)
  
  
  if (exp_relationship == "linear") {
    w <- 0.5 * v + 1.5 * s + rnorm(n, sd = sqrt(5)) + 5
  } else if (exp_relationship == "nonlinear") {
    w <- 3 * cos(s) + exp(0.1 * v) + 5 + rnorm(n, sd = 2)
  } else {
    stop("Must specify confounder-exposure relationship as linear or nonlinear")
  }
  
  #w <- rescale(w, to = c(0, 10))
  
  #w <- rescale(40*(1+exp(-0.8+0.1*s + 0.1*v))^(-1) - 18 + rnorm(n, sd = 2), to = c(0, 10))
  
  # Generate outcome
  outcome <- s^2 - 3 * sqrt(abs(v)) +  w*(0.1 - s + 0.1*v) + 
    0.01*w^3 + rnorm(n, sd = sqrt(10)) + 10
  
  # Add the variables w and outcome to the data frame
  dat$w <- w[, ]
  dat$outcome <- outcome[, ]
  dat$s <- s[, ]
  
  # Now I make some plots to assess how correlated things are
  ggplot(dat) + 
    geom_point(aes(x = x, y = y, color = s)) + 
    scale_color_viridis()
  
  # Now plot to see what true exposure response looks like
  map_dfr(seq(0, 10, length = 100), function(pot_exp){
    df <- dat
    df$w <- pot_exp
    truth = mean(s^2 - 3 * sqrt(abs(df$s)) + 
      df$w*(0.1 - df$s + 0.1*df$v) + 
               0.01*(df$w)^3 + 10)
    return(tibble(w = pot_exp, truth = truth))
  }) %>% 
  ggplot() + 
  geom_point(aes(x = w, y = truth))
  

  return(dat)
}

fit_models <- function(dat, n_core, spatial_hidden, w_values) {
  
  # Should the spatial confounder be included or hidden from the model
  if (spatial_hidden) {
    cov_names = "v" # data does not see "s" the spatial covariate
  } else {
    cov_names = c("s", "v")
  }
  
  lm_fit <- lm(reformulate(c(cov_names, "w"), response = "outcome"), data = dat)
  gam_fit <- mgcv::gam(reformulate(c(cov_names, "s(w, k = 6)"), response = "outcome"), data = dat)
  
  # Calculate weights for causal inference method
  ent_weight <- calculate_entropy_weighting(dat, cov_names)
  
  # Make correlation table for entropy weighting
  post_cor <- cov.wt(dat[, c("w", cov_names)], wt = ent_weight$weights, cor = T)$cor
  post_cor <- abs(post_cor[-1, 1])
  
  pre_cor <- cov.wt(dat[, c("w", cov_names)], cor = T)$cor
  pre_cor <- abs(pre_cor[-1, 1])
  
  correlation_table <- tibble(covariate = cov_names, 
                              pre_cor = pre_cor, 
                              post_cor = post_cor, 
                              method = "energy")
  
  # Fit weighted linear model
  lm_ent_weights <- lm(reformulate(c(cov_names, "w"), response = "outcome"), data = dat, weights = ent_weight$weights)
  
  # Fit weighted GAM model
  gam_ent_weights <- mgcv::gam(reformulate(c(cov_names, "s(w, k = 6)"), response = "outcome"), 
                               data = dat, weights = ent_weight$weights)
  
  # Get the gpcerf result by calling the separated function
  gpcerf <- fit_gpcerf(dat, n_core, cov_names, w_values = w_values)
  
  # Add GPCERF covariate balance to correlation table
  cor_table_gpcerf <- 
    tibble(covariate = cov_names, 
           pre_cor = gpcerf$cb_org, 
           post_cor = gpcerf$cb, 
           method = "gpcerf")
  
  correlation_table <- rbind(correlation_table, cor_table_gpcerf)
  
  return(list(lm_fit = lm_fit, gam_fit = gam_fit, 
              lm_ent_weights = lm_ent_weights, gam_ent_weights = gam_ent_weights,
              gpcerf = gpcerf, correlation_table = correlation_table))
}


# Function for calculating entropy weighting 
calculate_entropy_weighting <- function(dat, cov_names) {
  ent_weight <- WeightIt::weightit(reformulate(cov_names, response = "w"), method = "energy", data = dat)
  ent_weight$weights[ent_weight$weights > quantile(ent_weight$weights, 0.99)] <- quantile(ent_weight$weights, 0.99)
  
  # Optional visualization code
  # Histogram ent weights
  #hist(ent_weight$weights)
  
  # Check balance
  # ent_balance_check <- 
  #   bal.tab(w ~ s + v,
  #           data = dat,
  #           weights = ent_weight$weights,
  #           method = "weighting",
  #           stats = c("cor"),
  #           un = T,
  #           abs = T,
  #           thresholds = c(cor = .1), poly = 1)
  # 
  # plot(ent_balance_check)
  
  return(ent_weight)
}

# helper functions for gpcerf function
m_xgboost <- function(nthread = 10, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}

m_ranger <- function(num.threads = 10, ...){
  SuperLearner::SL.ranger(num.threads = num.threads, ...)
}

m_lm <- function(num.threads = 10, ...){
  SuperLearner::SL.lm(num.threads = num.threads, ...)
}

# Function to fit GPCERF method
fit_gpcerf <- function(dat, n_core, cov_names,
                       w_values = seq(min(dat$w), max(dat$w), length.out = 100),
                       nn = F) {
  
  gps_m <- estimate_gps(cov_mt = dat[, cov_names, drop = F],
                        w_all = dat$w,
                        sl_lib = c("m_xgboost", "m_ranger", "m_lm"),
                        dnorm_log = TRUE)
  
  # Set hyperparameters
  if (nn) {
    params_lst <- list(alpha = seq(1, 4, length.out = 8),
                       beta = seq(1, 4, length.out = 8),
                       g_sigma = c(0.1, 1, 10),
                       tune_app = "all",
                       n_neighbor = 20,
                       block_size = 1e4)
    
    # Fit nngpcerf (approximation helpful for now)
    gpcerf <- 
      estimate_cerf_nngp(data = data.frame(select(dat, outcome, w, !!!cov_names)),
                         w = w_values,
                         gps_m = gps_m,
                         params = params_lst,
                         outcome_col = "outcome", 
                         treatment_col = "w", 
                         covariates_col = cov_names,
                         nthread = n_core)
    
  } else if (nn == F) {
    params_lst <- list(alpha = seq(1, 4, length.out = 8),
                       beta = seq(1, 4, length.out = 8),
                       g_sigma = c(0.1, 1, 10),
                       tune_app = "all")
    
    gpcerf <- 
      estimate_cerf_gp(data = data.frame(select(dat, outcome, w, !!!cov_names)),
                         w = w_values,
                         gps_m = gps_m,
                         params = params_lst,
                         outcome_col = "outcome", 
                         treatment_col = "w", 
                         covariates_col = cov_names,
                         nthread = n_core)
  } else {
    stop("Must specify either nearest neighbor (nn = T) or GP (nn = F)")
  }
    
 
  return(gpcerf)
}

# Helper function to predict potential outcomes
predict_outcome <- function(dat, lm_fit, gam_fit, lm_ent_weights, 
                            gam_ent_weights, w_values) {
  
  potential_outcomes <- bind_rows(lapply(w_values, function(pot_exp) {
    
    # Keep both s and v since true relationship is based on both, 
    # But if spatial hidden the models will not use that data
    potential_data <- select(dat, s, v) %>%
      mutate(w = pot_exp)
    
    potential_outcome <- potential_data %>%
      mutate(
        linear_model = predict(lm_fit, potential_data, type = "response"),
        gam_model = predict(gam_fit, newdata = potential_data, type = "response"),
        linear_ent = predict(lm_ent_weights, potential_data, type = "response"),
        gam_ent = predict(gam_ent_weights, newdata = potential_data, type = "response"),
        true_relationship = potential_data$s^2 - 
          3 * sqrt(abs(potential_data$v)) + 
          potential_data$w*(0.1 - potential_data$s + 0.1*potential_data$v) + 
          0.01*(potential_data$w)^3 + 10
      ) %>%
      dplyr::select(w, true_relationship, linear_model, linear_ent, gam_model, gam_ent) %>% 
      summarize_all(mean)
    
    return(potential_outcome)
  }))
  
  return(potential_outcomes)
}

# Function to calculate metrics from fit
calculate_metrics <- function(dat, lm_fit, gam_fit, lm_ent_weights, 
                              gam_ent_weights, gpcerf, 
                              w_values = seq(min(dat$w), max(dat$w), length.out = 100)) {
  
  # Calculate predictions
  data_prediction <- predict_outcome(dat, lm_fit, gam_fit, lm_ent_weights, gam_ent_weights, 
                                     w_values)
  data_prediction$gpcerf <- gpcerf$posterior$mean # Add gpcerf
  
  # Calculate midpoint prediction
  # prediction_midpoint <- predict_outcome(dat, lm_fit, gam_fit, lm_ent_weights, gam_ent_weights,
  #                                        5)
  # prediction_midpoint$gpcerf <- median(gpcerf$posterior$mean) # Add median for gpcerf midpoint
  
  # Adjust predictions based on midpoint
  # prediction_adjusted <- bind_rows(apply(data_prediction[,-1], 1, function(x) {
  #   x - prediction_midpoint[,-1]
  # }))
  # prediction_adjusted <- cbind(w = data_prediction$w, prediction_adjusted)
  
  model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "gpcerf")
  
  # Calculate data metrics on non adjusted prediction for now (after looking into it)
  data_metrics <- 
    data_prediction %>%
    tidyr::gather(all_of(model_types), key = "model", value = "prediction") %>%
    group_by(model, w) %>%
    summarise(
      bias = prediction - true_relationship,
      mse = (prediction - true_relationship) ^ 2
    )
  
  return(list(data_metrics = data_metrics, data_prediction = data_prediction))
}

# Function to interpolate on a common exposure grid for plotting
interpolate_on_grid <- function(df_list, common_exposure_grid) {
  df_list %>%
    map_dfr(function(df) {
      df %>%
        group_by(model) %>%
        do(data.frame(
          w = common_exposure_grid,
          model = unique(.$model),
          prediction = approx(.$w, .$response, xout = common_exposure_grid)$y))
    })
}
