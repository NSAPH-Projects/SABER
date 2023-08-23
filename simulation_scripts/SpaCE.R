rm(list = ls())
library(reticulate)
library(geojsonsf) # convert between geojson and sf
library(sf)
library(zipcodeR)
library(tidyverse)
library(WeightIt)
library(cobalt)
library(data.table)
library(CausalGPS)
library(parallel)
library(GPCERF)

#py_install("git+https://github.com/NSAPH-Projects/space") # Install space bench package

# Name directories
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/GPCERF_spatial/"
data_folder <- paste0(proj_dir, "/GPCERF_spatial/cork_data/")

# Import the SpaceEnv class from the spacebench module
spacebench <- import("spacebench")

# Create an instance of the SpaceEnv class
env <- spacebench$SpaceEnv('healthd_pollutn_mortality_cont')

# Call the make() method
dataset <- env$make()

# Now extract from the dataset to see how this works
# treatment <- dataset$treatment # Not PM anymore but it is a treatment value
# outcome <- dataset$outcome # Outcome is number of deaths
# covariates <- dataset$covariates # Covariates is covariate values for outcome, there are 100 values
# counterfactual <- dataset$counterfactuals # Counterfactual for each with 100 different values
treatment_values <- dataset$treatment_values # this is the values for the counterfactual observed

# adjacency_matrix <- dataset$adjacency_matrix()
# erf <- dataset$erf()

# Use the synthetic ones and then use the prediction ones 
# Data from 2010 from CDC wonder to see data
# Confounding score is legacy one, the confounding score erf 

# # Make a plot of the true counterfactual
# counter_values2 <- dataset$erf() # Makes an ERF 
# counter_values <- colMeans(counterfactual)
# data.frame(pm = treatment_values, mort = counter_values) %>% 
#   ggplot() + 
#   geom_line(aes(x = pm, y =mort )) + 
#   labs(x = "PM2.5", y = "Mortality")
# 
# # Plot data used for association
# data.frame(treatment = treatment, outcome = outcome) %>% 
#   ggplot(aes(x = treatment, y =outcome )) + 
#   geom_point() + 
#   geom_smooth()
#   labs(x = "PM2.5", y = "Mortality")


## Now add the data from Mauricio, full dataset without dropped confoudners (yet)
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/GPCERF_spatial/"
data_folder <- paste0(proj_dir, "/GPCERF_spatial/cork_data/")

# Load in DF
df <- read_csv(paste0(data_folder, "synthetic_data_space.csv"))

# Names of the counterfactual columns in the dataset
counterfactural_col <- paste0("Y_synth_", c("00", "01", "02", "03", "04", 
                                            "05", "06", "07", "08", "09", 10:99))

# grab counterfactuals (over 100 concentrations)
counterfactual <- df %>% select(!!!counterfactural_col)

# Get mean of counterfactuals for ERF
counter_values <- colMeans(counterfactual)

# Create plot of true ERF
data.frame(pm = treatment_values, mort = counter_values) %>% 
  ggplot() + 
  geom_line(aes(x = pm, y =mort )) + 
  labs(x = "treatment", y = "outcome")

# Now calculate the treatment, outcome, and covariates used in study
outcome <- df$Y_synth
treatment <- df$qd_mean_pm25
covariates <- df[, 3:35]
cov_names <- names(covariates)

# Histogram of the outcome
hist(outcome)

# Create new dataset with outcome, treatment, and covariates
df2 <- cbind(outcome, treatment, covariates)

#df2_trimmed <- df2 %>% filter(between(treatment, quantile(treatment, 0.05), quantile(treatment, 0.95)))
# Perform entropy weighting
# Create data entropy matrix
data_ent_matrix <- model.matrix(reformulate(c("-1", cov_names)),
                                data = df2)

# Now scale to be between -1 and 1 for model to fit
c_mat <- scale(data_ent_matrix,
               center = T,
               scale = apply(data_ent_matrix, 2, function(x) ifelse(max(abs(x)) == 0, 1, max(abs(x)))))

# Run entropy weighting algorithm
source("/n/dominici_nsaph_l3/Lab/projects/ERC_simulation/Simulation_studies/data_application/functions/entropy_wt_functions.R")
e <- ebal(df2$treatment, c_mat)
ent_weights <- e$weights

# Truncate 99%
ent_weights[ent_weights > quantile(ent_weights, 0.95)] <- quantile(ent_weights, 0.95)

# Run through one more iteration to try to get rid of extreme weights
e <- ebal(df2$treatment, c_mat, base_weights = ent_weights)
ent_weights <- e$weights
ent_weights[ent_weights > quantile(ent_weights, 0.95)] <- quantile(ent_weights, 0.95)

ent_balance_check <- 
  bal.tab(reformulate(cov_names, response = "treatment"),
          data = df2,
          weights = ent_weights,
          method = "weighting",
          stats = c("cor"),
          un = T,
          continuous = "std",
          s.d.denom = "weighted",
          abs = T,
          thresholds = c(cor = .1), poly = 1)
plot(ent_balance_check)

# Add columns of entropy weights
df2$ent_weights <- ent_weights

number_of_bins <- 10
quantiles <- quantile(df2$treatment, probs = seq(0, 1, length.out = number_of_bins + 1))
bin_labels <- paste(round(quantiles[-length(quantiles)], 1), round(quantiles[-1], 1), sep = " - ")
df2$bins <- cut(df2$treatment, breaks = quantiles, labels = bin_labels, include.lowest = TRUE)

data_summary <- 
  df2 %>%
  group_by(bins) %>%
  summarise(sum_weight = sum(ent_weights),
            avg_outcome = mean(outcome))

ggplot(data_summary, aes(x = factor(bins))) +
  geom_bar(aes(y = sum_weight), stat = "identity", fill = "blue", alpha = 0.7) +
  geom_line(aes(y = avg_outcome * 10), group = 1, color = "red") +
  geom_point(aes(y = avg_outcome * 10), color = "red") +
  scale_y_continuous(sec.axis = sec_axis(~./10, name = "Average Outcome")) +
  labs(x = "Binned Treatment Values", y = "Sum of Weights", title = "Sum of Weights and Average Outcome by Treatment Bins") +
  scale_x_discrete(labels = bin_labels) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Now fit the entropy based weights
linear_fit <- 
  lm(reformulate(c("treatment", cov_names), response = "outcome"),
     data = df2)

# Fit gam entropy
gam_fit <-
  mgcv::bam(reformulate(c("s(treatment, bs = 'cr', k = 4)", cov_names), response = "outcome"),
            data = df2)

gam_fit <- 
  mgcv::gam(reformulate(c("s(treatment, k = 4)", cov_names), response = "outcome"),
          data = df2)

# Fit linear model 
linear_ent <- 
  lm(reformulate(c("treatment", cov_names), response = "outcome"),
     data = df2, weight = df2$ent_weights)

data_prediction <- 
  rbindlist(lapply(treatment_values, function(pot_exp) {
    
    # Get potential data if all had same potential exposure
    potential_data <- 
      select(df2, all_of(cov_names)) %>% 
      mutate(treatment = pot_exp)
    
    # Fit each model and take mean for ERC
    potential_outcome <- 
      potential_data %>% 
      mutate(
        linear_model = predict(linear_fit, potential_data, type = "response"),
        linear_ent = predict(linear_ent, potential_data, type = "response")) %>%
      dplyr::select(treatment, linear_model, linear_ent) %>% 
      summarize_all(mean) %>% 
      data.table()
    return(potential_outcome)
  }))

# Specify model types
model_types <- c("linear_model", "linear_ent")

# plot all data together
cbind(data_prediction, true_counterfactual = counter_values) %>% 
  pivot_longer(c(all_of(model_types), "true_counterfactual"), names_to = "model", values_to = "prediction") %>%
  arrange(treatment) %>%
  ggplot() + 
  geom_line(aes(x = treatment, y = prediction, color = model)) + 
  labs(x = "PM2.5", y = "Mortality")


influence_measures <- influence.measures(linear_ent)
summary(influence_measures) # Summary of influence measures
influence_measures$infmat   # Detailed influence measures
plot(influence_measures, "cook") # Plot Cook's distance

cooks_d <- cooks.distance(linear_ent)
plot(cooks_d, ylab="Cook's distance", xlab="Observation", main="Influence Plot of Cook's Distance")
abline(h = 4/length(cooks_d), col="red")

# Look into the top 5 influential points
top_5_influential <- order(cooks_d, decreasing = TRUE)[1:10]
View(df2[top_5_influential, ] )



high_cook <- which(influence_measures$infmat[, "cook.d"] > 4 * mean(influence_measures$infmat[, "cook.d"]))
df2[high_cook, ]
plot(influence_measures, "cook") 

# Fit GAM model
gam_ent <- 
  mgcv::gam(reformulate(c("s(treatment, bs = 'cr', k = 4)", cov_names), response = "outcome"),
            data = df2, weight = df2$ent_weights)

gam_ent <- 
  mgcv::gam(reformulate(c("s(treatment, k = 4)", cov_names), response = "outcome"),
            data = df2, weight = df2$ent_weights)


# now fit causalGPS 
tune_grid <- expand.grid(
  nrounds = c(100),
  eta = c(0.025, 0.05, 0.075, 0.1, 0.15),
  max_depth = c(2),
  delta = seq(0.5, 1.5, by = 0.05)
)

# Make list to pass to parLapply
tune_grid_list <- as.list(as.data.frame(t(tune_grid)))

causal_df <- cbind(id = 1:3109, df2)

# Wrapper function for running causalGPS
wrapper_func <- function(tune_param){
  pseudo_pop_tune <- generate_pseudo_pop(Y = causal_df[, c("id", "outcome")],
                                         w = causal_df[, c("id", "treatment")],
                                         c = data.frame(select(causal_df, !!! c("id", cov_names))),
                                         ci_appr = "matching",
                                         pred_model = "sl",
                                         gps_model = "parametric",
                                         use_cov_transform = FALSE,
                                         transformers = list("pow2", "pow3", "abs", "scale"),
                                         expos_trim_qlts = c(0,1),
                                         gps_trim_qlts = c(0, 0.99),
                                         optimized_compile = T,
                                         sl_lib = c("m_xgboost"),
                                         params = list(xgb_rounds = tune_param[[1]],
                                                       xgb_eta = tune_param[[2]],
                                                       xgb_max_depth = tune_param[[3]]),
                                         covar_bl_method = "absolute",
                                         covar_bl_trs = 0.1,
                                         covar_bl_trs_type = "mean",
                                         max_attempt = 1,
                                         dist_measure = "l1",
                                         delta_n = tune_param[[4]],
                                         scale = 1,
                                         nthread = 1)
  
  
  matched_pop_tune <- pseudo_pop_tune$pseudo_pop
  
  # Truncate upper 1%
  matched_pop_tune <- 
    matched_pop_tune %>%
    mutate(counter_weight = if_else(counter_weight > quantile(counter_weight, 0.99),
                                    quantile(counter_weight, 0.99), 
                                    counter_weight))
  
  # Now generate covariate balance tab
  balance_table <- 
    bal.tab(reformulate(cov_names, response = "treatment"),
            data = matched_pop_tune,
            weights = matched_pop_tune$counter_weight,
            method = "weighting",
            stats = c("cor"),
            un = T,
            continuous = "std",
            s.d.denom = "weighted",
            abs = T,
            thresholds = c(cor = .1), poly = 1)
  
  mean_post_cor <- mean(balance_table$Balance$Corr.Adj)
  
  results <- 
    data.frame(nrounds = tune_param[[1]], 
               eta = tune_param[[2]],
               max_depth = tune_param[[3]],
               delta = tune_param[[4]],
               scale = 1,
               mean_corr = mean_post_cor,
               max_corr = max(balance_table$Balance$Corr.Adj))
  return(list(results, matched_pop_tune))
  
}

# Now iterate through simulation function, bind together results
pseudo_pop_list <- mclapply(tune_grid_list, mc.cores = 10, wrapper_func)
corr_search <- do.call("rbind", (lapply(pseudo_pop_list, function(x) x[[1]])))
min_corr = corr_search %>% filter(mean_corr == min(mean_corr)) %>% slice_head()

# Extract minimum as your result
pseudo_pop_tuned <- wrapper_func(min_corr[1:4])[[2]]

# Check that correlation
post_cor <- cov.wt(select(pseudo_pop_tuned, !!! c("treatment", cov_names)), wt = (pseudo_pop_tuned$counter_weight), cor = T)$cor
post_cor <- abs(post_cor[-1, 1])

# Now fit semi-parametric here 
causal_gps_fit <- 
  mgcv::gam(formula = reformulate(c("s(treatment, bs = 'cr', k = 4)", cov_names), response = "outcome"),
            family = "gaussian",
            data = data.frame(pseudo_pop_tuned),
            weights = counter_weight)


data_prediction <- 
  rbindlist(lapply(treatment_values, function(pot_exp) {
    
    # Get potential data if all had same potential exposure
    potential_data <- 
      select(df2, all_of(cov_names)) %>% 
      mutate(treatment = pot_exp)
    
    # Fit each model and take mean for ERC
    potential_outcome <- 
      potential_data %>% 
      mutate(
        linear_model = predict(linear_fit, potential_data, type = "response"),
        gam_model = predict(gam_fit, newdata = potential_data, type = "response"),
        linear_ent = predict(linear_ent, potential_data, type = "response"),
        gam_ent = predict(gam_ent, newdata = potential_data, type = "response")
        #causal_gps = predict(causal_gps_fit, newdata = potential_data, type = "response")
      ) %>% 
      dplyr::select(treatment, linear_model, linear_ent, gam_model, gam_ent) %>% 
      summarize_all(mean) %>% 
      data.table()

    return(potential_outcome)
  }))


# Specify model types
model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "gpcerf")

# plot all data together
cbind(data_prediction, true_counterfactual = counter_values, gpcerf = cerf_nngp_obj$posterior$mean) %>% 
  pivot_longer(c(all_of(model_types), "true_counterfactual"), names_to = "model", values_to = "prediction") %>%
  arrange(treatment) %>%
  ggplot() + 
  geom_line(aes(x = treatment, y = prediction, color = model)) + 
  labs(x = "PM2.5", y = "Mortality")


# Now work on fitting Boyu's approach
n_core <- 12

# Define SuperLearner wrappers for xgboost and ranger with thread number specified
m_xgboost <- function(nthread = n_core, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}

m_ranger <- function(num.threads = n_core, ...) {
  SuperLearner::SL.ranger(num.threads = num.threads, ...)
}


# Estimate the GPS used in Boyu approach (could also be using similar to CausalGPS package)
new_data <- F
run_tag <- "first_run"
data_dir <- pasteo(proj_dir, "cork_data/SpaCE/", run_tag)

# Load GPS estimates (or run again if new dataset)

if (new_data) {
  # Estimation of Generalized Propensity Score using SuperLearner
  gps_m <- GPCERF::estimate_gps(
    cov_mt = select(df2, !!!cov_names), 
    w_all = df2$treatment,
    sl_lib = c("m_xgboost", "m_ranger"),
    dnorm_log = TRUE
  )
  
  # Safe because this step takes awhile
  saveRDS(gps_m, paste0(data_dir, "gps_m.RDS"))
} else {
  gps_m <- readRDS(paste0(data_dir, "gps_m.RDS"))
}



# Set hyperparameters
params_lst <- list(alpha = 10 ^ seq(0, 1.4, length.out = 10),
                   beta = 10 ^ seq(0, 1.4, length.out = 10),
                   g_sigma = c(0.1, 1, 10),
                   tune_app = "all",
                   n_neighbor = 50,
                   block_size = 1e3)

# Fit nnGPS
cerf_nngp_obj <- 
  estimate_cerf_nngp(data.frame(df2),
                     treatment_values,
                     gps_m,
                     params = params_lst,
                     nthread = 12)
saveRDS(cerf_nngp_obj, paste0(data_dir, "cerf_nngp_obj"))

# Now look into plot
summary(cerf_nngp_obj)
plot(cerf_nngp_obj)



##### Extra code for when getting to SF files
# Read in shapefile for the united state, this is from dataverse 
sf <- geojson_sf(paste0(data_folder, "counties.geojson"), expand_geometries = T)
ggplot() +
  geom_sf(data = sf, fill = NA, color = "black") +
  theme_minimal()  

# # Zip codes arent aligned with the zip code package, but no matter. 
# zipcode_states <- reverse_zipcode(unique(sf$GEOID10))
# 
# # Now fit some naive estimators to the data
# covariates = data.frame(covariates)
# cov_names = paste0("X", 1:33)
# dataset <- cbind(outcome, treatment, covariates)
# 
# ent_weight <- weightit(reformulate(cov_names, response = "treatment"), data = dataset, 
#                        method = "ebal", moments = 1, verbose = T)
# 
# # Trim to 99th percentile
# ent_weight.trim <- trim(ent_weight, at = .99)
# 
# # Refit with starting weights
# ent_weight <- weightit(reformulate(cov_names, response = "treatment"), data = dataset,
#                        weights = ent_weight.trim$weights, 
#                        method = "ebal", moments = 1, verbose = T)
# 
# # Trim to 99th percentile
# ent_weight.trim <- trim(ent_weight, at = .99)
# 
# gam_ent_fit <-
#   mgcv::bam(reformulate(c("s()", covs), response = "Y"),
#             data = ca_data_trim,
#             weights = ent_weight.trim$weights)
# plot(gam_ent_fit)





