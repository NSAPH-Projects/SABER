# Load necessary libraries
library(tidyverse)
library(zipcodeR)
library(GPCERF)
library(SuperLearner)
library(WeightIt)
library(data.table)

# Clear workspace
rm(list = ls())

# Name directories
proj_dir <- "/n/dominici_nsaph_l3/Lab/projects/GPCERF_spatial/"
data_folder <- paste0(proj_dir, "/GPCERF_spatial/cork_data/")
data_prep_tag <- "california"
data_dir <- paste0(data_folder, data_prep_tag, "/")

# Define covariates
covs <- c("year", "mean_bmi", "smoke_rate", "hispanic", "pct_blk", 
          "medhouseholdincome", "medianhousevalue", "poverty", "education",
          "popdensity", "pct_owner_occ", "summer_tmmx", "winter_tmmx", "summer_rmax",
          "winter_rmax")

# Are you creating new dataset? If not preload
new_data <- F

if (new_data) {
  # Create directory to store data and some results (for now, consider moving)
  dir.create(data_dir)
  
  # Load Medicare data from NEJM source 
  load("/n/dominici_nsaph_l3/Lab/projects/analytic/erc_strata/aggregate_data.RData")
  
  # Convert to tibble (tidyverse equivalent of data.table)
  data <- as_tibble(aggregate_data)
  
  # Subset to year 2000 and NORTHEAST region
  data <- data %>% filter(region == "WEST")
  
  # Include only the covariates and the exposure, removing stratification by age/sex
  data <- 
    data %>% 
    group_by(across(all_of(c("zip", "pm25", covs)))) %>%
    summarize(dead = sum(dead), time_count = sum(time_count), .groups = "drop") %>%
    select(any_of(c("zip", "dead", "time_count", "pm25", covs)))
  
  # Map zip codes to states 
  zip_code_states <- reverse_zipcode(unique(data$zip))
  zip_code_states <- zip_code_states %>% select(zipcode, state, lat, lng)
  
  # Merge with original data
  data_state <- 
    data %>%
    rename(zipcode = zip) %>% 
    left_join(zip_code_states)
  
  # Filter to CA state only
  ca_data <- data_state %>% filter(state == "CA")
  
  # # Include only the covariates and the exposure, removing stratification by age/sex
  # ma_data <- 
  #   ma_data %>% 
  #   group_by(across(all_of(c("zipcode", "year", "pm25", covs)))) %>%
  #   summarize(dead = sum(dead), time_count = sum(time_count), .groups = "drop") %>%
  #   select(any_of(c("zipcode", "dead", "time_count", "pm25", covs)))
  
  # Create mortality rate and rename pm25 as w
  ca_data <- 
    ca_data %>% 
    mutate(mort_rate = dead / time_count, 
           w = pm25)
  
  # Replace zero value with half of minimum observed
  min_value_obs <- min(ca_data$mort_rate[ca_data$mort_rate > 0])
  ca_data$mort_rate[ca_data$mort_rate == 0] <- min_value_obs / 2
  
  # Make outcome the log of the mortality rate 
  ca_data <- ca_data %>% mutate(Y = log(mort_rate))
  
  # Subset to what you need to create CERF
  ca_data_cerf <- ca_data %>% select(any_of(c("Y", "w", covs)))
  
  # Check for any missing columns
  missing_columns <- 
    ca_data_cerf %>%
    summarise(across(everything(), ~any(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "column", values_to = "has_missing") %>%
    filter(has_missing)
  
  missing_columns
  
  # Save ma_data_cerf for quicker use (can run from here later)
  saveRDS(ca_data_cerf, paste0(data_dir, "ca_data_cerf.RDS"))
} else {
  ca_data_cerf <- readRDS(paste0(data_dir, "ca_data_cerf.RDS"))
}

# not stratified...does that seem reasonable?
# We need to specify a kernel as well for this use case (Guassian, Matern, etc) 
# and I'm not sure where that is estimated
ca_data_cerf <- readRDS(paste0(data_dir, "ca_data_cerf.RDS"))

# Subset to first 16 years? 
ca_data_cerf <- ca_data_cerf %>% filter(year %in% 2006:2007)

# Define exposure values to estimate
q1 <- quantile(ca_data_cerf$w, 0.05)
q2 <- quantile(ca_data_cerf$w, 0.95)

# Trim data for fitting entropy balancing 
ca_data_trim <- 
  ca_data_cerf %>% 
  filter(between(w, q1, q2))

# Now going to fit entropy weighting to this to compare to what I get from CERF
ent_weight <- weightit(reformulate(covs, response = "w"), data = ca_data_trim, 
                       method = "ebal", moments = 1, verbose = T)

# Trim to 99th percentile
ent_weight.trim <- trim(ent_weight, at = .99)

# # Truncate the 99th percentile
# upper_bound = quantile(ent_weight$weights, 0.995)
# ent_weight$weights[ent_weight$weights > upper_bound] <- upper_bound
gam_ent_fit <-
  mgcv::bam(reformulate(c("s(w)", covs), response = "Y"),
            data = ca_data_trim,
            weights = ent_weight.trim$weights)
plot(gam_ent_fit)

# gam_ent_fit2 <- mgcv::bam(reformulate(c("s(w, k = 6)", covs), response = "Y"),
#                          data = ca_data_trim,
#                          weights = ent_weight.trim$weights)
# plot(gam_ent_fit2)

plot(mgcv::bam(reformulate(c("s(w)", covs), response = "Y"),
          data = ca_data_trim))

data_prediction <- 
  rbindlist(lapply(seq(q1, q2, length.out = 20), function(pot_exp) {
    
    # Get potential data if all had same potential exposure
    potential_data <- 
      select(ca_data_trim, all_of(covs)) %>% 
      mutate(w = pot_exp)
    
    # Fit each model and take mean for ERC
    potential_outcome <- 
      potential_data %>% 
      mutate(
        gam_ent = predict(gam_ent_fit, potential_data, type = "response"),
      ) %>% 
      dplyr::select(w, gam_ent) %>% 
      summarize_all(mean) %>%
      data.table()
    return(potential_outcome)
  }))

# Specify model types
model_types <- c("gam_ent")

# Create mortality estimates by model type
mortality_data <- 
  data_prediction %>% 
  pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
  arrange(w) %>% 
  mutate(prediction = exp(prediction))

# Create plot of ERC by model
mortality_plot <- 
  mortality_data %>%
  ggplot(aes(x = w, y = prediction, color = model, linetype = model)) +
  geom_line() +
  labs(x = "Exposure concentration", y = "Log mortality rate") + 
  theme_bw()


#### Now estimate using nearest 
# Define number of cores
n_core <- 12

# Define SuperLearner wrappers for xgboost and ranger with thread number specified
m_xgboost <- function(nthread = n_core, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}

m_ranger <- function(num.threads = n_core, ...) {
  SuperLearner::SL.ranger(num.threads = num.threads, ...)
}

new_data <- F
# Load GPS estimates (or run again if new dataset)
if (new_data) {
  # Estimation of Generalized Propensity Score using SuperLearner
  gps_m <- estimate_gps(
    cov_mt = ca_data_trim[,-c(1,2)], 
    w_all = ca_data_trim$w,
    sl_lib = c("m_xgboost", "m_ranger"),
    dnorm_log = TRUE
  )
  
  saveRDS(gps_m, paste0(data_dir, "gps_m.RDS"))
} else {
  gps_m <- readRDS(paste0(data_dir, "gps_m.RDS"))
}

# First fit nearest neighbor to get estimates 
# Define exposure values to estimate
# q1 <- quantile(ca_data_cerf$w, 0.05)
# q2 <- quantile(ca_data_cerf$w, 0.95)
w_all <- seq(q1, q2, length.out = 20)

# Define parameters list for search across grid space
# Not sure how to define number of neighbors and block size 
params_lst <- list(alpha = 10 ^ seq(-2, 2, length.out = 10),
                   beta = 10 ^ seq(-2, 2, length.out = 10),
                   g_sigma = c(0.1, 1, 10),
                   tune_app = "all",
                   n_neighbor = 50,
                   block_size = 1e3)

params_lst2 <- list(alpha = 10 ^ seq(0, 2, length.out = 10),
                    beta = 10 ^ seq(0, 2, length.out = 10),
                    g_sigma = c(0.1, 1, 10),
                    tune_app = "all",
                    n_neighbor = 50,
                    block_size = 1e3)

params_lst3 <- list(alpha = 6,
                    beta = 12,
                    g_sigma = 0.1,
                    tune_app = "all",
                    n_neighbor = 50,
                    block_size = 1e3)

# Run nearest neighbor on this data frame (this takes around ~10 minutes)
cerf_nngp_obj <- estimate_cerf_nngp(data.frame(ca_data_trim),
                                    w_all,
                                    gps_m,
                                    params = params_lst,
                                    nthread = 12)

cerf_nngp_obj2 <- estimate_cerf_nngp(data.frame(ca_data_trim),
                                     w_all,
                                     gps_m,
                                     params = params_lst2,
                                     nthread = 12)


cerf_nngp_obj3 <- estimate_cerf_nngp(data.frame(ca_data_trim),
                                     w_all,
                                     gps_m,
                                     params = params_lst3,
                                     nthread = 10)


summary(cerf_nngp_obj)
summary(cerf_nngp_obj3)

plot(cerf_nngp_obj)
plot(cerf_nngp_obj3)


# Save results
saveRDS(cerf_nngp_obj, paste0(data_dir, "cerf_nngp_obj"))

cerf_nngp_obj <- readRDS(paste0(data_dir, "cerf_nngp_obj"))


# Estimation of Exposure Response Function using full data, not nearest neighbor
params_lst <- list(alpha = 10 ^ seq(-2, 2, length.out = 10),
                   beta = 10 ^ seq(-2, 2, length.out = 10),
                   g_sigma = c(0.1, 1, 10),
                   tune_app = "all")

params_lst2 <- list(alpha = 6,
                    beta = 12,
                    g_sigma = 0.1,
                    tune_app = "all")


cerf_gp_obj <- 
  estimate_cerf_gp(data.frame(ca_data_trim),
                   w_all,
                   gps_m,
                   params_lst,
                   n_core)

cerf_gp_obj2 <- 
  estimate_cerf_gp(data.frame(ca_data_trim),
                   w_all,
                   gps_m,
                   params_lst2,
                   nthread = 10)



# Save results
saveRDS(cerf_gp_obj, paste0(data_dir, "cerf_gp_obj.RDS"))
cerf_gp_obj <- readRDS(paste0(data_dir, "cerf_gp_obj.RDS"))

# Summarize and plot the results
summary(cerf_gp_obj)
plot(cerf_gp_obj)

message("done")

# #### Here is example from GPCERF vignette which is helpful to reference ------------------
# set.seed(781)
# sim_data <- generate_synthetic_data(sample_size = 500, gps_spec = 1)
# 
# n_core <- 1
# 
# m_xgboost <- function(nthread = n_core, ...) {
#   SuperLearner::SL.xgboost(nthread = nthread, ...)
# }
# 
# m_ranger <- function(num.threads = n_core, ...){
#   SuperLearner::SL.ranger(num.threads = num.threads, ...)
# }
# 
# # Estimate GPS function
# gps_m <- estimate_gps(cov_mt = sim_data[,-(1:2)],
#                       w_all = sim_data$treat,
#                       sl_lib = c("m_xgboost", "m_ranger"),
#                       dnorm_log = TRUE)
# 
# # exposure values
# q1 <- stats::quantile(sim_data$treat, 0.05)
# q2 <- stats::quantile(sim_data$treat, 0.95)
# 
# w_all <- seq(q1, q2, 1)
# 
# params_lst <- list(alpha = 10 ^ seq(-2, 2, length.out = 10),
#                    beta = 10 ^ seq(-2, 2, length.out = 10),
#                    g_sigma = c(0.1, 1, 10),
#                    tune_app = "all")
# 
# cerf_gp_obj <- estimate_cerf_gp(sim_data,
#                                 w_all,
#                                 gps_m,
#                                 params = params_lst,
#                                 nthread = n_core)
# summary(cerf_gp_obj)
# plot(cerf_gp_obj)
# 
# 
# # Now quantify the strength and scale of the spatial correlatio
# # First look at the spatial scale
# # Load the necessary libraries
# library(spdep)
# library(gstat)
# 
# # For Moran statistic you can use binary weighting scheme, might be most helpful way to do this
# 
# 
# # I'm assuming your data is a spatial object (like a SpatialPointsDataFrame) and that you have a neighbors list defined.
# # This is just a placeholder, replace it with your actual data and neighbors list
# nb <- nb
# 
# # Loop through each column in your data
# for(i in 1:ncol(data)) {
# 
#   # Compute Moran's I
#   moran_result <- moran.test(ma_data[, i], nb, style="B")
# 
#   print(paste("Moran's I for variable", colnames(data)[i], "is", moran_result$estimate["Moran I"]))
# 
#   # Compute the variogram
#   variogram_result <- variogram(data[, i]~1, data)
# 
#   # Plot the variogram
#   plot(variogram_result, main = paste("Variogram for variable", colnames(data)[i]))
# 
#   # Extract the range (spatial scale) from the fitted model
#   vgm_model <- fit.variogram(variogram_result, model = vgm(1, "Sph", 300, 1))
#   print(paste("Spatial scale for variable", colnames(data)[i], "is", vgm_model$range))
# 
# }
# 
