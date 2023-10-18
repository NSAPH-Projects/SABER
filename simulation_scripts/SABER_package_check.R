### Script used to test implementation of SABER methodology
# Make sure to enter into the project before doing anything on this regard

# Source helper functions from simulation functions
source("~/nsaph_projects/SABER/simulation_scripts/simulation_functions.R")

set.seed(23)

# Simulate data
data <- data.frame(simulate_data(50, sigma = 10, exp_relationship = "linear"))
w = seq(0, 10, length.out = 10) # values for w that

# Get GPS values
cov_names = c("s", "v")
gps_m <- estimate_gps(cov_mt = data[, cov_names, drop = F],
                      w_all = data$w,
                      sl_lib = c("m_xgboost", "m_ranger", "m_lm"),
                      dnorm_log = TRUE)


params <- list(alpha = 1,
               beta = 1,
               g_sigma = 1,
               tune_app = "all")
outcome_col = "outcome"
treatment_col = "w"
covariates_col = cov_names
nthread = 1
kernel_fn = function(x) exp(-x ^ 2)

# Extract spatial coordinates
spatial_coords = data[, c("x", "y")]

# Now work on devtolls and load all package and dependencies
library(devtools)
load_all()

# Make sure this is able to run through with data setup
results <-
  estimate_cerf_gp_spatial(data.frame(data), w, gps_m, params,
                           outcome_col, treatment_col, covariates_col)

# now go into estimate_cerf_gp_spatial function to test

# In the package, alpha is used to scale w and beta is used to scale
# GPS following the same convention in the method paper.
