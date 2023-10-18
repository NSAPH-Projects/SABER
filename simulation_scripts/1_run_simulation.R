# Function to perform small spatial simulation to test GPCERF package
library("GPCERF", lib.loc = "/n/home_fasse/mcork/R/ifxrstudio/RELEASE_3_16")

# Source helper functions from simulation functions
source("~/nsaph_projects/SABER/simulation_scripts/simulation_functions.R")

set.seed(23)

# Main function for simulation
simulate_models <- function(n = 1000, sigma = 5, n_core = 12, 
                            exp_relationship = "linear", spatial_hidden = F,
                            interaction = T) {
  dat <- simulate_data(n, sigma, exp_relationship = exp_relationship, interaction = interaction)
  model_fits <- fit_models(dat, n_core = n_core, 
                           spatial_hidden, w_values = seq(0, 10, length.out = 100))
  metrics <- calculate_metrics(dat, lm_fit = model_fits$lm_fit, 
                               gam_fit = model_fits$gam_fit,
                               lm_ent_weights = model_fits$lm_ent_weights, 
                               gam_ent_weights = model_fits$gam_ent_weights,
                               gpcerf = model_fits$gpcerf,
                               w_values = seq(0, 10, length.out = 100),
                               interaction = interaction)
  
  return(list(
    metrics = metrics$data_metrics,
    data_prediction_adjusted = metrics$prediction_adjusted,
    data_prediction = metrics$data_prediction,
    correlation_table = model_fits$correlation_table
  ))
}

# time_taken <- system.time({
#   sim_result <- simulate_models(n = 500, sigma = 10, n_core = 1, 
#                                 spatial_hidden = F, interaction = F)
# })
 
# print(time_taken)

# # # # This is for testing it out before running 100 simulations
# sim_result <- simulate_models(n = 200, sigma = 10, n_core = 1, spatial_hidden = F)
# sim_result2 <- simulate_models(n = 1000, sigma = 5, n_core = 12, spatial_hidden = T)
 
# model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "gpcerf")
# 
# # # plot all data together, centering at w = 5
# sim_result$data_prediction %>%
#   pivot_longer(c(all_of(model_types), "true_relationship"), names_to = "model", values_to = "prediction") %>%
#   arrange(w) %>%
#   ggplot() +
#   geom_line(aes(x = w, y = prediction, color = model)) +
#   labs(x = "PM2.5", y = "Mortality")

# sim_result2$metrics %>% 
#   group_by(model) %>% 
#   summarize(mse = mean(mse))

# Now run the model simulations
model_tag <- "no_interaction"
results_dir <- paste0("~/nsaph_projects/SABER/model_runs/", model_tag, "/")
dir.create(results_dir)

# Run over several times to get through simulation
results <- lapply(1:50, function(x) {
  sim_result <- simulate_models(n = 1000, sigma = 10, n_core = 1, 
                                exp_relationship = "linear", spatial_hidden = F,
                                interaction = F)
  print(paste("Finished simulation", x))
  # Save for now so if it fails you can still recover
  saveRDS(sim_result, paste0(results_dir, "results_", x, ".RDS"))
  return(sim_result)
})

saveRDS(results, paste0(results_dir, "results_finished.RDS"))

# run with hiding the spatial covariate
results_hidden <- lapply(1:50, function(x) {
  sim_result <- simulate_models(n = 1000, sigma = 10, n_core = 1, 
                                exp_relationship = "linear", spatial_hidden = T,
                                interaction = F)
  print(paste("Finished hidden simulation", x))
  # Save for now so if it fails you can still recover
  saveRDS(sim_result, paste0(results_dir, "results_hidden_", x, ".RDS"))
  return(sim_result)
})

saveRDS(results_hidden, paste0(results_dir, "results_hidden_finished.RDS"))

# Now go on to plotting simulations file for this next portion

# # Now read in results after you have fit the model
# simulations <- readRDS(paste0(results_dir, "results_finished.RDS"))
# 
# # Define common grid to work over
# common_exposure_grid <- seq(0, 10, length.out = 100)
# 
# # Get data predictions adjusted and pivot long by model type
# simulations_long <- 
#   simulations %>%
#   map(~ .[["data_prediction_adjusted"]] %>% 
#         pivot_longer(cols = -w, names_to = "model", values_to = "response"))
# 
# # Create interpolated dataset for each exposure value
# interpolated_data <- 
#   simulations_long %>%
#   map_dfr(function(df) {
#     df %>%
#       group_by(model) %>%
#       do(data.frame(
#         w = common_exposure_grid,
#         model = unique(.$model),
#         prediction = approx(.$w, .$response, xout = common_exposure_grid)$y
#       ))
#   })
# 
# # Summarize interpolated data 
# summary_data <- 
#   interpolated_data %>%
#   group_by(model, w) %>%
#   summarise(
#     mean = mean(prediction, na.rm = TRUE),
#     lower = quantile(prediction, 0.025, na.rm = TRUE),
#     upper = quantile(prediction, 0.975, na.rm = TRUE)
#   ) %>%
#   ungroup()
# 
# 
# data_prediction <- bind_rows(lapply(results, function(df) df$data_prediction))
# correlation_table <- bind_rows(lapply(results, function(df) df$correlation_table))
# 
# # plot all data together, centering at w = 5
# model_types <- c("linear_model", "gam_model", "linear_ent", "gam_ent", "gpcerf")
# 
# data_adjusted_summary <-
#   data_prediction_adjusted %>%
#   pivot_longer(c(all_of(model_types), "true_relationship"), names_to = "model", values_to = "prediction") %>%
#   group_by(w, model) %>%
#   summarize(mean = mean(prediction),
#             upper = quantile(prediction, 0.975),
#             lower = quantile(prediction, 0.025)) %>%
#   arrange(w) %>%
#   ungroup()
# 
# data_adjusted_summary %>%
#   ggplot() +
#   geom_line(aes(x = w, y = mean, color = model)) +
#   geom_ribbon(aes(x = w, ymin = lower, ymax = upper, fill = model), color = NA, alpha = 0.1) +
#   coord_equal() +
#   labs(x = "PM2.5", y = "Mortality")
# 
# # Extract the true Exposure-response curve
# true_rel <-
#   filter(summary_data, model == "true_relationship") %>%
#   select(w, true_relationship = mean)
# 
# 
# summary_data %>% 
#   filter(model != "true_relationship") %>% 
#   ggplot() +
#   geom_ribbon(aes(x = w, ymin = lower, ymax = upper, fill = model), color = NA, alpha = 0.3) +
#   geom_line(aes(x = w, y = mean, color = model)) +
#   geom_line(data = true_rel, aes(x = w, y = true_relationship), linetype = "dashed", color = "black") +
#   labs(x = "Exposure", y = "Outcome") +
#   coord_cartesian(xlim = c(0, 10), ylim = c(-10, 10)) +
#   theme_bw() +
#   facet_wrap(~ model)
# 
# 
# # Plot all vs the truth, facet by model
# data_prediction_adjusted %>%
#   select(-true_relationship) %>%
#   pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
#   group_by(w, model) %>%
#   summarize(mean = mean(prediction),
#             upper = quantile(prediction, 0.975),
#             lower = quantile(prediction, 0.025)) %>%
#   arrange(w) %>%
#   ungroup() %>%
#   ggplot() +
#   geom_ribbon(aes(x = w, ymin = lower, ymax = upper, fill = model), color = NA, alpha = 0.1) +
#   geom_line(aes(x = w, y = mean, color = model)) +
#   geom_line(data = true_rel, aes(x = w, y = true_relationship), linetype = "dashed", color = "black") +
#   labs(x = "Exposure", y = "Outcome") +
#   coord_cartesian(xlim = c(0, 10), ylim = c(-10, 10)) +
#   theme_bw() +
#   facet_wrap(~ model)
# 
# data_prediction_adjusted %>%
#   select(-true_relationship) %>%
#   pivot_longer(c(all_of(model_types)), names_to = "model", values_to = "prediction") %>%
#   group_by(w, model) %>%
#   summarize(mean = mean(prediction),
#             upper = quantile(prediction, 0.975),
#             lower = quantile(prediction, 0.025)) %>%
#   arrange(w) %>%
#   ungroup() %>%
#   ggplot() +
#   geom_ribbon(aes(x = w, ymin = lower, ymax = upper, fill = model), color = NA, alpha = 0.1) +
#   geom_line(aes(x = w, y = mean, color = model)) +
#   geom_line(data = true_rel, aes(x = w, y = true_relationship, color = "True Relationship", linetype = "True Relationship")) +
#   labs(x = "Exposure", y = "Outcome") +
#   scale_color_manual(values = c("True Relationship" = "black")) +
#   scale_linetype_manual(values = c("True Relationship" = "dashed")) +
#   theme_minimal() +
#   facet_wrap(~ model)
# 
