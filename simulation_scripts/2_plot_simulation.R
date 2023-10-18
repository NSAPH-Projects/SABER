# Source helper functions from simulation functions
library(tidyverse)
source("~/nsaph_projects/SABER/simulation_scripts/simulation_functions.R", echo=TRUE)

# Now run the model simulations
model_tag <- "increase_nngp"
results_dir <- paste0("~/nsaph_projects/SABER/model_runs/", model_tag, "/")
plot_dir <- paste0(results_dir, "/plots/")
dir.create(plot_dir)

# Create correlation table 
results <- readRDS(paste0(results_dir, "results_finished.RDS"))
correlation_table <- map_dfr(results, function(df) df$correlation_table)

correlation_summary <- 
  correlation_table %>% 
  pivot_longer(cols = c("pre_cor", "post_cor"), names_to = "type", values_to = "value") %>% 
  group_by(covariate, method, type) %>% 
  summarize(mean = mean(value),
            lower = quantile(value, 0.05),
            upper = quantile(value, 0.95)) %>% 
  ungroup()

# Now plot the correlation plot 
gg_correlation <- 
  ggplot(correlation_summary, aes(x = covariate, y = mean, ymin = lower, ymax = upper, color = type, group = type)) +
  geom_hline(aes(yintercept = 0.1), alpha = 0.5, color = "black", linetype = "dashed") + 
  geom_errorbar(width = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7)) + 
  scale_color_discrete("", breaks = c("pre_cor", "post_cor"), labels=c("Pre-weighting", "Post-weighting")) + 
  theme_bw(base_size = 14) + 
  theme(text = element_text(size=14),
        axis.text.x = element_text(size = 12)) + 
  coord_flip() + 
  facet_wrap(~ method) + 
  labs(y = "Absolute correlation", x = "Covariate")

ggsave(paste0(plot_dir, "correlation.pdf"), gg_correlation, width = 8, height = 8)

# Now work on metrics -----------------------------------------
metrics <- map_dfr(1:50, function(i){
  results[[i]]$metrics %>% mutate(sim = i)
})

metrics %>% 
  #filter(model != "linear_ent") %>% 
  group_by(model) %>%
  summarize(mean_mse = mean(mse),
            lower_mse = quantile(mse, 0.05),
            upper_mse = quantile(mse, 0.95)) %>% 
  ggplot(aes(x = model, y = mean_mse, ymin = lower_mse, ymax = upper_mse)) + 
  geom_point() + 
  coord_cartesian(ylim = c(0, 100)) + 
  geom_errorbar(width = 0.5) + 
  theme_bw() + 
  labs(y = "MSE")

# Plot absolute bias how you report it in figures 
metrics %>% 
  group_by(w, model) %>%
  summarize(abs_bias = abs(mean(bias))) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  summarize(abs_bias = mean(abs_bias)) %>% 
  ungroup() %>% 
  ggplot(aes(x = model, y = abs_bias)) + 
  geom_point() + 
  coord_cartesian(ylim = c(0, 10)) + 
  #geom_errorbar(width = 0.5) + 
  theme_bw() + 
  labs(y = "Absolute Bias")

# Now read in results after you have fit the model
simulations <- readRDS(paste0(results_dir, "results_finished.RDS"))

# Define common grid to work over
common_exposure_grid <- seq(0, 10, length.out = 100)

# Get data predictions adjusted and pivot long by model type
simulations_long <- 
  simulations %>%
  map(~ .[["data_prediction"]] %>% 
        pivot_longer(cols = -w, names_to = "model", values_to = "response"))

# Use the interpolate_on_grid function
interpolated_data <- interpolate_on_grid(simulations_long, common_exposure_grid)

# Summarize interpolated data 
summary_data <- 
  interpolated_data %>%
  group_by(model, w) %>%
  summarise(
    mean = mean(prediction, na.rm = TRUE),
    lower = quantile(prediction, 0.025, na.rm = TRUE),
    upper = quantile(prediction, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

# Extract the true Exposure-response curve to plot
true_rel <-
  filter(summary_data, model == "true_relationship") %>%
  select(w, true_relationship = mean)

# now make plot of data using "s"
gg_erf <- 
  summary_data %>% 
  filter(model != "true_relationship") %>% 
  ggplot() +
  geom_ribbon(aes(x = w, ymin = lower, ymax = upper, fill = model), color = NA, alpha = 0.3) +
  geom_line(aes(x = w, y = mean, color = model)) +
  geom_line(data = true_rel, aes(x = w, y = true_relationship), linetype = "dashed", color = "black") +
  labs(x = "Exposure", y = "Outcome", title = "Exposure-Response with no hidden variables") +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 35)) +
  #coord_cartesian(xlim = c(0, 10), ylim = c(-7.5, 12.5)) +
  theme_bw() +
  theme(legend.position = "") + 
  facet_wrap(~ model)

gg_erf

ggsave(paste0(plot_dir, "erf_transparent.pdf"), gg_erf, width = 10, height = 8)

# Now read in results after you have fit the model ------------------------------------------
# Define common grid to work over
common_exposure_grid <- seq(0, 10, length.out = 100)
simulations_hidden <- readRDS(paste0(results_dir, "results_hidden_finished.RDS"))

# simulations_hidden <- map(1:97, function(i){
#   readRDS(paste0(results_dir, "results_hidden_", i, ".RDS"))
# })


# Get data predictions adjusted and pivot long by model type
simulations_long <- 
  simulations_hidden %>%
  map(~ .[["data_prediction"]] %>% 
        pivot_longer(cols = -w, names_to = "model", values_to = "response"))

# Use the interpolate_on_grid function
interpolated_data <- interpolate_on_grid(simulations_long, common_exposure_grid)

# Summarize interpolated data 
summary_data_hidden <- 
  interpolated_data %>%
  group_by(model, w) %>%
  summarise(
    mean = mean(prediction, na.rm = TRUE),
    lower = quantile(prediction, 0.025, na.rm = TRUE),
    upper = quantile(prediction, 0.975, na.rm = TRUE)
  ) %>%
  ungroup()

# Extract the true Exposure-response curve to plot
true_rel <-
  filter(summary_data_hidden, model == "true_relationship") %>%
  select(w, true_relationship = mean)

# now make plot of data using "s"
gg_erf_hidden <- 
  summary_data_hidden %>% 
  filter(model != "true_relationship") %>% 
  ggplot() +
  geom_ribbon(aes(x = w, ymin = lower, ymax = upper, fill = model), color = NA, alpha = 0.3) +
  geom_line(aes(x = w, y = mean, color = model)) +
  geom_line(data = true_rel, aes(x = w, y = true_relationship), linetype = "dashed", color = "black") +
  labs(x = "Exposure", y = "Outcome", title = "Exposure-Response with hidden variables") +
  #coord_cartesian(xlim = c(0, 10), ylim = c(-2.5, 20)) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 25)) +
  theme_bw() +
  theme(legend.position = "") + 
  facet_wrap(~ model)

gg_erf_hidden
ggsave(paste0(plot_dir, "erf_hidden.pdf"), gg_erf_hidden, width = 10, height = 8)

# Now work on metrics -----------------------------------------
metrics_hidden <- map_dfr(1:50, function(i){
  simulations_hidden[[i]]$metrics %>% mutate(sim = i)
})

metrics_hidden %>% 
  filter(model != "linear_ent") %>% 
  group_by(model) %>%
  summarize(mean_mse = mean(mse),
            lower_mse = quantile(mse, 0.10),
            upper_mse = quantile(mse, 0.9)) %>% 
  ggplot(aes(x = model, y = mean_mse, ymin = lower_mse, ymax = upper_mse)) + 
  geom_point() + 
  coord_cartesian(ylim = c(0, 20)) + 
  geom_errorbar(width = 0.5) + 
  theme_bw() + 
  labs(y = "MSE")

# Plot absolute bias how you report it in figures 
metrics_hidden %>% 
  group_by(w, model) %>%
  summarize(abs_bias = abs(mean(bias))) %>% 
  ungroup() %>% 
  group_by(model) %>% 
  summarize(abs_bias = mean(abs_bias)) %>% 
  ungroup() %>% 
  ggplot(aes(x = model, y = abs_bias)) + 
  geom_point() + 
  coord_cartesian(ylim = c(0, 20)) + 
  #geom_errorbar(width = 0.5) + 
  theme_bw() + 
  labs(y = "Absolute Bias")

# Now I want to plot the hidden and not hidden on the same plot so 
# that I can make sure it's doing what I want it to do
gg_bias_compare <- 
  rbind(metrics %>% mutate(setting = "transparent"),
        metrics_hidden %>% mutate(setting = "hidden")) %>%
  mutate(setting = factor(setting, levels = c("transparent", "hidden"))) %>% 
  # filter(model != "linear_ent") %>% 
  group_by(w, model, setting) %>%
  summarize(abs_bias = abs(mean(bias))) %>% 
  ungroup() %>% 
  group_by(model, setting) %>% 
  summarize(abs_bias = mean(abs_bias)) %>% 
  ungroup() %>% 
  ggplot(aes(x = model, y = abs_bias, color = setting)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 10)) + 
  theme_bw() + 
  labs(y = "Absolute bias") +
  scale_color_manual(values = c("transparent" = "deepskyblue", "hidden" = "darkorange4"))

gg_bias_compare

ggsave(paste0(plot_dir, "bias_compare.pdf"), gg_bias_compare, width = 8, height = 8)

gg_mse_compare <- 
  rbind(metrics %>% mutate(setting = "transparent"),
        metrics_hidden %>% mutate(setting = "hidden")) %>%
  mutate(setting = factor(setting, levels = c("transparent", "hidden"))) %>% 
  # filter(model != "linear_ent") %>% 
  group_by(model, setting) %>%
  summarize(mean_mse = mean(mse),
            lower_mse = quantile(mse, 0.05),
            upper_mse = quantile(mse, 0.95)) %>% 
  ggplot(aes(x = model, y = mean_mse, ymin = lower_mse, ymax = upper_mse, color = setting)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  geom_errorbar(width = 0.5, position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 50)) + 
  theme_bw() + 
  labs(y = "MSE") +
  scale_color_manual(values = c("transparent" = "deepskyblue", "hidden" = "darkorange4"))

gg_mse_compare
ggsave(paste0(plot_dir, "bias_compare.pdf"), gg_mse_compare, width = 8, height = 8)



