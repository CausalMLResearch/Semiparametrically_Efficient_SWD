library(dplyr)
library(ggplot2)
library(patchwork)

process_params <- function(cohort) {
  results_dir <- sprintf("../../results/supp_figure_2/%s_para", cohort)
  n_sims <- 2000
  all_results <- list()
  
  for (seed in 1:n_sims) {
    filename <- file.path(results_dir, sprintf("corr_params_seed_%04d.rds", seed))
    if (file.exists(filename)) {
      all_results[[seed]] <- readRDS(filename)
    }
  }
  
  combined_results <- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
      seed = x$seed,
      r = x$r,
      sigma_epsilon = x$sigma_epsilon,
      sigma_residual = x$sigma_residual,
      status = x$status,
      stringsAsFactors = FALSE
    )
  }))
  
  combined_results
}

patient_data <- process_params("patient")
caregiver_data <- process_params("caregiver")

combined_data <- rbind(
  patient_data %>% mutate(type = "Patient"),
  caregiver_data %>% mutate(type = "Caregiver")
) %>%
  filter(status == "success") %>%
  mutate(
    sigma_epsilon_sq = sigma_epsilon^2,
    sigma_residual_sq = sigma_residual^2
  )

combined_data$type <- factor(combined_data$type, levels = c("Patient", "Caregiver"))

p_r <- ggplot(combined_data, aes(x = type, y = r, fill = type)) +
  geom_boxplot(width = 0.4, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(0.75, 1), expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "", y = expression(italic(r))) +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3, color = "gray90"),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.position = "none",
    aspect.ratio = 1
  )

p_sigma_gamma <- ggplot(combined_data, aes(x = type, y = sigma_epsilon_sq, fill = type)) +
  geom_boxplot(width = 0.4, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(0.4, 1.4), expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "", y = expression(sigma[gamma]^2)) +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3, color = "gray90"),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.position = "none",
    aspect.ratio = 1
  )

p_sigma_epsilon <- ggplot(combined_data, aes(x = type, y = sigma_residual_sq, fill = type)) +
  geom_boxplot(width = 0.4, position = position_dodge(width = 0.6)) +
  scale_y_continuous(limits = c(0., 0.5), expand = expansion(mult = c(0.02, 0.02))) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "", y = expression(sigma[epsilon]^2)) +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3, color = "gray90"),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    legend.position = "none",
    aspect.ratio = 1
  )

combined_plot <- p_r | p_sigma_gamma | p_sigma_epsilon

ggsave("../../figures/supp_figure_2/mobile_parameters.pdf", combined_plot, 
       width = 12, height = 4, device = cairo_pdf)