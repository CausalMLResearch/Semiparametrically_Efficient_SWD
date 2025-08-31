library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
source("util_IRT_sim.R")

j_values <- c(8, 10, 12, 14, 16)
true_effect <- 2

process_data <- function(cohort, ell) {
  all_mse <- list()
  for (J in j_values) {
    results_dir <- sprintf("../../results/figure_5/%s_ell_%d_%d", cohort, ell, J)
    results <- load_all_design_results(results_dir, seeds = 1:2000)
    aipw_data <- do.call(rbind, lapply(results, function(res) {
      if (!is.null(res$combined_summary)) {
        res$combined_summary %>% filter(method == "aipw")
      }
    }))
    
    mse_results <- aipw_data %>%
      group_by(design) %>%
      summarise(
        mse = mean((estimate - true_effect)^2),
        sd_mse = sd((estimate - true_effect)^2),
        n = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        J = J,
        se = sd_mse / sqrt(n),
        t_crit = qt(0.975, df = n - 1),
        ci_lower = mse - t_crit * se,
        ci_upper = mse + t_crit * se
      )
    all_mse[[as.character(J)]] <- mse_results
  }
  do.call(rbind, all_mse) %>%
    filter(design %in% c("t-opt", "ff", "ffba", "efficient"))
}

create_plot <- function(data, plot_title) {
  mse_order <- data %>%
    group_by(design) %>%
    summarise(mean_mse = mean(mse)) %>%
    arrange(desc(mean_mse)) %>%
    pull(design)
  
  data$design <- factor(data$design, levels = mse_order)
  
  design_labels <- c(
    "t-opt" = expression(bold(Z)["t-opt"]),
    "ff" = expression(bold(Z)[ff]),
    "ffba" = expression(bold(Z)[ffba]),
    "efficient" = expression(bold(Z)[efficient])
  )
  
  design_colors <- setNames(rev(brewer.pal(6, "Set3"))[1:length(mse_order)], mse_order)
  
  ggplot(data, aes(x = J, y = mse, color = design, fill = design)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linetype = 1) +
    geom_line(size = 1) +
    geom_point(size = 2, stroke = 0.2) +
    scale_color_manual(values = design_colors, labels = design_labels) +
    scale_fill_manual(values = design_colors, labels = design_labels) +
    scale_x_continuous(breaks = j_values) +
    scale_y_continuous(limits = c(0.038, 0.112), 
                       breaks = seq(0.04, 0.11, length.out = 8),
                       expand = c(0, 0)) +
    labs(x = "J", y = "MSE", title = plot_title) +
    theme_bw(base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3, color = "gray90"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      plot.title = element_text(size = 16, hjust = 0.5),
      plot.margin = margin(5, 5, 5, 5),
      legend.position = "none",
      aspect.ratio = 1
    )
}

p_0_data <- process_data("patient", 0)
p_1_data <- process_data("patient", 1)
c_0_data <- process_data("caregiver", 0)
c_1_data <- process_data("caregiver", 1)

p1 <- create_plot(p_0_data, expression(paste("Patient Cohort: ", ℓ == 0)))
p2 <- create_plot(p_1_data, expression(paste("Patient Cohort: ", ℓ == 1)))
p3 <- create_plot(c_0_data, expression(paste("Caregiver Cohort: ", ℓ == 0)))
p4 <- create_plot(c_1_data, expression(paste("Caregiver Cohort: ", ℓ == 1)))

legend_data <- p_0_data
mse_order <- legend_data %>%
  group_by(design) %>%
  summarise(mean_mse = mean(mse)) %>%
  arrange(desc(mean_mse)) %>%
  pull(design)
legend_data$design <- factor(legend_data$design, levels = mse_order)

design_labels <- c(
  "t-opt" = expression(bold(Z)[linear]),
  "ff" = expression(bold(Z)[ff]),
  "ffba" = expression(bold(Z)[ffba]),
  "efficient" = expression(bold(Z)[efficient])
)
design_colors <- setNames(rev(brewer.pal(6, "Set3"))[1:length(mse_order)], mse_order)

p_legend <- ggplot(legend_data, aes(x = J, y = mse, color = design, fill = design)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2) +
  geom_line() +
  scale_color_manual(values = design_colors, labels = design_labels) +
  scale_fill_manual(values = design_colors, labels = design_labels) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 18, family = "serif"),
    legend.text.align = 0
  )

legend <- cowplot::get_legend(p_legend)

final_plot <- grid.arrange(
  p1, p2, p3, p4, legend,
  ncol = 5,
  widths = c(1, 1, 1, 1, 0.4)
)

ggsave("../../figures/figure_5/mobile.pdf", final_plot, 
       width = 13, height = 3, device = cairo_pdf)