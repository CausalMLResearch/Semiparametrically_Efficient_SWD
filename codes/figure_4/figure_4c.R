library(dplyr)
library(ggplot2)
library(RColorBrewer)
source("util_IRT_sim.R")

j_values <- c(8, 10, 12, 14, 16)
true_effect <- 2
all_mse <- list()

for (J in j_values) {
  results_dir <- sprintf("../../results/figure_4/r_high_%d", J)
  results <- load_all_design_results(results_dir, seeds = 1:2000)
  aipw_data <- do.call(rbind, lapply(results, function(res) {
    if (!is.null(res$combined_summary)) {
      res$combined_summary %>% filter(method == "aipw")
    }
  }))
  
  if (J == 10) {
    aipw_data <- aipw_data %>%
      filter(design %in% c("ba", "ff", "ffba", "t-opt", "efficient_9")) %>%
      mutate(design = ifelse(design == "efficient_9", "efficient", design))
  }
  
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

combined_data <- do.call(rbind, all_mse) %>%
  filter(design %in% c("t-opt", "ff", "ffba", "efficient"))

mse_order <- combined_data %>%
  group_by(design) %>%
  summarise(mean_mse = mean(mse)) %>%
  arrange(desc(mean_mse)) %>%
  pull(design)

combined_data$design <- factor(combined_data$design, levels = mse_order)

design_labels <- c(
  "t-opt" = expression(bold(Z)[linear]),
  "ff" = expression(bold(Z)[ff]),
  "ffba" = expression(bold(Z)[ffba]),
  "efficient" = expression(bold(Z)[efficient])
)

design_colors <- setNames(rev(brewer.pal(6, "Set3"))[1:length(mse_order)], mse_order)

p_mse_line <- ggplot(combined_data, aes(x = J, y = mse, color = design, fill = design)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linetype = 1) +
  geom_line(size = 1) +
  geom_point(size = 2, stroke = 0.2) +
  scale_color_manual(values = design_colors, labels = design_labels) +
  scale_fill_manual(values = design_colors, labels = design_labels) +
  scale_x_continuous(breaks = j_values) +
  scale_y_continuous(limits = c(0.06, 0.124), 
                     breaks = seq(0.06, 0.12, length.out = 6),
                     expand = c(0, 0)) +
  labs(x = "J", y = "MSE") +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3, color = "gray90"),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    legend.text.align = 0
  )

ggsave("../../figures/figure_4/vary_J.pdf", p_mse_line, width = 6, height = 4, dpi = 300, device = cairo_pdf)