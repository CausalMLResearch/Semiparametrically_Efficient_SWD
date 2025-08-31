library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
source("util_empirical_cluster.R")

j_values <- c(8, 10, 12, 14, 16)
true_effect <- 2

process_data <- function(ell) {
  all_mse <- list()
  for (J in j_values) {
    results_dir <- sprintf("../../results/supp_figure_3/cluster_ell_%d_%d", ell, J)
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
  
  combined <- do.call(rbind, all_mse) %>%
    filter(design %in% c("t-opt", "ff", "ffba", "efficient"))
  combined$design[combined$design == "t-opt"] <- "linear"
  combined
}

create_plot <- function(data, plot_title) {
  mse_order <- data %>%
    group_by(design) %>%
    summarise(mean_mse = mean(mse)) %>%
    arrange(desc(mean_mse)) %>%
    pull(design)
  
  data$design <- factor(data$design, levels = mse_order)
  
  design_labels <- c(
    "linear" = expression(bold(Z)[linear]),
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
    scale_y_continuous(limits = c(0.038, 0.072), 
                       breaks = seq(0.04, 0.07, length.out = 4),
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

ell_0_data <- process_data(0)
ell_1_data <- process_data(1)

p1 <- create_plot(ell_0_data, expression(paste(ℓ == 0)))
p2 <- create_plot(ell_1_data, expression(paste(ℓ == 1)))

legend_data <- ell_0_data
mse_order <- legend_data %>%
  group_by(design) %>%
  summarise(mean_mse = mean(mse)) %>%
  arrange(desc(mean_mse)) %>%
  pull(design)
legend_data$design <- factor(legend_data$design, levels = mse_order)

design_labels <- c(
  "linear" = expression(bold(Z)[linear]),
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
  p1, p2, legend,
  ncol = 3,
  widths = c(1, 1, 0.4)
)

ggsave("../../figures/supp_figure_3/mobile_cluster.pdf", final_plot, 
       width = 7, height = 3, device = cairo_pdf)