library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggpattern)

results_dir <- "../../results/figure_4/r_med_10"
seeds <- 1:2000
true_effect <- 2
all_results <- list()

for (seed in seeds) {
  filename <- file.path(results_dir, sprintf("sim_results_seed_%04d.rds", seed))
  if (file.exists(filename)) {
    all_results[[seed]] <- readRDS(filename)
  }
}

all_data <- do.call(rbind, lapply(all_results, function(res) {
  if (!is.null(res$combined_summary)) res$combined_summary
}))

aipw_data <- all_data[all_data$method == "aipw" & all_data$design != "ba", ]
designs <- unique(aipw_data$design)
stats <- data.frame()

for (d in designs) {
  subset_data <- aipw_data[aipw_data$design == d, ]
  mse_vals <- (subset_data$estimate - true_effect)^2
  
  mse_mean <- mean(mse_vals)
  mse_se <- sd(mse_vals) / sqrt(length(mse_vals))
  
  row <- data.frame(
    design = d,
    mse = mse_mean,
    mse_se = mse_se,
    mse_lower = mse_mean - mse_se,
    mse_upper = mse_mean + mse_se
  )
  stats <- rbind(stats, row)
}

design_order <- c("t-opt", "ff", "ffba", "efficient_9", "efficient_7", "efficient_5")
stats <- stats[match(design_order, stats$design), ]
stats$x <- 1:6
stats$pattern <- ifelse(stats$design %in% c("efficient_9", "efficient_7"), 
                        "stripe", "none")

labels <- c(
  expression(bold(Z)[linear]),
  expression(bold(Z)[ff]),
  expression(bold(Z)[ffba]),
  expression(bold(Z)[efficient]^{r==0.9}),
  expression(bold(Z)[efficient]^{r==0.7}),
  expression(bold(Z)[efficient]^{r==0.5})
)

design_colors_base <- rev(brewer.pal(6, "Set3"))[1:4]
stats$fill_color <- c(
  design_colors_base[1],  # t-opt
  design_colors_base[2],  # ff
  design_colors_base[3],  # ffba
  "white",
  "white", 
  design_colors_base[4]
)

theme_paper <- theme_bw(base_family = "serif") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 18, vjust = 0.5),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5, family = "serif"),  # title styling
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "none"
  )

p_single <- ggplot(stats, aes(x = factor(x), y = mse)) +
  geom_bar_pattern(
    aes(pattern = pattern, fill = fill_color),
    stat = "identity", 
    color = "black", 
    width = 0.7, 
    size = 0.3,
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_density = 0.005,
    pattern_spacing = 0.04,
    pattern_angle = 45
  ) +
  geom_errorbar(aes(ymin = mse_lower, ymax = mse_upper), width = 0.2, size = 0.4) +
  geom_hline(yintercept = min(stats$mse), linetype = "dashed", color = "black", size = 0.5) +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe")) +
  scale_fill_identity() +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(limits = c(0, 0.062), expand = c(0, 0)) +
  labs(x = "Design", 
       y = "MSE",
       title = "True Model: r = 0.5") +
  theme_paper

ggsave("../../figures/figure_4/vary_r_0.5.pdf", 
       p_single, width = 6, height = 4, 
       device = cairo_pdf)