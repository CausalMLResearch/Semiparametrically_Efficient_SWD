library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpattern)

results_dir <- "../../results/figure_4/multiple_estimands"
seeds <- 1:2000

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

aipw_data <- all_data %>%
  filter(method == "aipw", design != "ba")

calculate_mse_stats <- function(data, true_val) {
  designs <- unique(data$design)
  stats <- data.frame()
  
  for (d in designs) {
    subset_data <- data[data$design == d, ]
    squared_errors <- (subset_data$estimate - true_val)^2
    
    mse_mean <- mean(squared_errors)
    mse_se <- sd(squared_errors) / sqrt(length(squared_errors))
    
    row <- data.frame(
      design = d,
      mse = mse_mean,
      mse_se = mse_se,
      mse_lower = mse_mean - mse_se,
      mse_upper = mse_mean + mse_se
    )
    stats <- rbind(stats, row)
  }
  
  return(stats)
}

gate_data <- aipw_data %>% filter(effect_type == "gate")
gate_stats <- calculate_mse_stats(gate_data, true_val = 2)

cte_data <- aipw_data %>% filter(effect_type == "cte")
cte_stats <- calculate_mse_stats(cte_data, true_val = 1)

combined_stats <- data.frame()
designs <- unique(aipw_data$design)

for (d in designs) {
  design_data <- aipw_data[aipw_data$design == d, ]
  
  gate_subset <- design_data[design_data$effect_type == "gate", ]
  cte_subset <- design_data[design_data$effect_type == "cte", ]
  
  gate_squared_errors <- (gate_subset$estimate - 2)^2  # true effect for gate = 2
  cte_squared_errors <- (cte_subset$estimate - 1)^2   # true effect for cte = 1
  
  all_squared_errors <- c(gate_squared_errors, cte_squared_errors)
  
  combined_mse <- mean(all_squared_errors)
  combined_se <- sd(all_squared_errors) / sqrt(length(all_squared_errors))
  
  row <- data.frame(
    design = d,
    mse = combined_mse,
    mse_se = combined_se,
    mse_lower = combined_mse - combined_se,
    mse_upper = combined_mse + combined_se
  )
  combined_stats <- rbind(combined_stats, row)
}

prepare_plot_data <- function(stats) {
  design_order <- c("t-opt", "ff", "ffba", "efficient_0", "efficient_1")
  stats <- stats[match(design_order, stats$design), ]
  stats$x <- 1:5
  
  stats$pattern <- ifelse(stats$design == "efficient_0", "stripe", "none")
  
  design_colors <- rev(brewer.pal(6, "Set3"))[1:5]
  stats$fill_color <- c(
    design_colors[1],  # t-opt
    design_colors[2],  # ff
    design_colors[3],  # ffba
    "white",           # efficient_0 (striped)
    design_colors[4]   # efficient_1
  )
  
  return(stats)
}

gate_plot_data <- prepare_plot_data(gate_stats)
cte_plot_data <- prepare_plot_data(cte_stats)
combined_plot_data <- prepare_plot_data(combined_stats)

labels <- c(
  expression(bold(Z)[linear]),
  expression(bold(Z)[ff]),
  expression(bold(Z)[ffba]),
  expression(bold(Z)[efficient]^{τ^{gate}}),
  expression(bold(Z)[efficient]^{τ[0]^{cte} + τ^{gate}})
)

theme_paper <- theme_bw(base_family = "serif") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 18, vjust = 0.5),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 20, hjust = 0.5, family = "serif"),  # title styling
    plot.margin = margin(10, 10, 10, 18),
    legend.position = "none"
  )

create_mse_plot <- function(plot_data, y_label, plot_title, y_max = NULL) {
  if (is.null(y_max)) {
    y_max <- max(plot_data$mse_upper) * 1.1
  }
  
  p <- ggplot(plot_data, aes(x = factor(x), y = mse)) +
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
    geom_errorbar(aes(ymin = mse_lower, ymax = mse_upper), 
                  width = 0.2, size = 0.4) +
    geom_hline(yintercept = min(plot_data$mse), 
               linetype = "dashed", color = "black", size = 0.5) +
    scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe")) +
    scale_fill_identity() +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0), n.break = 8) +
    labs(x = "Design", y = y_label, title = plot_title) +  # add title
    theme_paper
  
  return(p)
}

p_gate <- create_mse_plot(gate_plot_data, 
                          expression("MSE for"~τ^{gate}),
                          expression("Estimand:"~τ^{gate}),  # title
                          y_max = 0.142)

p_cte <- create_mse_plot(cte_plot_data, 
                         expression("MSE for"~τ[0]^{cte}),
                         expression("Estimand:"~τ[0]^{cte}),  # title
                         y_max = 0.142)

p_combined <- create_mse_plot(combined_plot_data, 
                              expression("MSE for"~τ[0]^{cte}~"and"~τ^{gate}),
                              expression("Estimands:"~τ[0]^{cte}~"and"~τ^{gate}),  # title
                              y_max = 0.142)

ggsave("../../figures/figure_4/multiple_gate.pdf", 
       p_gate, width = 6, height = 4, 
       device = cairo_pdf)
ggsave("../../figures/figure_4/multiple_cte.pdf", 
       p_cte, width = 6, height = 4, 
       device = cairo_pdf)
ggsave("../../figures/figure_4/multiple_both.pdf", 
       p_combined, width = 6, height = 4, 
       device = cairo_pdf)