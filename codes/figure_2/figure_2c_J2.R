library(ggplot2)
library(patchwork)
library(RColorBrewer)

J <- 2
periods <- 1:J

pi1_cte <- (sqrt(5) - 1)/4
pi2_cte <- (sqrt(5) + 3)/8

pi1_combined <- 0.3976559
pi2_combined <- 0.5837904

dark2_color <- brewer.pal(3, "Dark2")[3]

make_plot <- function(pi, title_expr) {
  x_step <- numeric()
  y_control <- numeric()
  y_treated <- numeric()
  
  for(i in 1:J) {
    x_step <- c(x_step, i - 0.5, i + 0.5)
    y_control <- c(y_control, 1 - pi[i], 1 - pi[i])
    y_treated <- c(y_treated, pi[i], pi[i])
  }
  
  df <- data.frame(
    x = rep(x_step, 2),
    y = c(y_control, y_treated),
    group = factor(rep(c("Control", "Treated"), each = length(x_step)),
                   levels = c("Control", "Treated"))
  )
  
  ggplot(df, aes(x = x, y = y, fill = group)) +
    geom_area(position = "stack") +
    scale_fill_manual(values = c("Control" = "white", "Treated" = dark2_color)) +
    geom_hline(yintercept = 0.5, color = "grey", linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = 1.5, color = "grey", linewidth = 1, linetype = "dashed") +
    scale_x_continuous(breaks = 1:2, expand = c(0, 0), limits = c(0.5, 2.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = title_expr, x = NULL, y = NULL) +
    theme_bw(base_family = "serif") +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          panel.grid.minor = element_blank(),
          plot.margin = margin(10, 10, 10, 10),
          aspect.ratio = 1)
}

plot1 <- make_plot(
  c(0.5, 0.5),
  expression(bold(Var(tau^gate)))
)

plot2 <- make_plot(
  c(pi1_cte, pi2_cte),
  expression(bold(Var(tau[0]^cte)))
)

plot3 <- make_plot(
  c(pi1_combined, pi2_combined),
  expression(bold(Var(tau[0]^cte) + Var(tau^gate)))
)

combined_plot <- plot1 + plot2 + plot3 + 
  plot_layout(nrow = 1) +
  plot_annotation(
    title = "Two-Period Experiments with ℓ = 1",
    theme = theme(
      plot.title = element_text(size = 20, hjust = 0.5, family = "serif"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave("../../figures/figure_2/multiple_estimands.pdf", 
       plot = combined_plot,
       width = 6,
       height = 3,
       dpi = 1200,
       device = "pdf")