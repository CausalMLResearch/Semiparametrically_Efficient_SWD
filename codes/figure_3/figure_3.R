library(ggplot2)
library(patchwork)
library(RColorBrewer)

# setup
J <- 10
periods <- 1:J

designs <- list(
  ff = rep(0.5, J),
  ba = c(rep(0, J/2), rep(1, J/2)),
  ffba = c(rep(0, J/2), rep(0.5, J/2)),
  `linear` = (2*(1:J) - 1)/(2*J),
  efficient = c(0.1638731, 0.2765203, 0.3531312, 0.4159630, 0.4725966, 
                0.5274034, 0.5840370, 0.6468688, 0.7234797, 0.8361270)
)

dark2_color <- brewer.pal(3, "Dark2")[3]

make_plot <- function(pi, title) {
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
    scale_x_continuous(breaks = 1:10, expand = c(0, 0), limits = c(0.5, 10.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = bquote(bold(Z)[.(title)]), x = NULL, y = NULL) +
    theme_bw(base_family = "serif") +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 18, family = "serif"),
          axis.text.x = element_text(size = 12, family = "serif"),
          panel.grid.minor = element_blank(),
          plot.margin = margin(5, 5, 5, 5))
}

plots <- Map(make_plot, designs, names(designs))

combined_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + 
  plot_layout(nrow = 1)

print(combined_plot)

ggsave("../../figures/figure_3/designs.pdf", 
       plot = combined_plot,
       width = 7, 
       height = 2,
       dpi = 1200,
       device = "pdf")