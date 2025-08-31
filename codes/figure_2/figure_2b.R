library(nloptr)
library(ggplot2)
library(RSSthemes)
set_rss_palette("signif_qual")

minimize_variance_bound_nloptr <- function(J, ell, r) {
  
  objective_fn <- function(pi) {
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + 1 / pi[j - ell] + 1 / (1 - pi[j])
    
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * r ** abs(jp - j) * (1 / pi[jp - ell] + 1 / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * r ** abs(jp - j) * (1 / pi[jp - ell] + 1 / (1 - pi[j]) - (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
          }
        }
      }
    }
    (term1 + term2)
  }
  
  ineq_constraints <- function(pi) pi[-length(pi)] - pi[-1]
  
  pi_init <- seq(0.1, 0.9, length.out = J)
  
  opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-25, maxeval = 2e6)
  
  result <- nloptr(x0 = pi_init, eval_f = objective_fn, eval_g_ineq = ineq_constraints,
                   lb = rep(1e-10, J), ub = rep(1 - 1e-10, J), opts = opts)
  
  list(pi = result$solution, objective_value = result$objective,
       status = result$status, iterations = result$iterations)
}

# set up
J <- 10
r <- 0.9

ell_values <- c(0, 1, 2, 3)
results <- data.frame()

for (ell in ell_values) {
  res <- minimize_variance_bound_nloptr(J, ell, r)
  df <- data.frame(
    j = 1:J,
    pi_j = res$pi,
    ell = as.factor(ell)
  )
  results <- rbind(results, df)
}

p <- ggplot(results, aes(x = j, y = pi_j, color = ell, group = ell)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 1:J) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_color_rss_d(palette = "signif_qual", 
                    direction = -1,
                    name = NULL,
                    labels = c("\u2113 = 0", "\u2113 = 1", "\u2113 = 2", "\u2113 = 3")) +
  labs(x = expression("Period "* j), 
       y = expression(pi[j])) +
  theme_bw(base_size = 15) +
  theme(
    text = element_text(family = "serif"),
    
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.title = element_text(size = 14, family = "serif"),
    legend.text = element_text(size = 14, family = "serif"),
    
    axis.title.x = element_text(size = 20, family = "serif"),  # Larger x-axis label
    axis.title.y = element_text(size = 20, angle = 0, vjust = 0.5, family = "serif"),  # Larger y-axis label
    axis.text = element_text(size = 16, color = "black", family = "serif"),

    panel.grid.minor = element_blank()
  )

ggsave("../../figures/figure_2/treated_fraction_ell.pdf", plot = p, 
       width = 4.5, height = 4.5, dpi = 1200, device = cairo_pdf)