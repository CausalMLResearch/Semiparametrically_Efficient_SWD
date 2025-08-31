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
  
  pi_init <- seq(0.4, 0.6, length.out = J)
  
  opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-25, maxeval = 2e6)
  
  result <- nloptr(x0 = pi_init, eval_f = objective_fn, eval_g_ineq = ineq_constraints,
                   lb = rep(1e-10, J), ub = rep(1 - 1e-10, J), opts = opts)
  
  list(pi = result$solution, objective_value = result$objective,
       status = result$status, iterations = result$iterations)
}

# set up
J <- 2
ell <- 0

r_seq <- seq(0, 1, length.out = 100)
results <- data.frame()

for (r_val in r_seq) {
  opt <- minimize_variance_bound_nloptr(J, ell, r_val)
  df <- data.frame(
    r = r_val,
    pi_1 = opt$pi[1],
    pi_2 = opt$pi[2]
  )
  results <- rbind(results, df)
}

results_long <- data.frame(
  r = rep(results$r, 2),
  pi_value = c(results$pi_1, results$pi_2),
  pi_type = factor(rep(c("pi1", "pi2"), each = nrow(results)), 
                   levels = c("pi1", "pi2"))
)

p <- ggplot(results_long, aes(x = r, y = pi_value, color = pi_type, group = pi_type)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_color_rss_d(palette = "signif_qual", 
                    name = NULL,
                    direction = -1,
                    labels = c(expression(pi[1]), expression(pi[2]))) +
  labs(x = expression("Decay Rate "* r), 
       y = "Treated Fraction") +
  theme_bw(base_size = 15) +
  theme(

    text = element_text(family = "serif"),
    
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "black", size = 0.3),
    legend.title = element_text(size = 14, family = "serif"),  # legend title size 14
    legend.text = element_text(size = 14, family = "serif"),   # legend text size 14
    
    axis.title.x = element_text(size = 20, family = "serif"),  # larger x-axis label
    axis.title.y = element_text(size = 20, vjust = 0.5, family = "serif"),  # larger y-axis label
    axis.text = element_text(size = 16, color = "black", family = "serif"),
    
    panel.grid.minor = element_blank()
  )

ggsave("../../figures/figure_2/J_2.pdf", plot = p, 
       width = 4.5, height = 4.5, dpi = 1200, device = cairo_pdf)