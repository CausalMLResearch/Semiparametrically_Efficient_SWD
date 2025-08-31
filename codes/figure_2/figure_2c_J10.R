library(nloptr)
library(ggplot2)
library(RSSthemes)
set_rss_palette("signif_qual")

# tau^gate only
minimize_variance_gate <- function(J, ell, r) {
  
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

# tau^gate + tau_0^cte
minimize_variance_gate_cte <- function(J, ell, r) {
  
  objective_fn <- function(pi) {
    for(i in 1:(J-1)) {
      if(pi[i+1] - pi[i] < 1e-8) {
        return(1e10)
      }
    }
    
    term1 <- 0
    for (j in (ell + 1):J) {
      term1 <- term1 + 1 / pi[j - ell] + 1 / (1 - pi[j])
    }
    
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * r^abs(jp - j) * (1 / pi[jp - ell] + 1 / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * r^abs(jp - j) * (1 / pi[jp - ell] + 1 / (1 - pi[j]) - (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
          }
        }
      }
    }
    
    term_cte <- 0
    if (ell > 0) {
      for (s in 0:(ell - 1)) {
        term3 <- 0
        for (j in (s + 1):J) {
          if (j - s - 1 == 0) {
            pi_prev <- 0
          } else {
            pi_prev <- pi[j - s - 1]
          }
          term3 <- term3 + 1 / (pi[j - s] - pi_prev) + 1 / (1 - pi[j])
        }
        
        term4 <- 0
        if ((s + 1) <= (J - 1)) {
          for (j in (s + 1):(J - 1)) {
            for (jp in (j + 1):J) {
              if (j >= jp - s) {
                term4 <- term4 + 2 * r^abs(jp - j) / (1 - pi[j])
              }
            }
          }
        }
        term_cte <- term_cte + (term3 + term4) / ((J - s)^2)
      }
    }
    
    return((term1 + term2) / ((J - ell)^2) + term_cte)
  }
  
  ineq_constraints <- function(pi) pi[-length(pi)] - pi[-1]
  
  pi_init <- seq(0.01, 0.9, length.out = J)
  
  opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-25, maxeval = 2e6)
  
  result <- nloptr(x0 = pi_init, eval_f = objective_fn, eval_g_ineq = ineq_constraints, 
                   lb = rep(1e-10, J), ub = rep(1 - 1e-10, J), opts = opts)
  
  list(pi = result$solution, objective_value = result$objective, 
       status = result$status, iterations = result$iterations)
}

# tau_0^cte only
minimize_variance_cte <- function(J, ell, r) {
  
  objective_fn <- function(pi) {
    for(i in 1:(J-1)) {
      if(pi[i+1] - pi[i] < 1e-8) {
        return(1e10)
      }
    }
    
    term_cte <- 0
    if (ell > 0) {
      for (s in 0:(ell - 1)) {
        term3 <- 0
        for (j in (s + 1):J) {
          if (j - s - 1 == 0) {
            pi_prev <- 0
          } else {
            pi_prev <- pi[j - s - 1]
          }
          term3 <- term3 + 1 / (pi[j - s] - pi_prev) + 1 / (1 - pi[j])
        }
        term4 <- 0
        if ((s + 1) <= (J - 1)) {
          for (j in (s + 1):(J - 1)) {
            for (jp in (j + 1):J) {
              if (j >= jp - s) {
                term4 <- term4 + 2 * r^abs(jp - j) / (1 - pi[j])
              }
            }
          }
        }
        term_cte <- term_cte + (term3 + term4) / ((J - s)^2)
      }
    }
    return(term_cte)
  }
  
  ineq_constraints <- function(pi) pi[-length(pi)] - pi[-1]
  
  pi_init <- seq(0.01, 0.9, length.out = J)
  
  opts <- list(algorithm = "NLOPT_LN_COBYLA", xtol_rel = 1e-25, maxeval = 2e6)
  
  result <- nloptr(x0 = pi_init, eval_f = objective_fn, eval_g_ineq = ineq_constraints, 
                   lb = rep(1e-10, J), ub = rep(1 - 1e-10, J), opts = opts)
  
  list(pi = result$solution, objective_value = result$objective, 
       status = result$status, iterations = result$iterations)
}

# set up
J <- 10
r <- 0.9
ell <- 1

res_gate <- minimize_variance_gate(J, ell, r)
res_gate_cte <- minimize_variance_gate_cte(J, ell, r)
res_cte <- minimize_variance_cte(J, ell, r)

results <- data.frame(
  j = rep(1:J, 3),
  pi_j = c(res_gate$pi, res_gate_cte$pi, res_cte$pi),
  estimator = factor(rep(c("τ^gate", "τ^gate + τ₀^cte", "τ₀^cte"), each = J),
                     levels = c("τ^gate", "τ^gate + τ₀^cte", "τ₀^cte"))
)

p <- ggplot(results, aes(x = j, y = pi_j, color = estimator, group = estimator)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 1:J, limits = c(1, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  scale_color_rss_d(palette = "signif_qual", 
                    direction = -1,
                    name = NULL,
                    labels = c(expression(Var(tau^gate)), 
                               expression(Var(tau^gate) + Var(tau[0]^cte)),
                               expression(Var(tau[0]^cte)))) +
  labs(x = expression("Period "* j), 
       y = expression(pi[j]),
       title = "J = 10, r = 0.9, ℓ = 1") +
  coord_fixed(ratio = 9) + 
  theme_bw(base_size = 15) +
  theme(
    text = element_text(family = "serif"),
    
    plot.title = element_text(size = 18, hjust = 0.5, family = "serif"),
    
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = alpha("white", 0.5), color = "black", size = 0.3),
    legend.key = element_rect(fill = alpha("white", 0.5)),
    legend.title = element_text(size = 14, family = "serif"),
    legend.text = element_text(size = 14, family = "serif"),
    
    axis.title.x = element_text(size = 20, family = "serif"),
    axis.title.y = element_text(size = 20, angle = 0, vjust = 0.5, family = "serif"),
    
    axis.text = element_text(size = 16, color = "black", family = "serif"),
    panel.grid.minor = element_blank()
  )

ggsave("../../figures/figure_2/treated_fraction_vary_estimand.pdf", plot = p, 
       width = 4.5, height = 4.5, dpi = 1200, device = cairo_pdf)