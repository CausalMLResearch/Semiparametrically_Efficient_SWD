# optimization_utils.R
# utility functions for stepped wedge design optimization

library(nloptr)

setwd("CausalMLResearch/Semiparametrically_Efficient_SWD/shiny_app")

# ============================================================================
# individually randomized trial (IRT) optimization functions
# ============================================================================

# minimize variance bound for gate
minimize_variance_bound_nloptr <- function(J, ell, r, quo) {
  objective_fn <- function(pi) {
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + (quo + 1) / pi[j - ell] + 1 / (1 - pi[j])
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * ((quo + r ** abs(jp - j)) / pi[jp - ell] + r ** abs(jp - j) / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * ((quo + r ** abs(jp - j)) / pi[jp - ell] + r ** abs(jp - j) / (1 - pi[j]) - r ** abs(jp - j) * (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
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

# minimize variance for cte
minimize_variance_cte <- function(J, ell, r, quo) {
  
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
          term3 <- term3 + (quo + 1) / (pi[j - s] - pi_prev) + 1 / (1 - pi[j])
        }
        term4 <- 0
        if ((s + 1) <= (J - 1)) {
          for (j in (s + 1):(J - 1)) {
            for (jp in (j + 1):J) {
              if (j >= jp - s) {
                term4 <- term4 + 2 * (quo + r ** abs(jp - j)) / (1 - pi[j])
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

# minimize variance for both gate and cte
minimize_variance_gate_cte <- function(J, ell, r, quo) {
  
  objective_fn <- function(pi) {
    for(i in 1:(J-1)) {
      if(pi[i+1] - pi[i] < 1e-8) {
        return(1e10)
      }
    }
    
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + (quo + 1) / pi[j - ell] + 1 / (1 - pi[j])
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * ((quo + r ** abs(jp - j)) / pi[jp - ell] + r ** abs(jp - j) / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * ((quo + r ** abs(jp - j)) / pi[jp - ell] + r ** abs(jp - j) / (1 - pi[j]) - r ** abs(jp - j) * (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
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
          term3 <- term3 + (quo + 1) / (pi[j - s] - pi_prev) + 1 / (1 - pi[j])
        }
        term4 <- 0
        if ((s + 1) <= (J - 1)) {
          for (j in (s + 1):(J - 1)) {
            for (jp in (j + 1):J) {
              if (j >= jp - s) {
                term4 <- term4 + 2 * (quo + r ** abs(jp - j)) / (1 - pi[j])
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

# ============================================================================
# cluster randomized trial (CRT) optimization functions
# ============================================================================

# minimize variance bound for crt cross-section design
minimize_variance_crt_cross_section <- function(J, ell, r, sigma_gamma_sq, sigma_epsilon_sq, sigma_upsilon_sq, K) {
  objective_fn <- function(pi) {
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + (sigma_gamma_sq + sigma_epsilon_sq / K + sigma_upsilon_sq) / pi[j - ell] + (sigma_gamma_sq + sigma_epsilon_sq / K) / (1 - pi[j])
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * ((sigma_gamma_sq * r ** abs(jp - j) + sigma_upsilon_sq) / pi[jp - ell] + (sigma_gamma_sq * r ** abs(jp - j)) / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * ((sigma_gamma_sq * r ** abs(jp - j) + sigma_upsilon_sq) / pi[jp - ell] + (sigma_gamma_sq * r ** abs(jp - j)) / (1 - pi[j]) - (sigma_gamma_sq * r ** abs(jp - j)) * (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
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

# minimize variance bound for crt closed-cohort design
minimize_variance_crt_closed_cohort <- function(J, ell, r, sigma_gamma_sq, sigma_epsilon_sq, sigma_upsilon_sq, sigma_phi_sq, K) {
  objective_fn <- function(pi) {
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + (sigma_gamma_sq + sigma_epsilon_sq / K + sigma_phi_sq / K + sigma_upsilon_sq) / pi[j - ell] + (sigma_gamma_sq + sigma_epsilon_sq / K + sigma_phi_sq / K) / (1 - pi[j])
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * ((sigma_gamma_sq * r ** abs(jp - j) + sigma_upsilon_sq + sigma_phi_sq / K) / pi[jp - ell] + (sigma_gamma_sq * r ** abs(jp - j) + sigma_phi_sq / K) / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * ((sigma_gamma_sq * r ** abs(jp - j) + sigma_upsilon_sq + sigma_phi_sq / K) / pi[jp - ell] + (sigma_gamma_sq * r ** abs(jp - j) + sigma_phi_sq / K) / (1 - pi[j]) - (sigma_gamma_sq * r ** abs(jp - j) + sigma_phi_sq / K) * (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
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

# ============================================================================
# wrapper function for irt optimization
# ============================================================================

optimize_irt <- function(J, ell, r, quo, estimand = "gate") {
  # validate inputs
  if(ell >= J - 1) stop("ell must be less than J - 1")
  
  # when ell = 0, all estimands are equivalent
  if(ell == 0) {
    return(minimize_variance_bound_nloptr(J, ell, r, quo))
  }
  
  # choose optimization based on estimand
  if(estimand == "gate") {
    minimize_variance_bound_nloptr(J, ell, r, quo)
  } else if(estimand == "cte") {
    minimize_variance_cte(J, ell, r, quo)
  } else if(estimand == "both") {
    minimize_variance_gate_cte(J, ell, r, quo)
  } else {
    stop("Invalid estimand. Choose 'gate', 'cte', or 'both'")
  }
}

# ============================================================================
# wrapper function for crt optimization
# ============================================================================

optimize_crt <- function(design_type = "cross_section", J, ell, r, sigma_gamma_sq, 
                         sigma_epsilon_sq, sigma_upsilon_sq, K, sigma_phi_sq = NULL) {
  # validate inputs
  if(ell >= J - 1) stop("ell must be less than J - 1")
  
  if(design_type == "cross_section") {
    minimize_variance_crt_cross_section(J, ell, r, sigma_gamma_sq, sigma_epsilon_sq, sigma_upsilon_sq, K)
  } else if(design_type == "closed_cohort") {
    if(is.null(sigma_phi_sq)) stop("sigma_phi_sq is required for closed cohort design")
    minimize_variance_crt_closed_cohort(J, ell, r, sigma_gamma_sq, sigma_epsilon_sq, sigma_upsilon_sq, sigma_phi_sq, K)
  } else {
    stop("Invalid design_type. Choose 'cross_section' or 'closed_cohort'")
  }
}