library(dplyr)
library(lme4)
library(tidyverse)
library(MASS)
library(xtable)
library(glmmTMB)
library(nloptr)

minimize_variance_bound_nloptr <- function(J, ell, r, sigma_epsilon, sigma_residual) {
  
  objective_fn <- function(pi) {
    term1 <- 0
    for (j in (ell + 1):J) term1 <- term1 + 1 / pi[j - ell] + 1 / (1 - pi[j])
    
    term2 <- 0
    if ((ell + 1) <= (J - 1)) {
      for (j in (ell + 1):(J - 1)) {
        for (jp in (j + 1):J) {
          if (j >= jp - ell) {
            term2 <- term2 + 2 * (sigma_epsilon ** 2 * r ** abs(jp - j)) / (sigma_epsilon ** 2 + sigma_residual ** 2) * (1 / pi[jp - ell] + 1 / (1 - pi[j]))
          }
        }
        if ((j + ell + 1) <= J) {
          for (jp in (j + ell + 1):J) {
            term2 <- term2 + 2 * (sigma_epsilon ** 2 * r ** abs(jp - j)) / (sigma_epsilon ** 2 + sigma_residual ** 2) * (1 / pi[jp - ell] + 1 / (1 - pi[j]) - (pi[jp - ell] - pi[j]) / ((1 - pi[j]) * pi[jp - ell]))
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

generate_design_matrix <- function(I, J,
                                   design_type = "t-opt",
                                   efficient_matrix = NULL) {
  
  # design matrix
  X <- matrix(0, nrow = I, ncol = J)
  
  if (design_type == "t-opt") {
    # t-opt
    for (j in 1:J) {
      treated_fraction <- (2*j - 1) / (2*J)
      n_treated <- round(treated_fraction * I)
      
      if (n_treated > 0) {
        X[1:n_treated, j] <- 1
      }
    }
    
  } else if (design_type == "ba") {
    # ba
    X[, (floor(J / 2) + 1):J] <- 1
    X[1, ] <- 1
    
  } else if (design_type == "ff") {
    # ff
    half_I <- floor(I / 2)
    X[(half_I + 1):I, ] <- 1
    
  } else if (design_type == "ffba") {
    # ffba
    half_J <- floor(J / 2)
    half_I <- floor(I / 2)
    X[(half_I + 1):I, (half_J + 1):J] <- 1
    
  } else if (grepl("^efficient", design_type)) {
    
    for (j in 1:J) {
      n_treated <- round(efficient_matrix[j] * I)
      
      if (n_treated > 0) {
        X[1:n_treated, j] <- 1
      }
    }
  }
  
  return(X)
}

simulate_synthetic_panel <- function(df,
                                     I,
                                     J_pre,
                                     J_exp,
                                     design_type = "t-opt",
                                     delta = c(0, 2),
                                     max_exposure_time = 1,
                                     seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # get total periods and available units
  max_time <- max(df$time)
  all_units <- unique(df$id_unit)
  
  # check if enough units available
  if (length(all_units) < I) {
    stop("not enough units in data for specified I")
  }
  
  # randomly select I units and reassign consecutive ids
  selected_units <- sample(all_units, I, replace = FALSE)
  
  # filter data to selected units and reassign ids from 1 to I
  df <- df %>%
    filter(id_unit %in% selected_units) %>%
    mutate(id_unit = match(id_unit, selected_units))
  
  # check if enough periods available
  if (max_time < J_pre + J_exp) {
    stop("not enough periods in data for J_pre + J_exp")
  }
  
  # randomly select starting time point
  start_time <- sample(1:(max_time - J_pre - J_exp + 1), 1)
  end_time <- start_time + J_pre + J_exp - 1
  
  # extract window of data
  df_window <- df %>%
    filter(time >= start_time & time <= end_time) %>%
    mutate(
      relative_time = time - start_time + 1,
      period = ifelse(relative_time <= J_pre, "historical", "synthetic")
    )
  
  # split into historical and synthetic
  df_historical <- df_window %>%
    filter(period == "historical") %>%
    mutate(time = as.factor(relative_time)) %>%
    mutate(id_unit= as.factor(id_unit))
  
  df_synthetic <- df_window %>%
    filter(period == "synthetic") %>%
    mutate(time = as.factor(relative_time - J_pre))
  
  efficient_matrix <- NULL
  if (design_type == "efficient") {
    fit <- glmmTMB(out_cont_1 ~ time + ar1(time + 0 | id_unit), data = df_historical)
    r <- unname(attr(VarCorr(fit)$cond$id_unit, "correlation")[1, 2])
    sigma_epsilon <- unname(attr(VarCorr(fit)$cond$id_unit, "stddev")[1])
    sigma_residual <- unname(attr(VarCorr(fit)$cond, "sc"))
    ell <- max_exposure_time - 1
    opt_result <- minimize_variance_bound_nloptr(
      J = J_exp, 
      ell = ell, 
      r = r, 
      sigma_epsilon = sigma_epsilon, 
      sigma_residual = sigma_residual
    )
    efficient_matrix <- opt_result$pi
  }
  
  # generate design matrix for synthetic periods
  X <- generate_design_matrix(I, J_exp, design_type, efficient_matrix)
  X_long <- as.vector(t(X))
  
  # add treatment to synthetic data
  df_synthetic <- df_synthetic %>%
    arrange(id_unit, time) %>%
    mutate(trt = X_long)
  
  # calculate treatment exposure time
  df_adoption <- df_synthetic %>%
    filter(trt == 1) %>%
    group_by(id_unit) %>%
    summarise(first_time = min(as.numeric(as.character(time))), .groups = 'drop')
  
  df_synthetic <- df_synthetic %>%
    left_join(df_adoption, by = "id_unit") %>%
    mutate(
      time_on_trt_raw = ifelse(is.na(first_time) | trt == 0, 0, 
                               as.numeric(as.character(time)) - first_time + 1),
      time_on_trt = pmin(time_on_trt_raw, max_exposure_time)
    ) %>%
    mutate(out_cont_1 = out_cont_1 + delta[time_on_trt + 1]) %>%
    mutate(time_on_trt = as.factor(time_on_trt)) %>%
    dplyr::select(-time_on_trt_raw, -first_time)
  
  df_synthetic$id_unit <- as.factor(df_synthetic$id_unit)
  df_synthetic$time <- as.factor(df_synthetic$time)
  df_synthetic$time_on_trt <- as.factor(df_synthetic$time_on_trt)
  df_synthetic$out_cont_1 <- as.numeric(df_synthetic$out_cont_1)
  
  # set attributes for synthetic data
  attr(df_synthetic, "design_type") <- design_type
  attr(df_synthetic, "design_matrix") <- X
  attr(df_synthetic, "parameters") <- list(
    I = I,
    J = J_exp,
    J_pre = J_pre,
    delta = delta,
    max_exposure_time = max_exposure_time,
    start_time = start_time,
    end_time = end_time
  )
  
  return(df_synthetic)
}

outcome_regression_gate <- function(data) {
  
  # extract max exposure time from data attributes
  max_exposure_time <- attr(data, "parameters")$max_exposure_time
  
  # fit the outcome regression model
  model <- lm(out_cont_1 ~ time_on_trt + time, data = data)
  
  data_treated <- data
  data_treated$time_on_trt <- factor(max_exposure_time, levels = levels(data$time_on_trt))
  data_control <- data
  data_control$time_on_trt <- factor(0, levels = levels(data$time_on_trt))
  
  pred_treated <- predict(model, newdata = data_treated)
  pred_control <- predict(model, newdata = data_control)
  
  gate_estimate <- mean((pred_treated - pred_control))
  
  result <- list(
    estimate = gate_estimate,
    model = model,
    method = "outcome_regression",
    max_exposure_time = max_exposure_time,
    all_coef = coef(model)
  )
  
  class(result) <- c("gate_estimate", "list")
  return(result)
}

ipw_gate <- function(data) {
  
  max_exposure_time <- attr(data, "parameters")$max_exposure_time
  ell <- max_exposure_time - 1
  J <- attr(data, "parameters")$J
  I <- attr(data, "parameters")$I
  design_matrix <- attr(data, "design_matrix")
  
  prop_scores <- colMeans(design_matrix)
  
  data <- data %>%
    mutate(time_numeric = as.numeric(as.character(time)),
           prop_score = prop_scores[time_numeric])
  
  ipw_sum_treated <- 0
  ipw_sum_control <- 0
  denom_treated <- 0
  denom_control <- 0
  
  for (j in (ell+1):J) {
    data_j <- data %>% filter(time_numeric == j)
    
    treated_units <- data_j %>% 
      filter(time_on_trt == max_exposure_time)
    
    if (nrow(treated_units) > 0) {
      pi_j_minus_ell <- prop_scores[j - ell]
      ipw_sum_treated <- ipw_sum_treated + sum(treated_units$out_cont_1 / pi_j_minus_ell)
      denom_treated <- denom_treated + nrow(treated_units) / pi_j_minus_ell
    }
  }
  
  for (j in 1:J) {
    data_j <- data %>% filter(time_numeric == j)
    
    control_units <- data_j %>% 
      filter(trt == 0)
    
    if (nrow(control_units) > 0) {
      pi_j <- prop_scores[j]
      ipw_sum_control <- ipw_sum_control + sum(control_units$out_cont_1 / (1 - pi_j))
      denom_control <- denom_control + nrow(control_units) / (1 - pi_j)
    }
  }
  
  gate_estimate <- ipw_sum_treated / denom_treated - ipw_sum_control / denom_control
  
  result <- list(
    estimate = gate_estimate,
    method = "ipw",
    max_exposure_time = max_exposure_time,
    ell = ell,
    prop_scores = prop_scores
  )
  
  class(result) <- c("gate_estimate", "list")
  return(result)
}

aipw_gate <- function(data) {
  
  max_exposure_time <- attr(data, "parameters")$max_exposure_time
  ell <- max_exposure_time - 1
  J <- attr(data, "parameters")$J
  I <- attr(data, "parameters")$I
  design_matrix <- attr(data, "design_matrix")
  
  # estimated propensity scores
  prop_scores <- colMeans(design_matrix)
  
  data <- data %>%
    mutate(time_numeric = as.numeric(as.character(time)),
           prop_score = prop_scores[time_numeric])
  
  # outcome regression model
  outcome_model <- lm(out_cont_1 ~ time_on_trt + time, data = data)
  
  # for treated: predict with time_on_trt = max_exposure_time
  data_treated_scenario <- data %>%
    mutate(time_on_trt = factor(max_exposure_time, levels = levels(data$time_on_trt)))
  data$pred_treated <- predict(outcome_model, newdata = data_treated_scenario)
  
  # for control: predict with time_on_trt = 0
  data_control_scenario <- data %>%
    mutate(time_on_trt = factor(0, levels = levels(data$time_on_trt)))
  data$pred_control <- predict(outcome_model, newdata = data_control_scenario)
  
  # AIPW components
  aipw_sum_treated <- 0
  aipw_sum_control <- 0
  denom_treated <- 0
  denom_control <- 0
  
  # treated component
  for (j in (ell+1):J) {
    data_j <- data %>% filter(time_numeric == j)
    
    # units with exposure time = max_exposure_time at time j
    treated_units <- data_j %>% 
      filter(time_on_trt == max_exposure_time)
    
    if (nrow(treated_units) > 0) {
      pi_j_minus_ell <- prop_scores[j - ell]
      
      # weighted residual
      aipw_sum_treated <- aipw_sum_treated + sum((treated_units$out_cont_1 - treated_units$pred_treated) / pi_j_minus_ell)
      denom_treated <- denom_treated + nrow(treated_units) / pi_j_minus_ell
    }
  }
  
  # control component
  for (j in 1:J) {
    data_j <- data %>% filter(time_numeric == j)
    
    control_units <- data_j %>% 
      filter(trt == 0)
    
    if (nrow(control_units) > 0) {
      pi_j <- prop_scores[j]
      
      # weighted residual
      aipw_sum_control <- aipw_sum_control + sum((control_units$out_cont_1 - control_units$pred_control) / (1 - pi_j))
      denom_control <- denom_control + nrow(control_units) / (1 - pi_j)
    }
  }
  
  gate_estimate <- mean((data$pred_treated - data$pred_control)) + aipw_sum_treated / denom_treated - aipw_sum_control / denom_control
  
  result <- list(
    estimate = gate_estimate,
    outcome_model = outcome_model,
    method = "aipw",
    max_exposure_time = max_exposure_time,
    ell = ell,
    prop_scores = prop_scores
  )
  
  class(result) <- c("gate_estimate", "list")
  return(result)
}

load_all_design_results <- function(results_dir = "simulation_results_all_designs", seeds = 1:2000) {
  
  all_results <- list()
  
  for (seed in seeds) {
    filename <- file.path(results_dir, sprintf("sim_results_seed_%04d.rds", seed))
    
    if (file.exists(filename)) {
      all_results[[as.character(seed)]] <- readRDS(filename)
    } else {
      warning(sprintf("File not found for seed %d", seed))
    }
  }
  
  return(all_results)
}

summarize_all_designs <- function(results_dir = "simulation_results_all_designs", seeds = 1:2000) {
  
  all_results <- load_all_design_results(results_dir, seeds)
  
  all_summaries <- do.call(rbind, lapply(names(all_results), function(seed) {
    res <- all_results[[seed]]
    if (!is.null(res$combined_summary)) {
      cbind(seed = as.integer(seed), res$combined_summary)
    }
  }))
  
  # summary statistics by design and method
  summary_stats <- all_summaries %>%
    group_by(design, method) %>%
    summarise(
      mean_estimate = mean(estimate, na.rm = TRUE),
      bias = mean(estimate - true_effect, na.rm = TRUE),
      variance = mean((mean(estimate, na.rm = TRUE) - estimate)^2, na.rm = TRUE),
      mse = mean((estimate - true_effect)^2, na.rm = TRUE),
      mae = mean(abs(estimate - true_effect), na.rm = TRUE),
      n_successful = sum(!is.na(estimate)),
      .groups = 'drop'
    )
  
  return(list(
    summary_stats = summary_stats,
    raw_summaries = all_summaries
  ))
}

generate_table <- function(summary_results,
                           output_file = NULL,
                           digits = 3,
                           caption = "Simulation Results",
                           label = "tab:results") {
  
  if (is.list(summary_results)) {
    summary_stats <- summary_results$summary_stats
  } else {
    summary_stats <- summary_results
  }
  
  formatted_data <- summary_stats %>%
    mutate(
      Design = design,
      Method = case_when(
        method == "outcome_regression" ~ "OR",
        method == "ipw" ~ "IPW",
        method == "aipw" ~ "AIPW",
        TRUE ~ method
      ),
      Estimate = mean_estimate,
      Bias = bias,
      Variance = variance,
      MSE = mse,
      MAE = mae,
      N = n_successful
    ) %>%
    dplyr::select(Design, Method, Estimate, Bias, Variance, MSE, MAE, N)
  
  xt <- xtable(formatted_data,
               caption = caption,
               label = label,
               digits = c(0, 0, digits, digits, digits, digits, digits, digits, 0))
  
  display(xt) <- c("s", "s", "e", "e", "e", "e", "e", "e", "d")
  
  if (!is.null(output_file)) {
    print(xt,
          file = output_file,
          include.rownames = FALSE,
          booktabs = TRUE,
          math.style.exponents = TRUE,
          hline.after = c(-1, 0, 
                          which(diff(as.numeric(factor(formatted_data$Design))) != 0),
                          nrow(formatted_data)))
    
    cat("Table saved to:", output_file, "\n")
  }
  
  return(xt)
}