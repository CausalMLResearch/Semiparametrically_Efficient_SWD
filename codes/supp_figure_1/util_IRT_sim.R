library(dplyr)
library(lme4)
library(tidyverse)
library(MASS)
library(xtable)
library(glmmTMB)

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

simulate_panel_data <- function(I = 100, 
                                J = 10,
                                design_type = NULL,
                                efficient_matrix = NULL,
                                delta = c(0, 2),
                                max_exposure_time = 1,
                                sigma_epsilon = 1,
                                r = 0.9,
                                seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # design matrix
  X <- generate_design_matrix(I, J, design_type, efficient_matrix)
  X_long <- as.vector(t(X))
  
  # residual
  structure_decay <- matrix(NA, ncol = J, nrow = J)
  for (j in 1:(J - 1)) {
    for (k in (j + 1):J) {
      structure_decay[j, k] <- sigma_epsilon * sigma_epsilon * r^abs(j - k)
      structure_decay[k, j] <- structure_decay[j, k]
    }
  }
  
  for (j in 1:J){
    structure_decay[j, j] <- sigma_epsilon * sigma_epsilon
  }
  
  dtd <- MASS::mvrnorm(n = I, mu = rep(0, J), Sigma = structure_decay)
  
  Y_ij <- c()
  for (i in 1:I) {
    for (j in 1:J) {
      temp <- (i - 1) * J + j
      Y_ij[temp] <- j/J + dtd[i, j]
    }
  }
  
  df <- data.frame(
    id_unit = as.factor(rep(1:I, each = J)),
    time = as.factor(rep(1:J, times = I)),
    out_cont_1 = Y_ij,
    trt = X_long
  )
  
  df_adoption <- df %>%
    filter(trt == 1) %>%
    group_by(id_unit) %>%
    summarise(first_time = min(as.numeric(as.character(time))), .groups = 'drop')
  
  df <- df %>%
    left_join(df_adoption, by = "id_unit") %>%
    mutate(
      # time_on_trt
      time_on_trt_raw = ifelse(is.na(first_time) | trt == 0, 0, 
                               as.numeric(as.character(time)) - first_time + 1),
      # time_on_trt is less than or equal to max_exposure_time
      time_on_trt = pmin(time_on_trt_raw, max_exposure_time)
    ) %>%
    mutate(out_cont_1 = out_cont_1 + delta[time_on_trt + 1]) %>%
    mutate(time_on_trt = as.factor(time_on_trt)) %>%
    dplyr::select(-time_on_trt_raw)
  
  attr(df, "design_type") <- design_type
  attr(df, "design_matrix") <- X
  attr(df, "parameters") <- list(
    I = I,
    J = J,
    delta = delta,
    max_exposure_time = max_exposure_time,
    sigma_epsilon = sigma_epsilon,
    r = r
  )
  
  return(df)
}

outcome_regression_gate <- function(data) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
  I                 <- params$I
  
  model <- lm(out_cont_1 ~ time_on_trt + time, data = data)
  
  data_trt_scenario <- data %>%
    mutate(time_on_trt = factor(max_exposure_time, levels = levels(data$time_on_trt)))
  data_ctrl_scenario <- data %>%
    mutate(time_on_trt = factor(0, levels = levels(data$time_on_trt)))
  data$mu_1 <- predict(model, newdata = data_trt_scenario)
  data$mu_0 <- predict(model, newdata = data_ctrl_scenario)
  
  data <- data %>%
    mutate(time_numeric = as.numeric(as.character(time)))
  
  tau_j <- numeric(J - ell)
  for (idx in seq_len(J - ell)) {
    j <- ell + idx
    data_j <- data %>% filter(time_numeric == j)
    tau_j[idx] <- mean(data_j$mu_1 - data_j$mu_0)   # = (1/I) * sum_i
  }
  
  gate_estimate <- mean(tau_j)
  
  result <- list(
    estimate          = gate_estimate,
    model             = model,
    method            = "outcome_regression",
    max_exposure_time = max_exposure_time,
    ell               = ell,
    all_coef          = coef(model)
  )
  
  class(result) <- c("gate_estimate", "list")
  return(result)
}

ipw_gate <- function(data) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
  I                 <- params$I
  design_matrix     <- attr(data, "design_matrix")
  
  prop_scores <- colMeans(design_matrix)
  
  data <- data %>%
    mutate(time_numeric = as.numeric(as.character(time)))
  
  tau_j <- numeric(J - ell)
  for (idx in seq_len(J - ell)) {
    j <- ell + idx
    data_j <- data %>% filter(time_numeric == j)
    
    pi_trt  <- prop_scores[j - ell]
    pi_ctrl <- 1 - prop_scores[j]
    
    I_trt  <- as.numeric(data_j$time_on_trt == max_exposure_time)
    I_ctrl <- as.numeric(data_j$trt == 0)
    
    trt_term  <- ifelse(I_trt  == 1, data_j$out_cont_1 / pi_trt,  0)
    ctrl_term <- ifelse(I_ctrl == 1, data_j$out_cont_1 / pi_ctrl, 0)
    
    summand <- trt_term - ctrl_term
    
    tau_j[idx] <- mean(summand)
  }
  
  gate_estimate <- mean(tau_j)
  
  result <- list(
    estimate          = gate_estimate,
    method            = "ipw",
    max_exposure_time = max_exposure_time,
    ell               = ell,
    prop_scores       = prop_scores
  )
  
  class(result) <- c("gate_estimate", "list")
  return(result)
}

aipw_gate <- function(data) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
  I                 <- params$I
  design_matrix     <- attr(data, "design_matrix")
  
  prop_scores <- colMeans(design_matrix)
  
  outcome_model <- lm(out_cont_1 ~ time_on_trt + time, data = data)
  
  data_trt_scenario <- data %>%
    mutate(time_on_trt = factor(max_exposure_time, levels = levels(data$time_on_trt)))
  data_ctrl_scenario <- data %>%
    mutate(time_on_trt = factor(0, levels = levels(data$time_on_trt)))
  data$mu_1 <- predict(outcome_model, newdata = data_trt_scenario)
  data$mu_0 <- predict(outcome_model, newdata = data_ctrl_scenario)
  
  data <- data %>%
    mutate(time_numeric = as.numeric(as.character(time)))
  
  tau_j <- numeric(J - ell)
  for (idx in seq_len(J - ell)) {
    j <- ell + idx
    data_j <- data %>% filter(time_numeric == j)
    
    pi_trt  <- prop_scores[j - ell]
    pi_ctrl <- 1 - prop_scores[j]
    
    I_trt  <- as.numeric(data_j$time_on_trt == max_exposure_time)
    I_ctrl <- as.numeric(data_j$trt == 0)
    
    trt_correction  <- ifelse(I_trt  == 1,
                              (data_j$out_cont_1 - data_j$mu_1) / pi_trt,
                              0)
    ctrl_correction <- ifelse(I_ctrl == 1,
                              (data_j$out_cont_1 - data_j$mu_0) / pi_ctrl,
                              0)
    
    summand <- data_j$mu_1 + trt_correction -
      data_j$mu_0 - ctrl_correction
    
    tau_j[idx] <- mean(summand) 
  }
  
  gate_estimate <- mean(tau_j)
  
  result <- list(
    estimate          = gate_estimate,
    outcome_model     = outcome_model,
    method            = "aipw",
    max_exposure_time = max_exposure_time,
    ell               = ell,
    prop_scores       = prop_scores
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