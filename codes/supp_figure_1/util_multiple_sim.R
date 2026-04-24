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
    tau_j[idx] <- mean(data_j$mu_1 - data_j$mu_0)
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

outcome_regression_cte <- function(data, s = 0) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
  
  if (s < 0 || s > ell) stop("s must satisfy 0 <= s <= ell")
  
  model <- lm(out_cont_1 ~ time_on_trt + time, data = data)
  
  # treated scenario: exposed for s+1 periods (time_on_trt = s+1)
  data_trt_scenario <- data %>%
    mutate(time_on_trt = factor(s + 1, levels = levels(data$time_on_trt)))
  data_ctrl_scenario <- data %>%
    mutate(time_on_trt = factor(0, levels = levels(data$time_on_trt)))
  data$mu_1 <- predict(model, newdata = data_trt_scenario)
  data$mu_0 <- predict(model, newdata = data_ctrl_scenario)
  
  data <- data %>%
    mutate(time_numeric = as.numeric(as.character(time)))
  
  # outer average over j = s+1, ..., J
  tau_j <- numeric(J - s)
  for (idx in seq_len(J - s)) {
    j <- s + idx
    data_j <- data %>% filter(time_numeric == j)
    tau_j[idx] <- mean(data_j$mu_1 - data_j$mu_0)
  }
  
  cte_estimate <- mean(tau_j)
  
  result <- list(
    estimate          = cte_estimate,
    model             = model,
    method            = "outcome_regression",
    max_exposure_time = max_exposure_time,
    ell               = ell,
    s                 = s,
    all_coef          = coef(model)
  )
  
  class(result) <- c("cte_estimate", "list")
  return(result)
}

ipw_gate <- function(data) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
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
    
    tau_j[idx] <- mean(trt_term - ctrl_term)
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

ipw_cte <- function(data, s = 0) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
  design_matrix     <- attr(data, "design_matrix")
  
  if (s < 0 || s > ell) stop("s must satisfy 0 <= s <= ell")
  
  prop_scores <- colMeans(design_matrix)
  
  data <- data %>%
    mutate(time_numeric = as.numeric(as.character(time)))
  
  # outer average over j = s+1, ..., J
  tau_j <- numeric(J - s)
  for (idx in seq_len(J - s)) {
    j <- s + idx
    
    # treated weight: pi_{j-s} - pi_{j-s-1}  (denominator of lambda * pi in paper eq)
    # convention: pi_0 = 0 (nobody treated before trial start)
    pi_jms  <- prop_scores[j - s]
    pi_jms1 <- if (j - s - 1 >= 1) prop_scores[j - s - 1] else 0
    pi_trt  <- pi_jms - pi_jms1
    pi_ctrl <- 1 - prop_scores[j]
    
    data_j <- data %>% filter(time_numeric == j)
    
    # treated arm: units with time_on_trt == s+1 (switchers who've been exposed for exactly s+1 periods)
    I_trt  <- as.numeric(data_j$time_on_trt == s + 1)
    I_ctrl <- as.numeric(data_j$trt == 0)
    
    trt_term  <- ifelse(I_trt  == 1, data_j$out_cont_1 / pi_trt,  0)
    ctrl_term <- ifelse(I_ctrl == 1, data_j$out_cont_1 / pi_ctrl, 0)
    
    tau_j[idx] <- mean(trt_term - ctrl_term)
  }
  
  cte_estimate <- mean(tau_j)
  
  result <- list(
    estimate          = cte_estimate,
    method            = "ipw",
    max_exposure_time = max_exposure_time,
    ell               = ell,
    s                 = s,
    prop_scores       = prop_scores
  )
  
  class(result) <- c("cte_estimate", "list")
  return(result)
}

aipw_gate <- function(data) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
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

aipw_cte <- function(data, s = 0) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
  design_matrix     <- attr(data, "design_matrix")
  
  if (s < 0 || s > ell) stop("s must satisfy 0 <= s <= ell")
  
  prop_scores <- colMeans(design_matrix)
  
  outcome_model <- lm(out_cont_1 ~ time_on_trt + time, data = data)
  
  # treated scenario: exposed for s+1 periods (time_on_trt = s+1)
  data_trt_scenario <- data %>%
    mutate(time_on_trt = factor(s + 1, levels = levels(data$time_on_trt)))
  data_ctrl_scenario <- data %>%
    mutate(time_on_trt = factor(0, levels = levels(data$time_on_trt)))
  data$mu_1 <- predict(outcome_model, newdata = data_trt_scenario)
  data$mu_0 <- predict(outcome_model, newdata = data_ctrl_scenario)
  
  data <- data %>%
    mutate(time_numeric = as.numeric(as.character(time)))
  
  # outer average over j = s+1, ..., J
  tau_j <- numeric(J - s)
  for (idx in seq_len(J - s)) {
    j <- s + idx
    
    pi_jms  <- prop_scores[j - s]
    pi_jms1 <- if (j - s - 1 >= 1) prop_scores[j - s - 1] else 0
    pi_trt  <- pi_jms - pi_jms1     # pi_{j-s} - pi_{j-s-1}
    pi_ctrl <- 1 - prop_scores[j]
    
    data_j <- data %>% filter(time_numeric == j)
    
    I_trt  <- as.numeric(data_j$time_on_trt == s + 1)
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
  
  cte_estimate <- mean(tau_j)
  
  result <- list(
    estimate          = cte_estimate,
    outcome_model     = outcome_model,
    method            = "aipw",
    max_exposure_time = max_exposure_time,
    ell               = ell,
    s                 = s,
    prop_scores       = prop_scores
  )
  
  class(result) <- c("cte_estimate", "list")
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

summarize_all_designs_by_type <- function(results_dir = "simulation_results_all_designs", seeds = 1:2000) {
  
  all_results <- load_all_design_results(results_dir, seeds)
  
  all_summaries <- do.call(rbind, lapply(names(all_results), function(seed) {
    res <- all_results[[seed]]
    if (!is.null(res$combined_summary)) {
      cbind(seed = as.integer(seed), res$combined_summary)
    }
  }))
  
  gate_stats <- all_summaries %>%
    filter(effect_type == "gate") %>%
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
  
  cte_stats <- all_summaries %>%
    filter(effect_type == "cte") %>%
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
    gate_stats = gate_stats,
    cte_stats = cte_stats,
    raw_summaries = all_summaries
  ))
}