library(dplyr)
library(lme4)
library(tidyverse)
library(MASS)
library(parallel)
library(doSNOW)
library(foreach)
library(doRNG)

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
  
  max_exposure_time <- attr(data, "parameters")$max_exposure_time
  
  model <- lm(out_cont_1 ~ time_on_trt, data = data)
  
  coef_names <- names(coef(model))
  max_exposure_time_coef_name <- paste0("time_on_trt", max_exposure_time)
  
  if (max_exposure_time_coef_name %in% coef_names) {
    gate_estimate <- coef(model)[max_exposure_time_coef_name]
  } else {
    stop(paste("Coefficient for time_on_trt =", max_exposure_time, "not found"))
  }
  
  result <- list(
    estimate = as.numeric(gate_estimate),
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

aipw_gate_1 <- function(data) {
  
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
  outcome_model <- lm(out_cont_1 ~ time_on_trt, data = data)
  
  # for treated scenario: predict with time_on_trt = max_exposure_time
  data_treated_scenario <- data %>%
    mutate(time_on_trt = factor(max_exposure_time, levels = levels(data$time_on_trt)))
  data$pred_treated <- predict(outcome_model, newdata = data_treated_scenario)
  
  # for control scenario: predict with time_on_trt = 0
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
      
      aipw_sum_treated <- aipw_sum_treated + sum(treated_units$out_cont_1 / pi_j_minus_ell) 
      + sum((1 - pi_j_minus_ell) * treated_units$pred_treated / pi_j_minus_ell)
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
      
      aipw_sum_control <- aipw_sum_control + sum(control_units$out_cont_1 / (1 - pi_j)) 
      + sum((0 - pi_j) * control_units$pred_control / (1 - pi_j))
      denom_control <- denom_control + nrow(control_units) / (1 - pi_j)
    }
  }
  
  gate_estimate <- aipw_sum_treated / denom_treated - aipw_sum_control / denom_control
  
  result <- list(
    estimate = gate_estimate,
    outcome_model = outcome_model,
    method = "aipw1",
    max_exposure_time = max_exposure_time,
    ell = ell,
    prop_scores = prop_scores
  )
  
  class(result) <- c("gate_estimate", "list")
  return(result)
}

aipw_gate_2 <- function(data) {
  
  max_exposure_time <- attr(data, "parameters")$max_exposure_time
  ell <- max_exposure_time - 1
  J <- attr(data, "parameters")$J
  I <- attr(data, "parameters")$I
  design_matrix <- attr(data, "design_matrix")
  
  # true propensity scores
  prop_scores <- colMeans(design_matrix)
  
  data <- data %>%
    mutate(time_numeric = as.numeric(as.character(time)),
           prop_score = prop_scores[time_numeric])
  
  # outcome regression model
  outcome_model <- lm(out_cont_1 ~ time_on_trt, data = data)
  
  coef_names <- names(coef(outcome_model))
  max_exposure_time_coef_name <- paste0("time_on_trt", max_exposure_time)
  
  if (max_exposure_time_coef_name %in% coef_names) {
    outcome_model_estimate <- as.numeric(coef(outcome_model)[max_exposure_time_coef_name])
  } else {
    stop(paste("Coefficient for time_on_trt =", max_exposure_time, "not found"))
  }
  
  # for treated scenario: predict with time_on_trt = max_exposure_time
  data_treated_scenario <- data %>%
    mutate(time_on_trt = factor(max_exposure_time, levels = levels(data$time_on_trt)))
  data$pred_treated <- predict(outcome_model, newdata = data_treated_scenario)
  
  # for control scenario: predict with time_on_trt = 0
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
  
  gate_estimate <- outcome_model_estimate + aipw_sum_treated / denom_treated - aipw_sum_control / denom_control
  
  result <- list(
    estimate = gate_estimate,
    outcome_model = outcome_model,
    method = "aipw2",
    max_exposure_time = max_exposure_time,
    ell = ell,
    prop_scores = prop_scores
  )
  
  class(result) <- c("gate_estimate", "list")
  return(result)
}

run_two_period_simulation <- function(pi_1_seq = seq(0.01, 0.99, by = 0.01),
                                      pi_2_seq = seq(0.01, 0.99, by = 0.01),
                                      n_sims = 200,
                                      I = 200,
                                      seed = 123,
                                      save_intermediate = TRUE,
                                      checkpoint_every = 100) {
  
  set.seed(seed)
  
  J <- 2
  delta <- c(0, 2)
  max_exposure_time <- 1
  sigma_epsilon <- 1
  r <- 0.9
  
  true_gate <- delta[max_exposure_time + 1]
  
  prop_grid <- expand.grid(pi_1 = pi_1_seq, pi_2 = pi_2_seq)
  prop_grid <- prop_grid[prop_grid$pi_1 <= prop_grid$pi_2, ]
  n_settings <- nrow(prop_grid)
  
  cat(sprintf("Running simulation on %d parameter combinations\n", n_settings))
  cat(sprintf("Total simulations to run: %d\n", n_settings * n_sims))
  
  n_cores <- min(detectCores() - 1, 15)
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(min = 0, max = n_settings, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  registerDoRNG(seed)
  
  clusterEvalQ(cl, {
    library(dplyr)
    library(tidyverse)
    library(MASS)
  })
  
  clusterExport(cl, c("generate_design_matrix", "simulate_panel_data", 
                      "outcome_regression_gate", "ipw_gate", 
                      "aipw_gate_1", "aipw_gate_2"))
  
  cat(sprintf("Using %d cores for parallel processing\n", n_cores))
  cat("Starting simulation...\n")
  start_time <- Sys.time()
  
  results <- foreach(i = 1:n_settings, 
                     .combine = rbind,
                     .packages = c("dplyr", "tidyverse", "MASS"),
                     .options.snow = opts) %dorng% {
                       
                       pi_1 <- prop_grid$pi_1[i]
                       pi_2 <- prop_grid$pi_2[i]
                       
                       est_or <- numeric(n_sims)
                       est_ipw <- numeric(n_sims)
                       est_aipw1 <- numeric(n_sims)
                       est_aipw2 <- numeric(n_sims)
                       
                       for (sim in 1:n_sims) {
                         tryCatch({
                           data <- simulate_panel_data(
                             I = I,
                             J = J,
                             design_type = "efficient",
                             efficient_matrix = c(pi_1, pi_2),
                             delta = delta,
                             max_exposure_time = max_exposure_time,
                             sigma_epsilon = sigma_epsilon,
                             r = r,
                             seed = NULL
                           )
                           
                           est_or[sim] <- outcome_regression_gate(data)$estimate
                           est_ipw[sim] <- ipw_gate(data)$estimate
                           est_aipw1[sim] <- aipw_gate_1(data)$estimate
                           est_aipw2[sim] <- aipw_gate_2(data)$estimate
                         }, error = function(e) {
                           est_or[sim] <- NA
                           est_ipw[sim] <- NA
                           est_aipw1[sim] <- NA
                           est_aipw2[sim] <- NA
                         })
                       }
                       
                       # save intermediate results
                       if (save_intermediate && (i %% checkpoint_every == 0)) {
                         checkpoint_data <- data.frame(
                           checkpoint = i,
                           time = Sys.time(),
                           completed = i,
                           total = n_settings
                         )
                         write.csv(checkpoint_data, 
                                   paste0("checkpoint_", i, ".csv"), 
                                   row.names = FALSE)
                       }
                       
                       data.frame(
                         pi_1 = pi_1,
                         pi_2 = pi_2,
                         mse_or = mean((est_or - true_gate)^2, na.rm = TRUE),
                         mse_ipw = mean((est_ipw - true_gate)^2, na.rm = TRUE),
                         mse_aipw1 = mean((est_aipw1 - true_gate)^2, na.rm = TRUE),
                         mse_aipw2 = mean((est_aipw2 - true_gate)^2, na.rm = TRUE),
                         abias_or = mean(abs(est_or - true_gate), na.rm = TRUE),
                         abias_ipw = mean(abs(est_ipw - true_gate), na.rm = TRUE),
                         abias_aipw1 = mean(abs(est_aipw1 - true_gate), na.rm = TRUE),
                         abias_aipw2 = mean(abs(est_aipw2 - true_gate), na.rm = TRUE),
                         var_or = var(est_or, na.rm = TRUE),
                         var_ipw = var(est_ipw, na.rm = TRUE),
                         var_aipw1 = var(est_aipw1, na.rm = TRUE),
                         var_aipw2 = var(est_aipw2, na.rm = TRUE)
                       )
                     }
  
  close(pb)
  stopCluster(cl)
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "hours")
  
  cat(sprintf("\nSimulation completed in %.2f hours!\n", as.numeric(runtime)))
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename <- paste0("../../results/figure_1/simulation_results_", timestamp, ".csv")
  write.csv(results, filename, row.names = FALSE)
  cat(sprintf("results saved to %s\n", filename))
  
  runtime_info <- data.frame(
    start_time = start_time,
    end_time = end_time,
    runtime_hours = as.numeric(runtime),
    n_settings = n_settings,
    n_sims = n_sims,
    total_sims = n_settings * n_sims,
    n_cores = n_cores
  )
  write.csv(runtime_info, paste0("runtime_info_", timestamp, ".csv"), row.names = FALSE)
  
  return(results)
}

results <- run_two_period_simulation(
  pi_1_seq = seq(0.10, 0.95, by = 0.01),
  pi_2_seq = seq(0.10, 0.95, by = 0.01),
  n_sims = 2000,
  I = 100,
  seed = 123,
  save_intermediate = TRUE,
  checkpoint_every = 20
)