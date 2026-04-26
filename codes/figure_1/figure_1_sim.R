library(dplyr)
library(tidyverse)
library(MASS)
library(parallel)
library(doSNOW)
library(foreach)
library(doRNG)

# =============================================================================
# Design matrix and data generation
# =============================================================================
generate_design_matrix <- function(I, J,
                                   design_type = "t-opt",
                                   efficient_matrix = NULL) {
  
  X <- matrix(0, nrow = I, ncol = J)
  
  if (design_type == "t-opt") {
    for (j in 1:J) {
      treated_fraction <- (2 * j - 1) / (2 * J)
      n_treated <- round(treated_fraction * I)
      if (n_treated > 0) X[1:n_treated, j] <- 1
    }
  } else if (design_type == "ba") {
    X[, (floor(J / 2) + 1):J] <- 1
    X[1, ] <- 1
  } else if (design_type == "ff") {
    half_I <- floor(I / 2)
    X[(half_I + 1):I, ] <- 1
  } else if (design_type == "ffba") {
    half_J <- floor(J / 2)
    half_I <- floor(I / 2)
    X[(half_I + 1):I, (half_J + 1):J] <- 1
  } else if (grepl("^efficient", design_type)) {
    for (j in 1:J) {
      n_treated <- round(efficient_matrix[j] * I)
      if (n_treated > 0) X[1:n_treated, j] <- 1
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
  
  if (!is.null(seed)) set.seed(seed)
  
  X <- generate_design_matrix(I, J, design_type, efficient_matrix)
  X_long <- as.vector(t(X))
  
  # Toeplitz residual covariance: Sigma_{jj'} = sigma^2 * r^|j-j'|
  structure_decay <- matrix(0, nrow = J, ncol = J)
  for (j in 1:J) for (k in 1:J) {
    structure_decay[j, k] <- sigma_epsilon^2 * r^abs(j - k)
  }
  
  dtd <- MASS::mvrnorm(n = I, mu = rep(0, J), Sigma = structure_decay)
  
  # True DGP (Figure 1 caption): Y_ij = beta_j + delta[time_on_trt+1] + eps_ij
  # with beta_j = j/J
  Y_ij <- numeric(I * J)
  for (i in 1:I) for (j in 1:J) {
    Y_ij[(i - 1) * J + j] <- j / J + dtd[i, j]
  }
  
  df <- data.frame(
    id_unit    = as.factor(rep(1:I, each = J)),
    time       = as.factor(rep(1:J, times = I)),
    out_cont_1 = Y_ij,
    trt        = X_long
  )
  
  df_adoption <- df %>%
    filter(trt == 1) %>%
    group_by(id_unit) %>%
    summarise(first_time = min(as.numeric(as.character(time))), .groups = "drop")
  
  df <- df %>%
    left_join(df_adoption, by = "id_unit") %>%
    mutate(
      time_on_trt_raw = ifelse(is.na(first_time) | trt == 0, 0,
                               as.numeric(as.character(time)) - first_time + 1),
      time_on_trt     = pmin(time_on_trt_raw, max_exposure_time)
    ) %>%
    mutate(out_cont_1 = out_cont_1 + delta[time_on_trt + 1]) %>%
    mutate(time_on_trt = as.factor(time_on_trt)) %>%
    dplyr::select(-time_on_trt_raw)
  
  attr(df, "design_type")   <- design_type
  attr(df, "design_matrix") <- X
  attr(df, "parameters")    <- list(
    I = I, J = J, delta = delta, max_exposure_time = max_exposure_time,
    sigma_epsilon = sigma_epsilon, r = r
  )
  
  return(df)
}


# =============================================================================
# Misspecified working outcome model per Figure 1 caption:
#     Y_hat_ij = mu_hat + tau_0_hat * Z_ij
# Fit by OLS on all (i, j) pooled.  Omits the period trend beta_j that exists
# in the true DGP -- this misspecification is the point of the figure.
# =============================================================================
.fit_working_model <- function(data) {
  lm(out_cont_1 ~ trt, data = data)
}


# =============================================================================
# Outcome regression estimator for tau_0 (ell = 0 here)
#   tau_0 = (1/J) sum_{j=1}^J (1/I) sum_i [mu_1(X_i, j) - mu_0(X_i, j)]
# Under the misspecified working model Y_hat = mu + tau_0 * Z,
#   mu_1 - mu_0 = tau_0_hat (constant), so this equals the OLS coefficient.
# =============================================================================
outcome_regression_gate <- function(data) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
  
  model <- .fit_working_model(data)
  
  data_trt_scenario  <- data %>% mutate(trt = 1)
  data_ctrl_scenario <- data %>% mutate(trt = 0)
  data$mu_1 <- predict(model, newdata = data_trt_scenario)
  data$mu_0 <- predict(model, newdata = data_ctrl_scenario)
  
  data <- data %>% mutate(time_numeric = as.numeric(as.character(time)))
  
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
    ell               = ell
  )
  class(result) <- c("gate_estimate", "list")
  return(result)
}


# =============================================================================
# IPW estimator for tau_0 (Horvitz-Thompson)
#   tau_{j, 0} = (1/I) sum_i { Z_{ij} / pi_j * Y_ij - (1-Z_{ij})/(1-pi_j) * Y_ij }
# averaged over j = 1, ..., J.
# =============================================================================
ipw_gate <- function(data) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
  design_matrix     <- attr(data, "design_matrix")
  
  prop_scores <- colMeans(design_matrix)
  
  data <- data %>% mutate(time_numeric = as.numeric(as.character(time)))
  
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


# =============================================================================
# AIPW estimator for tau_0 with MISSPECIFIED working outcome model
# (paper Section 2.4 formula, line 370, plus Figure 1 working-model spec)
# =============================================================================
aipw_gate <- function(data) {
  
  params            <- attr(data, "parameters")
  max_exposure_time <- params$max_exposure_time
  ell               <- max_exposure_time - 1
  J                 <- params$J
  design_matrix     <- attr(data, "design_matrix")
  
  prop_scores <- colMeans(design_matrix)
  
  outcome_model <- .fit_working_model(data)
  
  data_trt_scenario  <- data %>% mutate(trt = 1)
  data_ctrl_scenario <- data %>% mutate(trt = 0)
  data$mu_1 <- predict(outcome_model, newdata = data_trt_scenario)
  data$mu_0 <- predict(outcome_model, newdata = data_ctrl_scenario)
  
  data <- data %>% mutate(time_numeric = as.numeric(as.character(time)))
  
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


# =============================================================================
# Two-period grid simulation driver
# =============================================================================
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
  max_exposure_time <- 1       # paper's ell = 0 (no carryover; target is tau_0)
  sigma_epsilon <- 1
  r <- 0.9
  
  true_gate <- delta[max_exposure_time + 1]   # = delta[2] = 2
  
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
                      ".fit_working_model",
                      "outcome_regression_gate", "ipw_gate", "aipw_gate"))
  
  cat(sprintf("Using %d cores for parallel processing\n", n_cores))
  cat("Starting simulation...\n")
  start_time <- Sys.time()
  
  results <- foreach(i = 1:n_settings,
                     .combine = rbind,
                     .packages = c("dplyr", "tidyverse", "MASS"),
                     .options.snow = opts) %dorng% {
                       
                       pi_1 <- prop_grid$pi_1[i]
                       pi_2 <- prop_grid$pi_2[i]
                       
                       est_or   <- numeric(n_sims)
                       est_ipw  <- numeric(n_sims)
                       est_aipw <- numeric(n_sims)
                       
                       for (sim in 1:n_sims) {
                         tryCatch({
                           data <- simulate_panel_data(
                             I                 = I,
                             J                 = J,
                             design_type       = "efficient",
                             efficient_matrix  = c(pi_1, pi_2),
                             delta             = delta,
                             max_exposure_time = max_exposure_time,
                             sigma_epsilon     = sigma_epsilon,
                             r                 = r,
                             seed              = NULL
                           )
                           
                           est_or[sim]   <- outcome_regression_gate(data)$estimate
                           est_ipw[sim]  <- ipw_gate(data)$estimate
                           est_aipw[sim] <- aipw_gate(data)$estimate
                         }, error = function(e) {
                           est_or[sim]   <<- NA
                           est_ipw[sim]  <<- NA
                           est_aipw[sim] <<- NA
                         })
                       }
                       
                       if (save_intermediate && (i %% checkpoint_every == 0)) {
                         checkpoint_data <- data.frame(
                           checkpoint = i, time = Sys.time(),
                           completed  = i, total = n_settings
                         )
                         write.csv(checkpoint_data,
                                   paste0("checkpoint_", i, ".csv"),
                                   row.names = FALSE)
                       }
                       
                       data.frame(
                         pi_1       = pi_1,
                         pi_2       = pi_2,
                         mse_or     = mean((est_or   - true_gate)^2, na.rm = TRUE),
                         mse_ipw    = mean((est_ipw  - true_gate)^2, na.rm = TRUE),
                         mse_aipw   = mean((est_aipw - true_gate)^2, na.rm = TRUE),
                         abias_or   = mean(abs(est_or   - true_gate), na.rm = TRUE),
                         abias_ipw  = mean(abs(est_ipw  - true_gate), na.rm = TRUE),
                         abias_aipw = mean(abs(est_aipw - true_gate), na.rm = TRUE),
                         var_or     = var(est_or,   na.rm = TRUE),
                         var_ipw    = var(est_ipw,  na.rm = TRUE),
                         var_aipw   = var(est_aipw, na.rm = TRUE)
                       )
                     }
  
  close(pb)
  stopCluster(cl)
  
  end_time <- Sys.time()
  runtime  <- difftime(end_time, start_time, units = "hours")
  
  cat(sprintf("\nSimulation completed in %.2f hours!\n", as.numeric(runtime)))
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  filename  <- paste0("simulation_results_", timestamp, ".csv")
  write.csv(results, filename, row.names = FALSE)
  cat(sprintf("results saved to %s\n", filename))
  
  runtime_info <- data.frame(
    start_time    = start_time,
    end_time      = end_time,
    runtime_hours = as.numeric(runtime),
    n_settings    = n_settings,
    n_sims        = n_sims,
    total_sims    = n_settings * n_sims,
    n_cores       = n_cores
  )
  write.csv(runtime_info,
            paste0("runtime_info_", timestamp, ".csv"),
            row.names = FALSE)
  
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