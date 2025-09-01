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
                                K = c(5, 10),
                                design_type = NULL,
                                efficient_matrix = NULL,
                                delta_small = c(0, 2),
                                delta_large = c(0, 2),
                                max_exposure_time = 1,
                                sigma_gamma = 1,
                                sigma_epsilon = 1,
                                r = 0.9,
                                omega_type = c("individual", "cluster"),
                                seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # match omega type
  omega_type <- match.arg(omega_type)
  
  # randomly assign cluster sizes
  half_I <- floor(I / 2)
  cluster_sizes <- rep(NA, I)
  
  # randomly select which clusters get K[1] vs K[2]
  small_clusters <- sample(1:I, half_I, replace = FALSE)
  cluster_sizes[small_clusters] <- K[1]
  cluster_sizes[-small_clusters] <- K[2]
  
  # design matrix
  X <- generate_design_matrix(I, J, design_type, efficient_matrix)
  
  # dtd structure for gamma (unit-period level)
  structure_decay <- matrix(NA, ncol = J, nrow = J)
  for (j in 1:(J - 1)) {
    for (k in (j + 1):J) {
      structure_decay[j, k] <- sigma_gamma * sigma_gamma * r^abs(j - k)
      structure_decay[k, j] <- structure_decay[j, k]
    }
  }
  
  for (j in 1:J){
    structure_decay[j, j] <- sigma_gamma * sigma_gamma
  }
  
  dtd <- MASS::mvrnorm(n = I, mu = rep(0, J), Sigma = structure_decay)
  
  # create long format data
  df_list <- list()
  row_idx <- 1
  
  # storage for omega aggregations
  omega_ij_matrix <- matrix(0, nrow = I, ncol = J)
  omega_j_vector <- rep(0, J)
  
  for (i in 1:I) {
    K_i <- cluster_sizes[i]
    for (j in 1:J) {
      omega_ij_sum <- 0
      for (k in 1:K_i) {
        # calculate omega_ijk based on type
        if (omega_type == "individual") {
          omega_ijk <- 1
        } else if (omega_type == "cluster") {
          omega_ijk <- 1 / K_i
        }
        
        omega_ij_sum <- omega_ij_sum + omega_ijk
        
        df_list[[row_idx]] <- data.frame(
          id_unit = as.factor(i),
          time = as.factor(j),
          id_ind = as.factor(k),
          out_cont_1 = j/J + dtd[i, j] + rnorm(1, 0, sigma_epsilon),
          trt = X[i, j],
          cluster_size = K_i,
          omega_ijk = omega_ijk
        )
        row_idx <- row_idx + 1
      }
      omega_ij_matrix[i, j] <- omega_ij_sum
      omega_j_vector[j] <- omega_j_vector[j] + omega_ij_sum
    }
  }
  
  df <- bind_rows(df_list)
  
  # add omega_ij and omega_j to dataframe
  df <- df %>%
    mutate(
      omega_ij = NA_real_,
      omega_j = NA_real_
    )
  
  for (i in 1:I) {
    for (j in 1:J) {
      df[df$id_unit == i & df$time == j, "omega_ij"] <- omega_ij_matrix[i, j]
      df[df$time == j, "omega_j"] <- omega_j_vector[j]
    }
  }
  
  # add treatment timing effects
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
      time_on_trt = pmin(time_on_trt_raw, max_exposure_time),
      # apply treatment effect based on cluster size
      trt_effect = ifelse(cluster_size == K[1], 
                          delta_small[time_on_trt + 1],
                          delta_large[time_on_trt + 1])
    ) %>%
    mutate(out_cont_1 = out_cont_1 + trt_effect) %>%
    mutate(time_on_trt = as.factor(time_on_trt)) %>%
    dplyr::select(-time_on_trt_raw, -trt_effect)
  
  attr(df, "design_type") <- design_type
  attr(df, "design_matrix") <- X
  attr(df, "cluster_sizes") <- cluster_sizes
  attr(df, "omega_type") <- omega_type
  attr(df, "omega_ij_matrix") <- omega_ij_matrix
  attr(df, "omega_j_vector") <- omega_j_vector
  attr(df, "parameters") <- list(
    I = I,
    J = J,
    K = K,
    delta_small = delta_small,
    delta_large = delta_large,
    max_exposure_time = max_exposure_time,
    sigma_gamma = sigma_gamma,
    sigma_epsilon = sigma_epsilon,
    r = r,
    omega_type = omega_type
  )
  
  return(df)
}

aipw_gate <- function(data) {
  
  max_exposure_time <- attr(data, "parameters")$max_exposure_time
  ell <- max_exposure_time - 1
  J <- attr(data, "parameters")$J
  I <- attr(data, "parameters")$I
  design_matrix <- attr(data, "design_matrix")
  
  # fit outcome regression model
  outcome_model <- lm(out_cont_1 ~ time_on_trt + time + cluster_size, 
                      data = data, weights = omega_ijk)
  
  # predictions
  data_treated_scenario <- data %>%
    mutate(time_on_trt = factor(max_exposure_time, levels = levels(data$time_on_trt)))
  data$pred_treated <- predict(outcome_model, newdata = data_treated_scenario)
  
  data_control_scenario <- data %>%
    mutate(time_on_trt = factor(0, levels = levels(data$time_on_trt)))
  data$pred_control <- predict(outcome_model, newdata = data_control_scenario)
  
  # cluster-period aggregations
  cluster_period_avg <- data %>%
    group_by(id_unit, time) %>%
    summarise(
      y_bar_ij = sum(omega_ijk * out_cont_1) / sum(omega_ijk),
      pred_treated_bar_ij = sum(omega_ijk * pred_treated) / sum(omega_ijk),
      pred_control_bar_ij = sum(omega_ijk * pred_control) / sum(omega_ijk),
      omega_ij = first(omega_ij),
      omega_j = first(omega_j),
      trt = first(trt),
      time_on_trt = first(time_on_trt),
      .groups = 'drop'
    ) %>%
    mutate(time_numeric = as.numeric(as.character(time)))
  
  # outcome regression component - treated
  or_treated_num <- 0
  or_treated_den <- 0
  
  for (j in 1:J) {
    data_j <- cluster_period_avg %>% 
      filter(time_numeric == j)
    
    if (nrow(data_j) > 0) {
      omega_j_val <- data_j$omega_j[1]
      pred_treated_j <- sum(data_j$omega_ij * data_j$pred_treated_bar_ij) / sum(data_j$omega_ij)
      
      or_treated_num <- or_treated_num + pred_treated_j * omega_j_val
      or_treated_den <- or_treated_den + omega_j_val
    }
  }
  
  # outcome regression component - control
  or_control_num <- 0
  or_control_den <- 0
  
  for (j in 1:J) {
    data_j <- cluster_period_avg %>% 
      filter(time_numeric == j)
    
    if (nrow(data_j) > 0) {
      omega_j_val <- data_j$omega_j[1]
      pred_control_j <- sum(data_j$omega_ij * data_j$pred_control_bar_ij) / sum(data_j$omega_ij)
      
      or_control_num <- or_control_num + pred_control_j * omega_j_val
      or_control_den <- or_control_den + omega_j_val
    }
  }
  
  or_estimate <- (or_treated_num / or_treated_den) - (or_control_num / or_control_den)
  
  # augmentation terms for treated units
  aug_treated_num <- 0
  aug_treated_den <- 0
  
  for (j in (ell+1):J) {
    data_j <- cluster_period_avg %>% 
      filter(time_numeric == j)
    
    treated_units <- data_j %>% 
      filter(time_on_trt == max_exposure_time)
    
    if (nrow(treated_units) > 0) {
      omega_j_val <- treated_units$omega_j[1]
      residuals <- treated_units$y_bar_ij - treated_units$pred_treated_bar_ij
      
      aug_treated_num <- aug_treated_num + 
        sum(treated_units$omega_ij * residuals) * omega_j_val / sum(treated_units$omega_ij)
      aug_treated_den <- aug_treated_den + omega_j_val
    }
  }
  
  # augmentation terms for control units
  aug_control_num <- 0
  aug_control_den <- 0
  
  for (j in 1:J) {
    data_j <- cluster_period_avg %>% 
      filter(time_numeric == j)
    
    control_units <- data_j %>% 
      filter(trt == 0)
    
    if (nrow(control_units) > 0) {
      omega_j_val <- control_units$omega_j[1]
      residuals <- control_units$y_bar_ij - control_units$pred_control_bar_ij
      
      aug_control_num <- aug_control_num + 
        sum(control_units$omega_ij * residuals) * omega_j_val / sum(control_units$omega_ij)
      aug_control_den <- aug_control_den + omega_j_val
    }
  }
  
  # AIPW estimate
  gate_estimate <- or_estimate + 
    (aug_treated_num / aug_treated_den) - 
    (aug_control_num / aug_control_den)
  
  result <- list(
    estimate = gate_estimate,
    outcome_model = outcome_model,
    method = "aipw",
    max_exposure_time = max_exposure_time,
    ell = ell,
    or_component = or_estimate,
    aug_treated = aug_treated_num / aug_treated_den,
    aug_control = aug_control_num / aug_control_den
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