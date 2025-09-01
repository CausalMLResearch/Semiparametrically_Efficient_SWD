source("util_multiple_sim.R")
results_dir <- "../../results/figure_4/multiple_estimands"

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
base_params <- list(
  I = 1000,
  J = 10,
  delta = c(0, 1, 2),
  max_exposure_time = 2,
  sigma_epsilon = 1,
  r = 0.9
)
design_types <- c("t-opt", "ba", "ff", "ffba", "efficient_0", "efficient_1")
efficient_matrix_0 <- c(0.3012258, 0.3012258, 0.3522969, 0.4158044, 0.4725896, 
                        0.5274104, 0.5841956, 0.6477031, 0.6987741, 0.6987741)
efficient_matrix_1 <- c(0.1323045, 0.2212771, 0.3024931, 0.3790078, 0.4531267, 
                        0.5265982, 0.6010467, 0.6784067, 0.7616822, 0.8572062)

n_sims <- 2000

true_effect_gate <- base_params$delta[base_params$max_exposure_time + 1]  # delta[3] = 2
true_effect_cte <- base_params$delta[2]  # delta[2] = 1

for (sim_iter in 1:n_sims) {
  
  current_seed <- sim_iter
  
  cat(sprintf("Running simulation %d of %d (seed = %d)\n", sim_iter, n_sims, current_seed))
  
  seed_results <- list()
  seed_results$designs <- list()
  
  for (design in design_types) {
    
    cat(sprintf(" - Running design: %s\n", design))
    
    design_results <- list()
    design_results$design_type <- design
    
    tryCatch({
      df <- simulate_panel_data(
        I = base_params$I,
        J = base_params$J,
        design_type = design,
        efficient_matrix = if(grepl("efficient", design)) {
          switch(design,
                 "efficient_0" = efficient_matrix_0,
                 "efficient_1" = efficient_matrix_1
          )
        } else NULL,
        delta = base_params$delta,
        max_exposure_time = base_params$max_exposure_time,
        sigma_epsilon = base_params$sigma_epsilon,
        r = base_params$r,
        seed = current_seed
      )
      
      # GATE estimates
      or_gate_result <- outcome_regression_gate(df)
      design_results$outcome_regression_gate <- list(
        estimate = or_gate_result$estimate,
        method = "outcome_regression",
        effect_type = "gate"
      )
      
      ipw_gate_result <- ipw_gate(df)
      design_results$ipw_gate <- list(
        estimate = ipw_gate_result$estimate,
        method = "ipw",
        effect_type = "gate"
      )
      
      aipw_gate_result <- aipw_gate(df)
      design_results$aipw_gate <- list(
        estimate = aipw_gate_result$estimate,
        method = "aipw",
        effect_type = "gate"
      )
      
      # CTE estimates
      or_cte_result <- outcome_regression_cte(df)
      design_results$outcome_regression_cte <- list(
        estimate = or_cte_result$estimate,
        method = "outcome_regression",
        effect_type = "cte"
      )
      
      ipw_cte_result <- ipw_cte(df)
      design_results$ipw_cte <- list(
        estimate = ipw_cte_result$estimate,
        method = "ipw",
        effect_type = "cte"
      )
      
      aipw_cte_result <- aipw_cte(df)
      design_results$aipw_cte <- list(
        estimate = aipw_cte_result$estimate,
        method = "aipw",
        effect_type = "cte"
      )
      
      design_results$summary <- data.frame(
        design = design,
        method = rep(c("outcome_regression", "ipw", "aipw"), 2),
        effect_type = c(rep("gate", 3), rep("cte", 3)),
        estimate = c(or_gate_result$estimate, ipw_gate_result$estimate, aipw_gate_result$estimate,
                     or_cte_result$estimate, ipw_cte_result$estimate, aipw_cte_result$estimate),
        true_effect = c(rep(true_effect_gate, 3), rep(true_effect_cte, 3))
      )
      
      design_results$status <- "success"
      
      rm(df, or_gate_result, ipw_gate_result, aipw_gate_result,
         or_cte_result, ipw_cte_result, aipw_cte_result)
      
    }, error = function(e) {
      design_results$status <- "error"
      design_results$error_message <- as.character(e)
      warning(sprintf("Error in simulation %d, design %s: %s", sim_iter, design, e))
    })
    
    seed_results$designs[[design]] <- design_results
    
    gc(verbose = FALSE)
  }
  
  seed_results$combined_summary <- do.call(rbind, lapply(design_types, function(d) {
    if (!is.null(seed_results$designs[[d]]$summary)) {
      seed_results$designs[[d]]$summary
    }
  }))
  
  filename <- file.path(results_dir, sprintf("sim_results_seed_%04d.rds", current_seed))
  saveRDS(seed_results, file = filename)
  
  rm(seed_results)
  gc(verbose = FALSE)
  
  if (sim_iter %% 100 == 0) {
    cat(sprintf("completed %d simulations (%.1f%%)\n", sim_iter, 100 * sim_iter / n_sims))
  }
}