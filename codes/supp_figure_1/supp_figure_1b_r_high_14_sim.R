source("util_IRT_sim.R")
results_dir <- "../../results/supp_figure_1/r_high_14"

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
base_params <- list(
  I = 1000,
  J = 14,
  delta = c(0, 2),
  max_exposure_time = 1,
  sigma_epsilon = 1,
  r = 0.9
)
design_types <- c("t-opt", "ba", "ff", "ffba", "efficient")
efficient_matrix <- c(0.1469673, 0.2460769, 0.3111981, 0.3623029, 0.4058288, 
                      0.4450129, 0.4819100, 0.5180899, 0.5549871, 0.5941712, 
                      0.6376971, 0.6888019, 0.7539231, 0.8530326)
n_sims <- 2000
true_effect <- base_params$delta[length(base_params$delta)]
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
                 "efficient" = efficient_matrix
          )
        } else NULL,
        delta = base_params$delta,
        max_exposure_time = base_params$max_exposure_time,
        sigma_epsilon = base_params$sigma_epsilon,
        r = base_params$r,
        seed = current_seed
      )
      
      or_result <- outcome_regression_gate(df)
      design_results$outcome_regression <- list(
        estimate = or_result$estimate,
        method = "outcome_regression"
      )
      
      ipw_result <- ipw_gate(df)
      design_results$ipw <- list(
        estimate = ipw_result$estimate,
        method = "ipw"
      )
      
      aipw_result <- aipw_gate(df)
      design_results$aipw <- list(
        estimate = aipw_result$estimate,
        method = "aipw"
      )
      
      design_results$summary <- data.frame(
        design = design,
        method = c("outcome_regression", "ipw", "aipw"),
        estimate = c(or_result$estimate, ipw_result$estimate, aipw_result$estimate),
        true_effect = true_effect
      )
      
      design_results$status <- "success"
      
      rm(df, or_result, ipw_result, aipw_result)
      
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