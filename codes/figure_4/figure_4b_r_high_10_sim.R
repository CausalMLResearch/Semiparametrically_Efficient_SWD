source("util_IRT_sim.R")
results_dir <- "../../results/figure_4/r_high_10"

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
base_params <- list(
  I = 32,
  J = 10,
  delta = c(0, 2),
  max_exposure_time = 1,
  sigma_epsilon = 1,
  r = 0.9
)
design_types <- c("t-opt", "ba", "ff", "ffba", "efficient_9", "efficient_7", "efficient_5")
efficient_matrix_9 <- c(0.1638731, 0.2765203, 0.3531312, 0.4159630, 0.4725966, 
                        0.5274034, 0.5840370, 0.6468688, 0.7234797, 0.8361270)
efficient_matrix_7 <- c(0.2264541, 0.3420778, 0.4050425, 0.4485636, 0.4836648, 
                        0.5163352, 0.5514364, 0.5949575, 0.6579222, 0.7735461)
efficient_matrix_5 <- c(0.2972705, 0.4048493, 0.4507607, 0.4759553, 0.4927645, 
                        0.5072355, 0.5240447, 0.5492393, 0.5951507, 0.7027295)
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
                 "efficient_9" = efficient_matrix_9,
                 "efficient_7" = efficient_matrix_7,
                 "efficient_5" = efficient_matrix_5
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