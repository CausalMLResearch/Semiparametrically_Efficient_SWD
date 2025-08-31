source("util_empirical_cluster.R")

df <- read.csv("mobile_all.csv")
df$id_unit <- as.integer(df$id_unit)
df$time <- as.integer(df$time)
df$out_cont_1 <- as.numeric(df$out_cont_1)

results_dir <- "../../results/supp_figure_3/cluster_ell_1_12"
if (!dir.exists(results_dir)) dir.create(results_dir)

# simulation parameters
base_params <- list(
  I = 40,
  J_pre = 12,
  J_exp = 12,
  delta = c(0, 1, 2),
  max_exposure_time = 2
)

design_types <- c("t-opt", "ba", "ff", "ffba", "efficient")
n_sims <- 2000
true_effect <- base_params$delta[length(base_params$delta)]

# main simulation loop
for (sim_iter in 1:n_sims) {
  
  current_seed <- sim_iter
  cat(sprintf("Running simulation %d of %d (seed = %d)\n", sim_iter, n_sims, current_seed))
  
  seed_results <- list()
  seed_results$designs <- list()
  
  for (design in design_types) {
    
    cat(sprintf("  - Running design: %s\n", design))
    
    design_results <- list()
    design_results$design_type <- design
    
    tryCatch({
      # generate synthetic panel with treatment
      df_synthetic <- simulate_synthetic_panel(
        df = df,
        I = base_params$I,
        J_pre = base_params$J_pre,
        J_exp = base_params$J_exp,
        design_type = design,
        delta = base_params$delta,
        max_exposure_time = base_params$max_exposure_time,
        seed = current_seed
      )
      
      # estimate GATE using three methods
      or_result <- outcome_regression_gate(df_synthetic)
      design_results$outcome_regression <- list(
        estimate = or_result$estimate,
        method = "outcome_regression"
      )
      
      ipw_result <- ipw_gate(df_synthetic)
      design_results$ipw <- list(
        estimate = ipw_result$estimate,
        method = "ipw"
      )
      
      aipw_result <- aipw_gate(df_synthetic)
      design_results$aipw <- list(
        estimate = aipw_result$estimate,
        method = "aipw"
      )
      
      # create summary dataframe
      design_results$summary <- data.frame(
        design = design,
        method = c("outcome_regression", "ipw", "aipw"),
        estimate = c(or_result$estimate, ipw_result$estimate, aipw_result$estimate),
        true_effect = true_effect
      )
      
      design_results$status <- "success"
      design_results$parameters <- attr(df_synthetic, "parameters")
      
      # clean up memory
      rm(df_synthetic, or_result, ipw_result, aipw_result)
      
    }, error = function(e) {
      design_results$status <- "error"
      design_results$error_message <- as.character(e)
      warning(sprintf("Error in simulation %d, design %s: %s", sim_iter, design, e))
    })
    
    seed_results$designs[[design]] <- design_results
    gc(verbose = FALSE)
  }
  
  # combine summaries from all designs
  seed_results$combined_summary <- do.call(rbind, lapply(design_types, function(d) {
    if (!is.null(seed_results$designs[[d]]$summary)) {
      seed_results$designs[[d]]$summary
    }
  }))
  
  # save results for this seed
  filename <- file.path(results_dir, sprintf("sim_results_seed_%04d.rds", current_seed))
  saveRDS(seed_results, file = filename)
  
  rm(seed_results)
  gc(verbose = FALSE)
  
  # progress update every 100 simulations
  if (sim_iter %% 100 == 0) {
    cat(sprintf("Completed %d simulations (%.1f%%)\n", sim_iter, 100 * sim_iter / n_sims))
  }
}