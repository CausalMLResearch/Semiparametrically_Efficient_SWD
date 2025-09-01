source("util_unequal.R")

# SLURM array task handling
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == "") {
  task_id <- "1"
}
sim_id <- as.numeric(task_id)
seed <- sim_id

# output directory
output_dir <- "../../results/supp_table_6/unequal_ell_1"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(seed)
print(paste("Running simulation", sim_id, "with seed", seed))

# set up
base_params <- list(
  I = 32,
  J = 10,
  K = c(50, 100),
  delta_small = c(0, 1, 2),
  delta_large = c(0, 2, 4),
  max_exposure_time = 2,
  sigma_gamma = 1,
  sigma_epsilon = 1,
  r = 0.9,
  omega_type = "cluster"
)

design_types <- c("t-opt", "ba", "ff", "ffba", "efficient")

efficient_matrix <- c(0.3019901, 0.3019901, 0.3525803, 0.4159568, 0.4726380, 
                      0.5273620, 0.5840432, 0.6474198, 0.6980098, 0.6980098)

true_effect <- (base_params$delta_small[3] + base_params$delta_large[3]) / 2

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
      K = base_params$K,
      design_type = design,
      efficient_matrix = if(design == "efficient") efficient_matrix else NULL,
      delta_small = base_params$delta_small,
      delta_large = base_params$delta_large,
      max_exposure_time = base_params$max_exposure_time,
      sigma_gamma = base_params$sigma_gamma,
      sigma_epsilon = base_params$sigma_epsilon,
      r = base_params$r,
      omega_type = base_params$omega_type,
      seed = seed
    )
    
    aipw_result <- aipw_gate(df)
    design_results$aipw <- list(
      estimate = aipw_result$estimate,
      method = "aipw",
      or_component = aipw_result$or_component,
      aug_treated = aipw_result$aug_treated,
      aug_control = aipw_result$aug_control
    )
    
    design_results$summary <- data.frame(
      design = design,
      method = "aipw",
      estimate = aipw_result$estimate,
      true_effect = true_effect,
      bias = aipw_result$estimate - true_effect
    )
    
    design_results$status <- "success"
    
    rm(df, aipw_result)
    
  }, error = function(e) {
    design_results$status <- "error"
    design_results$error_message <- as.character(e)
    warning(sprintf("Error in simulation %d, design %s: %s", sim_id, design, e))
  })
  
  seed_results$designs[[design]] <- design_results
  
  gc(verbose = FALSE)
}

seed_results$combined_summary <- do.call(rbind, lapply(design_types, function(d) {
  if (!is.null(seed_results$designs[[d]]$summary)) {
    seed_results$designs[[d]]$summary
  }
}))

seed_results$simulation_id <- sim_id
seed_results$seed <- seed

filename <- file.path(output_dir, sprintf("sim_results_seed_%04d.rds", seed))
saveRDS(seed_results, file = filename)

print(paste("Simulation", sim_id, "finished and saved to", filename))