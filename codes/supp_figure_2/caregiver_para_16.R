source("util_para.R")

# load data
df <- read.csv("mobile_caregiver.csv")
df$id_unit <- as.integer(df$id_unit)
df$time <- as.integer(df$time)
df$out_cont_1 <- as.numeric(df$out_cont_1)

# create results directory
results_dir <- "../../results/supp_figure_2/caregiver_para"
if (!dir.exists(results_dir)) dir.create(results_dir)

# fixed parameters
I <- 40
J_pre <- 16
J_exp <- 16
delta <- c(0, 2)
max_exposure_time <- 1
n_sims <- 2000

# main simulation loop
for (sim_iter in 1:n_sims) {
  
  current_seed <- sim_iter
  cat(sprintf("running simulation %d of %d (seed = %d)\n", sim_iter, n_sims, current_seed))
  
  seed_results <- list()
  
  tryCatch({
    # generate synthetic panel with efficient design
    df_synthetic <- simulate_synthetic_panel(
      df = df,
      I = I,
      J_pre = J_pre,
      J_exp = J_exp,
      design_type = "efficient",
      delta = delta,
      max_exposure_time = max_exposure_time,
      seed = current_seed
    )
    
    # extract correlation parameters
    corr_params <- attr(df_synthetic, "correlation_params")
    
    seed_results$seed <- current_seed
    seed_results$r <- corr_params$r
    seed_results$sigma_epsilon <- corr_params$sigma_epsilon
    seed_results$sigma_residual <- corr_params$sigma_residual
    seed_results$status <- "success"
    
    # clean up memory
    rm(df_synthetic)
    
  }, error = function(e) {
    seed_results$seed <- current_seed
    seed_results$r <- NA
    seed_results$sigma_epsilon <- NA
    seed_results$sigma_residual <- NA
    seed_results$status <- "error"
    seed_results$error_message <- as.character(e)
    warning(sprintf("error in simulation %d: %s", sim_iter, e))
  })
  
  # save results for this seed
  filename <- file.path(results_dir, sprintf("corr_params_seed_%04d.rds", current_seed))
  saveRDS(seed_results, file = filename)
  
  gc(verbose = FALSE)
  
  # progress update every 100 simulations
  if (sim_iter %% 100 == 0) {
    cat(sprintf("completed %d simulations (%.1f%%)\n", sim_iter, 100 * sim_iter / n_sims))
  }
}