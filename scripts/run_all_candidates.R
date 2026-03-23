###############################################################################
# FILE: scripts/run_all_candidates.R
#
# Batch runner for all candidate ensemble component models. Sources each
# config file and runs the full forecast loop. Designed to be run from
# the ATSF2026 repo root directory.
#
# USAGE:
#   1. Open R in the ATSF2026 repo root directory
#   2. source("scripts/run_all_candidates.R")
#
# Or run a single model interactively:
#   source("scripts/forecast_engine.R")
#   source("configs/arima_bc_bs.R")
#   config$run_loop <- TRUE
#   result <- run_forecast(config)
#
# RUNTIME ESTIMATES (per model, 139 origin dates, n_sims=1000):
#   snaive_bc_bs  ~  5-10 min  (fast: no parameter estimation)
#   ets_log       ~ 15-30 min  (moderate: auto-selects ETS variant)
#   arima_bc_bs   ~ 30-60 min  (slow: searches ARIMA order space)
#   nnetar_bc_bs  ~ 45-90 min  (slow: fits 20 networks per location)
#   theta_bc_bs   ~ 10-20 min  (moderate: SES + seasonal decomposition)
#   tslm_fourier  ~  5-10 min  (fast: fixed regression, no iteration)
#   stl_arima_bc  ~ 30-60 min  (slow: STL decomposition + ARIMA search)
#   hist_week     ~  3-5  min  (fastest: pure resampling, no fitting)
#
# TOTAL: roughly 2.5 - 5 hours for all 8 models.
#
# TIPS FOR MANAGING RUNTIME:
#   - Set n_sims = 500 for the initial exploratory run (halves sim time)
#   - Set dev_mode = TRUE to process only the first 5 dates (quick test)
#   - Run multiple R sessions in parallel (each sourcing a different config)
#   - The skip-if-exists logic means you can interrupt and resume safely
#
# PARALLEL EXECUTION (recommended):
#   Open 3-4 terminal windows, each running R. In each, run:
#     source("scripts/forecast_engine.R")
#     source("configs/<model>.R")
#     config$run_loop <- TRUE
#     result <- run_forecast(config)
#   Assign the slowest models (arima, nnetar, stl_arima) to separate
#   sessions, and batch the fast ones together.
#
###############################################################################

cat("==========================================================\n")
cat("  BATCH RUNNER: All Candidate Ensemble Models\n")
cat("==========================================================\n\n")

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

# Set to TRUE for a quick test (5 dates, 500 sims per model)
QUICK_TEST <- FALSE

# List of all config files to run, in order of expected runtime (fast first)
# so you get early results quickly.
config_files <- c(
  "configs/hist_week.R",       # ~3 min   (custom: pure resampling)
  "configs/snaive_bc_bs.R",    # ~5 min   (no parameter estimation)
  "configs/tslm_fourier.R",    # ~5 min   (fixed regression)
  "configs/theta_bc_bs.R",     # ~15 min  (SES + decomposition)
  "configs/ets_log.R",         # ~20 min  (auto ETS selection)
  "configs/arima_bc_bs.R",     # ~45 min  (auto ARIMA search)
  "configs/stl_arima_bc.R",    # ~45 min  (STL + auto ARIMA)
  "configs/nnetar_bc_bs.R"     # ~60 min  (20 neural networks)
)

# ---------------------------------------------------------------------------
# SOURCE THE ENGINE (once)
# ---------------------------------------------------------------------------

source("scripts/forecast_engine.R")

# ---------------------------------------------------------------------------
# MAIN LOOP
# ---------------------------------------------------------------------------

batch_start <- Sys.time()
results <- list()

for (cfg_path in config_files) {

  cat(sprintf("\n>>> Loading config: %s\n", cfg_path))

  # Check that the config file exists

  if (!file.exists(cfg_path)) {
    cat(sprintf("    SKIPPED: file not found at %s\n", cfg_path))
    next
  }

  # Source the config (creates a `config` object in global env)
  source(cfg_path)

  # Override run control for batch execution
  config$run_loop <- TRUE

  # Apply quick-test overrides if requested
  if (QUICK_TEST) {
    config$dev_mode <- TRUE       # only first 5 dates
    config$n_sims   <- 500        # fewer simulation paths
    cat("    [QUICK_TEST mode: 5 dates, 500 sims]\n")
  }

  # Run the forecast pipeline
  model_id <- paste0(config$team_abbr, "-", config$model_abbr)
  cat(sprintf(">>> Running: %s\n", model_id))

  tryCatch({
    result <- run_forecast(config)
    results[[model_id]] <- result
    cat(sprintf(">>> DONE: %s\n", model_id))
  }, error = function(e) {
    cat(sprintf(">>> FAILED: %s -- %s\n", model_id, e$message))
  })
}

# ---------------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------------

batch_elapsed <- as.numeric(difftime(Sys.time(), batch_start, units = "mins"))

cat("\n==========================================================\n")
cat("  BATCH COMPLETE\n")
cat("==========================================================\n\n")
cat(sprintf("  Total time: %.1f minutes\n", batch_elapsed))
cat(sprintf("  Models run: %d / %d\n", length(results), length(config_files)))

# Count output CSVs per model
cat("\n  Output file counts:\n")
for (model_id in names(results)) {
  output_dir <- results[[model_id]]$output_dir
  if (dir.exists(output_dir)) {
    n_files <- length(list.files(output_dir, pattern = "\\.csv$"))
    cat(sprintf("    %s: %d files\n", model_id, n_files))
  }
}

cat("\n  Next step: run scripts/analyze_candidates.R for error correlation analysis\n")
