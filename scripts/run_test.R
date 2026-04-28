###############################################################################
# FILE: scripts/run_test.R
#
# PURPOSE:
#   Generic driver for running a config's TEST-period forecasts (Seasons 3-5,
#   2017-10-28 through 2020-02-29, 77 origin dates). Sibling to
#   scripts/run_validation.R, which handles the validation period.
#
#   Used during Phase 3 to fill gaps in ensemble-component coverage:
#   ets_bc and bsts_seasonal were generated for validation only during
#   Phase 1 candidate selection. To build the full 134-date ensemble in
#   Phase 3B, we need both models' test-period forecasts too.
#
#   USE THIS SCRIPT ONLY for component models that have already passed
#   Phase 2 selection. Running test-period forecasts for any candidate
#   that's not in the final ensemble is wasted compute and clutters the
#   model-output directory.
#
# INPUTS:
#   CONFIG_PATH — Path to the model config to run. Either:
#     (a) Set the variable before sourcing:
#         CONFIG_PATH <- "configs/ets_bc.R"
#         source("scripts/run_test.R")
#     (b) Pass as a Rscript argument:
#         Rscript scripts/run_test.R configs/ets_bc.R
#
# OUTPUTS:
#   - Up to 77 CSVs in model-output/<team>-<model>/
#   - Per-date timing log to console
#   - End-of-run summary
#
# RESUMABILITY:
#   Each config has overwrite=FALSE, so existing CSVs are skipped on rerun.
#   Existing validation CSVs (57 per model) will be skipped automatically
#   since they exist; only the 77 missing test-period CSVs will be written.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md §3.1 (test window),
#     Phase 3B (ensemble forecast generation, requires full 134 dates)
#   - scripts/run_validation.R (sibling for validation period)
###############################################################################


# ============================================================================
# RESOLVE CONFIG_PATH
# ============================================================================
if (!exists("CONFIG_PATH")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1L) {
    CONFIG_PATH <- args[[1]]
  } else {
    stop(
      "CONFIG_PATH not set. Provide one of:\n",
      "  (a) Before sourcing in R:\n",
      "      CONFIG_PATH <- \"configs/ets_bc.R\"; source(\"scripts/run_test.R\")\n",
      "  (b) From shell:\n",
      "      Rscript scripts/run_test.R configs/ets_bc.R"
    )
  }
}

if (!file.exists(CONFIG_PATH)) {
  stop(sprintf("Config not found: %s", CONFIG_PATH))
}


# ============================================================================
# CONFIGURATION CONSTANTS
# ============================================================================

# Validation/test split. Anything strictly AFTER this date is test.
# Same constant as scripts/run_validation.R — kept in sync manually.
VALIDATION_END_DATE <- as.Date("2017-05-06")


# ============================================================================
# LOAD ENGINE & CONFIG
# ============================================================================
source("scripts/forecast_engine.R")
source(CONFIG_PATH)

if (!exists("config") || !is.list(config)) {
  stop(sprintf(
    "Sourcing %s did not produce a valid `config` list.", CONFIG_PATH
  ))
}

cat(sprintf(
  "\n=== Test-period run: %s-%s ===\n",
  config$team_abbr, config$model_abbr
))
cat(sprintf("  Config: %s\n", CONFIG_PATH))


# ============================================================================
# INITIALIZE WITHOUT LOOPING
# ============================================================================
result <- run_forecast(config)


# ============================================================================
# FILTER TO TEST-PERIOD DATES
# ============================================================================
all_dates  <- result$origin_dates
test_dates <- all_dates[all_dates > VALIDATION_END_DATE]

# Same data-derived sanity checks as run_validation.R, just inverted to
# the test-period window. Nominal expected count is 77 (28 + 29 + 20
# in-season dates across Seasons 3-5).
if (length(test_dates) == 0L) {
  stop("No test dates found! Check VALIDATION_END_DATE and tasks.json.")
}
if (length(test_dates) > length(all_dates)) {
  stop("More test dates than total origin dates - impossible. Investigate.")
}
if (length(test_dates) < 77L) {
  warning(sprintf(
    "Only %d test dates found; expected at least 77 (28 + 29 + 20 in-season). ",
    length(test_dates)),
    "Investigate before relying on this output."
  )
}

cat(sprintf(
  "  Test set: %d origin dates from %s to %s\n",
  length(test_dates),
  format(min(test_dates)),
  format(max(test_dates))
))


# ============================================================================
# MAIN LOOP
# ============================================================================
loop_start <- Sys.time()

n_processed <- 0L
n_skipped   <- 0L
n_errored   <- 0L

for (i in seq_along(test_dates)) {
  d <- test_dates[i]

  cat(sprintf("[%d/%d] %s ... ", i, length(test_dates), format(d)))

  date_start <- Sys.time()

  out <- tryCatch(
    {
      result$forecast(d)
    },
    error = function(e) {
      message(sprintf("\n  ERROR on %s: %s", format(d), conditionMessage(e)))
      structure(list(), class = "test_error")
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), date_start, units = "secs"))

  if (inherits(out, "test_error")) {
    n_errored <- n_errored + 1L
    cat(sprintf("FAILED in %.1f sec\n", elapsed))
  } else if (is.null(out)) {
    n_skipped <- n_skipped + 1L
    cat(sprintf("skipped (CSV exists) in %.1f sec\n", elapsed))
  } else {
    n_processed <- n_processed + 1L
    cat(sprintf("done in %.1f sec\n", elapsed))
  }
}


# ============================================================================
# END-OF-RUN SUMMARY
# ============================================================================
total_elapsed <- as.numeric(difftime(Sys.time(), loop_start, units = "mins"))

cat("\n")
cat("========================================\n")
cat(sprintf("  Test Run Complete: %s-%s\n",
            config$team_abbr, config$model_abbr))
cat("========================================\n")
cat(sprintf("  Processed: %d new CSVs\n", n_processed))
cat(sprintf("  Skipped:   %d (already existed)\n", n_skipped))
cat(sprintf("  Errored:   %d\n", n_errored))
cat(sprintf("  Total time: %.1f minutes\n", total_elapsed))
cat("========================================\n")

output_dir <- file.path(config$hub_path, "model-output",
                        paste0(config$team_abbr, "-", config$model_abbr))
csv_count <- length(list.files(output_dir, pattern = "\\.csv$"))
cat(sprintf("\nCSVs in output dir: %d\n", csv_count))

# Expectation: validation CSVs (57) should already exist + test CSVs (77)
# = 134 total. Inform if we don't have full coverage.
if (csv_count < 134L) {
  cat(sprintf("\nNote: only %d of expected 134 total CSVs present. ", csv_count))
  cat("Re-run this script (or run_validation.R) to fill gaps.\n")
} else if (csv_count == 134L) {
  cat("\nFull 134-date coverage achieved.\n")
}
