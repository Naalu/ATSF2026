###############################################################################
# FILE: scripts/run_validation.R
#
# PURPOSE:
#   Generic driver for running ANY config's validation forecasts (seasons
#   1-2, 57 origin dates from 2015-10-24 through 2017-05-06). Replaces the
#   need for one-off run_<model>_validation.R scripts.
#
#   Phase 1 of FLUCAST_IMPLEMENTATION_PLAN.md introduces multiple new
#   candidate models. Rather than copy-pasting run_bsts_validation.R three
#   times (once each for arima_log, ets_bc, nnetar_log), this generic
#   driver takes the config path as a command-line argument or as a
#   variable set before sourcing.
#
#   The bsts driver (scripts/run_bsts_validation.R) predates this and
#   stays intact for backwards compatibility / git history clarity.
#
# INPUTS:
#   CONFIG_PATH — Path to the model config to run. Two ways to provide:
#     (a) Set the variable before sourcing:
#         CONFIG_PATH <- "configs/arima_log.R"
#         source("scripts/run_validation.R")
#     (b) Pass as a Rscript command-line argument:
#         Rscript scripts/run_validation.R configs/arima_log.R
#
# OUTPUTS:
#   - Up to 57 CSVs in model-output/<team>-<model>/
#   - Per-date timing log to console
#   - End-of-run summary
#
# RESUMABILITY:
#   Each config has overwrite=FALSE, so existing CSVs are skipped on rerun.
#   Safe to Ctrl-C and restart — picks up where it left off.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md §3.1 (validation window),
#     Checkpoint 1C (transform variants)
#   - scripts/run_bsts_validation.R (template this is generalized from)
###############################################################################


# ============================================================================
# RESOLVE CONFIG_PATH
# ============================================================================
# Prefer (a) a pre-set R variable; fall back to (b) commandArgs() for
# Rscript invocations. If neither is set, error with usage instructions.
if (!exists("CONFIG_PATH")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1L) {
    CONFIG_PATH <- args[[1]]
  } else {
    stop(
      "CONFIG_PATH not set. Provide one of:\n",
      "  (a) Before sourcing in R:\n",
      "      CONFIG_PATH <- \"configs/arima_log.R\"; source(\"scripts/run_validation.R\")\n",
      "  (b) From shell:\n",
      "      Rscript scripts/run_validation.R configs/arima_log.R"
    )
  }
}

if (!file.exists(CONFIG_PATH)) {
  stop(sprintf("Config not found: %s", CONFIG_PATH))
}


# ============================================================================
# CONFIGURATION CONSTANTS
# ============================================================================

# Validation/test split — same value bsts driver uses. Last validation
# origin date in tasks.json (last Saturday of Season 2016-17).
VALIDATION_END_DATE <- as.Date("2017-05-06")


# ============================================================================
# LOAD ENGINE & CONFIG
# ============================================================================
source("scripts/forecast_engine.R")
source(CONFIG_PATH)

# Defensive: confirm sourcing CONFIG_PATH actually populated `config`.
# Catches the case where someone names a stale or malformed file.
if (!exists("config") || !is.list(config)) {
  stop(sprintf(
    "Sourcing %s did not produce a valid `config` list.", CONFIG_PATH
  ))
}

cat(sprintf(
  "\n=== Validation run: %s-%s ===\n",
  config$team_abbr, config$model_abbr
))
cat(sprintf("  Config: %s\n", CONFIG_PATH))


# ============================================================================
# INITIALIZE WITHOUT LOOPING
# ============================================================================
result <- run_forecast(config)


# ============================================================================
# FILTER TO VALIDATION-SEASON DATES
# ============================================================================
all_dates        <- result$origin_dates
validation_dates <- all_dates[all_dates <= VALIDATION_END_DATE]

# Same data-derived sanity checks as bsts driver. We don't hardcode the
# expected count (57) because tasks.json defines it; we only enforce a
# defensive lower bound (56 = 29 + 27 in-season weeks across 2 seasons).
if (length(validation_dates) == 0L) {
  stop("No validation dates found! Check VALIDATION_END_DATE and tasks.json.")
}
if (length(validation_dates) > length(all_dates)) {
  stop("More validation dates than total origin dates — impossible. Investigate.")
}
if (length(validation_dates) < 56L) {
  warning(sprintf(
    "Only %d validation dates found; expected at least 56 (29 + 27 in-season). ",
    length(validation_dates)),
    "Investigate before relying on this output."
  )
}

cat(sprintf(
  "  Validation set: %d origin dates from %s to %s\n",
  length(validation_dates),
  format(min(validation_dates)),
  format(max(validation_dates))
))


# ============================================================================
# MAIN LOOP
# ============================================================================
loop_start <- Sys.time()

n_processed <- 0L
n_skipped   <- 0L
n_errored   <- 0L

for (i in seq_along(validation_dates)) {
  d <- validation_dates[i]

  cat(sprintf("[%d/%d] %s ... ", i, length(validation_dates), format(d)))

  date_start <- Sys.time()

  out <- tryCatch(
    {
      result$forecast(d)
    },
    error = function(e) {
      message(sprintf("\n  ERROR on %s: %s", format(d), conditionMessage(e)))
      structure(list(), class = "validation_error")
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), date_start, units = "secs"))

  if (inherits(out, "validation_error")) {
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
cat(sprintf("  Validation Run Complete: %s-%s\n",
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
cat(sprintf("\nCSVs in output dir: %d (expected: %d)\n",
            csv_count, length(validation_dates)))

if (csv_count < length(validation_dates)) {
  cat("\nSome dates didn't produce a CSV. Re-run this script to retry the missing dates.\n")
}
