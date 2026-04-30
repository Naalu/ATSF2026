###############################################################################
# FILE: scripts/run_bsts_validation.R
#
# PURPOSE:
#   Driver script that generates KReger-bsts_seasonal forecasts for the 56
#   validation-season origin dates ONLY (seasons 2015-16 and 2016-17,
#   2015-10-24 through 2017-05-06).
#
#   Why a separate driver instead of letting the engine auto-loop?
#   The engine's run_loop=TRUE path iterates over ALL origin dates in
#   tasks.json (134 dates across 5 seasons). For Phase 1 of the final-
#   submission plan (FLUCAST_IMPLEMENTATION_PLAN.md, Checkpoint 1A) we
#   want validation-only output so the test seasons (3-5) remain unseen
#   for downstream model selection. This driver enforces that constraint
#   without modifying the engine.
#
# INPUTS:
#   - scripts/forecast_engine.R    (the generic engine)
#   - configs/bsts_seasonal.R      (the bsts model config)
#   - hub-config/tasks.json        (origin date enumeration, read by engine)
#   - target-data/time-series.csv  (wILI training data, read by engine)
#
# OUTPUTS:
#   - 56 CSVs in model-output/KReger-bsts_seasonal/
#     (one per validation-season origin date)
#   - Console log of per-date timing
#
# HOW TO RUN:
#   # From the project root in an R session:
#   setwd("<path-to-ATSF2026-repo>")
#   source("scripts/run_bsts_validation.R")
#
#   # OR from the shell:
#   cd "<path-to-ATSF2026-repo>"
#   Rscript scripts/run_bsts_validation.R
#
# RESUMABILITY:
#   The bsts config has overwrite=FALSE, so existing CSVs are skipped on
#   rerun. Safe to Ctrl-C and restart — it'll pick up where it left off.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md §3.1 (season boundaries),
#     Checkpoint 1A (validation-only run)
###############################################################################


# ============================================================================
# CONFIGURATION CONSTANTS
# ============================================================================

# The cutoff date that separates validation (seasons 1-2) from test
# (seasons 3-5). Per FLUCAST_IMPLEMENTATION_PLAN.md §3.1, season 2016-17
# ends 2017-05-06; the next season starts 2017-10-28. Anything <= this
# date is validation; anything > is test.
VALIDATION_END_DATE <- as.Date("2017-05-06")


# ============================================================================
# LOAD ENGINE & CONFIG
# ============================================================================
# Engine first (defines run_forecast and helpers); config second (defines
# `config` and the BSTS_* constants for the custom_simulate_fn closure).
source("scripts/forecast_engine.R")
source("configs/bsts_seasonal.R")


# ============================================================================
# INITIALIZE WITHOUT LOOPING
#
# The bsts config has run_loop = FALSE. run_forecast(config) therefore does
# all the setup work — loads the wILI data, parses origin dates from
# tasks.json, validates the config, prepares the output directory — and
# returns a `result` object with a callable forecast() method, but does NOT
# iterate over dates. We drive the iteration ourselves below.
# ============================================================================
result <- run_forecast(config)


# ============================================================================
# FILTER TO VALIDATION-SEASON DATES
# ============================================================================
# result$origin_dates is the full vector of 134 origin dates from tasks.json.
# We subset to the 56 dates on or before VALIDATION_END_DATE.
all_dates        <- result$origin_dates
validation_dates <- all_dates[all_dates <= VALIDATION_END_DATE]

# Defensive sanity check. The implementation plan documents ~56 dates
# (29 + 27 from the two regular flu seasons), but tasks.json may include
# additional off-season dates between seasons. We don't hardcode the count;
# instead we sanity-check that we got *some* dates and they all fall within
# the expected calendar window.
if (length(validation_dates) == 0L) {
  stop("No validation dates found! Check VALIDATION_END_DATE and tasks.json.")
}
if (length(validation_dates) > length(all_dates)) {
  stop("More validation dates than total origin dates — impossible. Investigate.")
}
# Lower bound: the two flu seasons contribute 29 + 27 = 56 in-season dates.
# Anything less means we're missing in-season dates and shouldn't proceed.
if (length(validation_dates) < 56L) {
  warning(sprintf(
    "Only %d validation dates found; expected at least 56 (29 + 27 in-season). ",
    length(validation_dates)),
    "Investigate before relying on this output."
  )
}

cat(sprintf(
  "Validation set: %d origin dates from %s to %s\n",
  length(validation_dates),
  format(min(validation_dates)),
  format(max(validation_dates))
))


# ============================================================================
# RUNTIME ESTIMATE
# ============================================================================
# bsts at 1100 iterations x 11 locations runs in ~30-90 sec/date on a
# typical laptop. Print a budget so the user knows roughly when to come
# back. (We intentionally don't return early on long runs — this is just
# informational.)
cat(sprintf(
  "Estimated runtime: %.1f - %.1f minutes (varies by data size & hardware)\n",
  length(validation_dates) * 30 / 60,
  length(validation_dates) * 90 / 60
))


# ============================================================================
# MAIN LOOP
# ============================================================================
# For each validation date, call result$forecast(d). This invokes the
# engine's full per-date pipeline:
#   set per-date seed -> filter train_data -> custom_simulate_fn ->
#   compute quantiles -> format for hub -> validate -> write CSV
#
# The wall-clock time per iteration is logged so we can spot pathologically
# slow dates (e.g., late-season dates where bsts's data history is largest).
# ============================================================================

# Track total wall-clock time so we can print a final summary.
loop_start <- Sys.time()

# Counters for the end-of-run summary.
n_processed <- 0L
n_skipped   <- 0L
n_errored   <- 0L

for (i in seq_along(validation_dates)) {
  d <- validation_dates[i]

  cat(sprintf("[%d/%d] %s ... ", i, length(validation_dates), format(d)))

  date_start <- Sys.time()

  # Wrap each date in tryCatch so one bad date doesn't kill the whole run.
  # The engine's overwrite=FALSE means a rerun will resume from where this
  # one stopped (skipping any successfully-written CSVs).
  out <- tryCatch(
    {
      # forecast() returns NULL when it skips an existing CSV (overwrite=FALSE).
      # When it writes a new one, it returns the formatted forecast tibble.
      result$forecast(d)
    },
    error = function(e) {
      # Tag the error with the date so the log makes the failure attributable.
      message(sprintf("\n  ERROR on %s: %s", format(d), conditionMessage(e)))
      structure(list(), class = "bsts_error")
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), date_start, units = "secs"))

  if (inherits(out, "bsts_error")) {
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
cat("  BSTS Validation Run Complete\n")
cat("========================================\n")
cat(sprintf("  Processed: %d new CSVs\n", n_processed))
cat(sprintf("  Skipped:   %d (already existed)\n", n_skipped))
cat(sprintf("  Errored:   %d\n", n_errored))
cat(sprintf("  Total time: %.1f minutes\n", total_elapsed))
cat(sprintf("  Output:    model-output/KReger-bsts_seasonal/\n"))
cat("========================================\n")

# Final CSV count check. We expect exactly 56 files in the output dir
# after a successful run. Doesn't fail — just informational.
output_dir <- file.path(config$hub_path, "model-output",
                        paste0(config$team_abbr, "-", config$model_abbr))
csv_count <- length(list.files(output_dir, pattern = "\\.csv$"))
cat(sprintf("\nCSVs in output dir: %d (expected: %d)\n",
            csv_count, length(validation_dates)))

if (csv_count < length(validation_dates)) {
  cat("\nSome dates didn't produce a CSV. Re-run this script to retry the missing dates.\n")
}
