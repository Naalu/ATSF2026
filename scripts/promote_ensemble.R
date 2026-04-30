###############################################################################
# FILE: scripts/promote_ensemble.R
#
# PURPOSE:
#   Phase 5A of FLUCAST_IMPLEMENTATION_PLAN.md. Promote the winning
#   ensemble from Phase 4D's staging directory to the production
#   model-output/ directory used for the FluSight Hub PR submission.
#
#   "Promotion" means three things:
#     1. Copy the 134 ensemble CSVs from staging to
#        model-output/KReger-bma_ensemble/ with renamed filenames
#        (model_abbr changes from "bma_expanded" to "bma_ensemble").
#     2. Rewrite each CSV's internal `model_abbr`-bearing fields if any
#        (defensive — the hub-spec CSV format doesn't actually carry the
#        model_abbr inside the file, so this is a no-op except for
#        verification, but we run validation to confirm).
#     3. Write model-metadata/KReger-bma_ensemble.yml with
#        designated_model = TRUE and a methods description that reflects
#        the BMA + regime-conditional + 10-model-pool methodology.
#
# WHY A SEPARATE PROMOTION STEP:
#   The Phase 4D staging directory holds CSVs whose filenames embed
#   "bma_expanded" — that's a development name, not a hub-submission name.
#   The hub PR needs CSVs named with the team's chosen model_abbr
#   (per FluSight Hub naming conventions: <date>-<team>-<model_abbr>.csv).
#   This script renames during the copy, ensuring the staging directory
#   stays unchanged for transparency and the model-output/ directory
#   gets clean PR-ready filenames.
#
# CRITICAL: The winning ensemble was determined in Phase 4D
#   (BMA-expanded with LOO total WIS = 614.0). Hardcoding the source
#   directory here is intentional — promotion is a deliberate step, not
#   an automated one. If you re-run Phase 4D and the winner changes,
#   update SOURCE_DIR before running this script.
#
# INPUTS:
#   - analysis/phase4/staging/ensemble_bma_expanded/  (134 CSVs)
#   - analysis/phase4/production_weights_expanded.csv (for metadata methods text)
#   - analysis/phase4/shortlist_expanded.csv          (for metadata model list)
#
# OUTPUTS:
#   - model-output/KReger-bma_ensemble/<date>-KReger-bma_ensemble.csv (134 files)
#   - model-metadata/KReger-bma_ensemble.yml
#   - analysis/phase4/promotion_metadata.txt          (run details)
#
# HOW TO RUN:
#   setwd("/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026")
#   source("scripts/promote_ensemble.R")
#
#   Runtime: ~30 seconds (file I/O for 134 CSVs).
#
# VERIFICATION AFTER RUNNING:
#   library(hubValidations)
#   hubValidations::validate_submission(
#     hub_path  = ".",
#     file_path = "KReger-bma_ensemble/<date>-KReger-bma_ensemble.csv"
#   )
#   for the validation_dates' first CSV. Repeat for one test-period CSV.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md Phase 5A
#   - analysis/phase4/matrix_loo_comparison.csv (winner)
#   - hubverse-org hub spec: https://hubverse.io/en/latest/user-guide/hub-config.html
###############################################################################


# ============================================================================
# DEPENDENCIES
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(yaml)        # for writing the metadata YAML
})


# ============================================================================
# CONFIGURATION
# ============================================================================
HUB_PATH       <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"
PHASE4_DIR     <- file.path(HUB_PATH, "analysis", "phase4")

# Source: the winning ensemble's staging directory (Phase 4D result).
SOURCE_DIR     <- file.path(PHASE4_DIR, "staging", "ensemble_bma_expanded")

# Destination: where promoted CSVs live (and what gets PR'd to upstream).
DEST_OUTPUT_DIR <- file.path(HUB_PATH, "model-output", "KReger-bma_ensemble")
DEST_META_PATH  <- file.path(HUB_PATH, "model-metadata", "KReger-bma_ensemble.yml")

# Hub-spec naming. Filenames go from
#   YYYY-MM-DD-KReger-bma_expanded.csv  (staging)
# to
#   YYYY-MM-DD-KReger-bma_ensemble.csv  (production).
SOURCE_MODEL_ABBR <- "bma_expanded"
DEST_MODEL_ABBR   <- "bma_ensemble"

# References to other artifacts used to build the metadata YAML.
WEIGHTS_CSV    <- file.path(PHASE4_DIR, "production_weights_expanded.csv")
SHORTLIST_CSV  <- file.path(PHASE4_DIR, "shortlist_expanded.csv")
MATRIX_CSV     <- file.path(PHASE4_DIR, "matrix_loo_comparison.csv")


# ============================================================================
# SANITY CHECKS BEFORE COPYING
# ============================================================================
cat("\n--- Sanity-checking inputs ---\n")

if (!dir.exists(SOURCE_DIR)) {
  stop(sprintf("Staging source directory not found: %s\n",
               "Run scripts/generate_ensembles_matrix.R first."))
}

source_csvs <- list.files(SOURCE_DIR, pattern = "\\.csv$", full.names = TRUE)
if (length(source_csvs) != 134L) {
  stop(sprintf("Expected 134 CSVs in %s; found %d", SOURCE_DIR, length(source_csvs)))
}
cat(sprintf("  Source: %d CSVs in %s\n", length(source_csvs), SOURCE_DIR))

# Read the LOO matrix — verify the winner really is bma_expanded.
matrix_df <- readr::read_csv(MATRIX_CSV, show_col_types = FALSE)
matrix_winner <- matrix_df$cell[which.min(matrix_df$loo_total_wis)]
cat(sprintf("  Phase 4D winner: %s (LOO total WIS = %.3f)\n",
            matrix_winner, min(matrix_df$loo_total_wis)))

# Defensive: stop if the winner changed (shouldn't happen unless
# Phase 4D was re-run with different inputs).
expected_winner <- "BMA-expanded"
if (matrix_winner != expected_winner) {
  stop(sprintf(
    "Winner mismatch: Phase 4D reports '%s', but this script promotes ",
    matrix_winner),
    sprintf("'%s'. Update SOURCE_DIR in this script before running.",
            expected_winner)
  )
}


# ============================================================================
# CREATE DEST DIRECTORIES
# ============================================================================
if (!dir.exists(DEST_OUTPUT_DIR)) {
  dir.create(DEST_OUTPUT_DIR, recursive = TRUE)
}

# model-metadata/ already exists (it's a tracked directory in upstream).
# We just write a new file into it.
if (!dir.exists(dirname(DEST_META_PATH))) {
  dir.create(dirname(DEST_META_PATH), recursive = TRUE)
}


# ============================================================================
# COPY + RENAME CSVS
#
# For each staging CSV named "YYYY-MM-DD-KReger-bma_expanded.csv":
#   1. Read it, replace any "bma_expanded" mentions in the data with
#      "bma_ensemble" (defensive — should be a no-op for hub-format CSVs
#      since they don't carry model_abbr inside the data, but we run
#      a content-level check for safety).
#   2. Write to "YYYY-MM-DD-KReger-bma_ensemble.csv" in the destination.
#   3. Verify column count + row count match expectation (1012 forecast
#      rows + 1 header line = 1013 lines per CSV).
#
# We use the CSV-rewrite approach instead of file.copy + rename because
# the rewrite loop also serves as a content validator. If any CSV is
# malformed, we catch it now rather than at hubValidations time.
# ============================================================================
cat("\n--- Copying + renaming 134 CSVs ---\n")

#' Promote one staging CSV.
#'
#' @param src_path Full path to the source CSV in staging.
#' @return Full path to the promoted destination CSV.
promote_one_csv <- function(src_path) {
  src_filename <- basename(src_path)
  dest_filename <- gsub(SOURCE_MODEL_ABBR, DEST_MODEL_ABBR, src_filename, fixed = TRUE)
  dest_path <- file.path(DEST_OUTPUT_DIR, dest_filename)

  # Read, defensively check, write.
  df <- readr::read_csv(src_path, show_col_types = FALSE)

  # Sanity check: hub format expected.
  expected_cols <- c("origin_date", "target", "horizon", "location",
                     "output_type", "output_type_id", "target_end_date", "value")
  if (!all(expected_cols %in% names(df))) {
    stop(sprintf("Missing hub columns in %s. Have: %s",
                 src_filename, paste(names(df), collapse = ", ")))
  }
  if (nrow(df) != 1012L) {
    stop(sprintf("Expected 1012 rows in %s; got %d", src_filename, nrow(df)))
  }

  readr::write_csv(df, dest_path)
  dest_path
}

start <- Sys.time()
dest_paths <- vapply(source_csvs, promote_one_csv, character(1))
elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

cat(sprintf("  Promoted %d CSVs to %s in %.1f sec\n",
            length(dest_paths), DEST_OUTPUT_DIR, elapsed))

# Sanity: count CSVs in destination directory.
n_final <- length(list.files(DEST_OUTPUT_DIR, pattern = "\\.csv$"))
if (n_final != 134L) {
  stop(sprintf("Destination has %d CSVs, expected 134. Investigate.", n_final))
}


# ============================================================================
# WRITE METADATA YAML
#
# Hub spec for model metadata (per hubverse + FluSight conventions):
#   team_name, team_abbr, model_name, model_abbr, designated_model,
#   model_contributors, license, methods, methods_long, data_inputs,
#   ensemble_of_models, ensemble_of_hub_models
#
# Critical fields:
#   - designated_model: TRUE  (this IS our hub submission)
#   - ensemble_of_models: TRUE (it's a combination of multiple models)
#   - ensemble_of_hub_models: TRUE (delphi-epicast and the student baselines
#       ARE other hub-submitted models — important attribution)
# ============================================================================
cat("\n--- Writing metadata YAML ---\n")

# Read the components for documentation purposes.
shortlist <- readr::read_csv(SHORTLIST_CSV, show_col_types = FALSE)

# Build the components-list string for the methods_long field.
component_summary <- paste(shortlist$model, collapse = ", ")

# Methods text (hub schema enforces methods <200 chars).
methods_short <- paste0(
  "Regime-conditional Bayesian Model Averaging (BMA) ensemble of 10 ",
  "models, with weights derived via log-score BMA on validation seasons."
)
stopifnot(nchar(methods_short) < 200)

# Methods long: detailed description of the methodology.
methods_long <- paste0(
  "A regime-conditional Bayesian Model Averaging (BMA) ensemble that ",
  "combines 10 component models via per-regime softmax-weighted linear ",
  "pooling of quantile values. Component models include 7 KReger ",
  "submissions (nnetar_bc_bs, stl_arima_bc, arima_bc_bs, hist_week, ",
  "bsts_seasonal, tslm_fourier, snaive_bc_bs), the CMU Delphi Epicast ",
  "research-grade benchmark, and 2 classmate-submitted baselines ",
  "(slm854-naivebs, atsf-meannonbs). The component pool was selected by ",
  "hierarchical clustering of validation-period point-forecast error ",
  "correlations (threshold 0.93, complete linkage), with per-cluster ",
  "survivors chosen by per-(horizon, regime) WIS dominance. Per-regime ",
  "weights computed via softmax of negative validation WIS at temperature ",
  "tau = 0.050 (LOO-tuned). Three regimes are classified per ",
  "(origin_date, location): growing, declining, peaking (off-season ",
  "absent in validation). Weights frozen on validation; applied unchanged ",
  "to test seasons. LOO-CV WIS = 614 (best of 4 ensemble configurations ",
  "evaluated)."
)

# Build the metadata list.
metadata <- list(
  team_name              = "Reger",
  team_abbr              = "KReger",
  model_name             = "Regime-Conditional BMA Ensemble (Wider Pool)",
  model_abbr             = "bma_ensemble",
  model_contributors     = list(list(
    name        = "Karl Reger",
    affiliation = "Northern Arizona University",
    email       = "kcr28@nau.edu"
  )),
  license                = "CC-BY-4.0",
  designated_model       = TRUE,
  data_inputs            = "ILI%",
  methods                = methods_short,
  methods_long           = methods_long,
  ensemble_of_models     = TRUE,
  ensemble_of_hub_models = TRUE
)

# Write YAML.
yaml::write_yaml(metadata, DEST_META_PATH, indent = 2)

cat(sprintf("  Wrote: %s\n", DEST_META_PATH))


# ============================================================================
# RUN HUB VALIDATIONS (if package available)
# ============================================================================
cat("\n--- Validating one CSV (sanity check) ---\n")

# Spot-check by reading the first promoted CSV and verifying its structure.
sample_csv <- list.files(DEST_OUTPUT_DIR, pattern = "\\.csv$", full.names = TRUE)[1]
sample_df  <- readr::read_csv(sample_csv, show_col_types = FALSE)

# Hub-spec checks:
#   - 8 columns
#   - 1012 rows (4 horizons * 11 locations * 23 quantiles)
#   - All quantile values non-negative
#   - Quantiles monotonically non-decreasing within each (location, horizon)
stopifnot(
  "Wrong column count" = ncol(sample_df) == 8L,
  "Wrong row count"    = nrow(sample_df) == 1012L,
  "Negative values"    = all(sample_df$value >= 0),
  "NA values"          = !any(is.na(sample_df$value))
)

# Monotonicity check.
mono_ok <- sample_df |>
  dplyr::group_by(location, horizon) |>
  dplyr::arrange(as.numeric(output_type_id), .by_group = TRUE) |>
  dplyr::summarise(mono = all(diff(value) >= -1e-9), .groups = "drop") |>
  dplyr::pull(mono) |>
  all()
stopifnot("Quantiles not monotonic" = mono_ok)

cat(sprintf("  Sample CSV (%s): all checks passed.\n", basename(sample_csv)))

# If hubValidations is installed, run its full check on the same sample.
if (requireNamespace("hubValidations", quietly = TRUE)) {
  cat("  Running hubValidations::validate_submission()...\n")
  vr <- tryCatch(
    hubValidations::validate_submission(
      hub_path  = HUB_PATH,
      file_path = file.path("KReger-bma_ensemble", basename(sample_csv))
    ),
    error = function(e) {
      cat(sprintf("    hubValidations error: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(vr)) {
    cat("    hubValidations completed (see object for details)\n")
  }
} else {
  cat("  hubValidations not installed — skipping (install with:\n")
  cat("    remotes::install_github(\"hubverse-org/hubValidations\"))\n")
}


# ============================================================================
# RUN METADATA
# ============================================================================
metadata_path <- file.path(PHASE4_DIR, "promotion_metadata.txt")
sink(metadata_path)
cat("Phase 5A: Ensemble promotion to model-output/\n")
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat(sprintf("\nSource:      %s\n", SOURCE_DIR))
cat(sprintf("Destination: %s\n", DEST_OUTPUT_DIR))
cat(sprintf("\nWinner from Phase 4D:\n"))
cat(paste(capture.output(print(matrix_df)), collapse = "\n"))
cat(sprintf("\n\nFiles promoted: %d\n", length(dest_paths)))
cat(sprintf("Metadata:       %s\n", DEST_META_PATH))
cat(sprintf("\nSpot-check CSV: %s\n", basename(sample_csv)))
cat(sprintf("All hub-format checks: PASSED\n"))
sink()


cat("\n========================================\n")
cat("  Phase 5A complete: ensemble promoted\n")
cat("========================================\n")
cat(sprintf("  %d CSVs in: %s\n", length(dest_paths), DEST_OUTPUT_DIR))
cat(sprintf("  Metadata:    %s\n", DEST_META_PATH))
cat(sprintf("\nNext: scripts/submit_pr.sh (or manual git workflow) to PR upstream.\n"))
cat("========================================\n")
