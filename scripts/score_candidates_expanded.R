###############################################################################
# FILE: scripts/score_candidates_expanded.R
#
# PURPOSE:
#   Phase 4A of FLUCAST_IMPLEMENTATION_PLAN.md. Score the EXPANDED candidate
#   pool on validation data only and compute pairwise correlation matrices.
#
#   The expanded pool extends the Phase 2 pool (11 KReger models) by adding
#   9 non-KReger models pulled from upstream sjfox/ATSF2026: 8 student
#   baseline submissions plus delphi-epicast (CMU's research-grade FluSight
#   benchmark). Total: 11 + 9 = 20 models, but we exclude KReger-theta_bc_bs
#   (not in original analysis pool) and hist-avg (only 132/134 CSVs), so
#   the actual scored set is 17 models.
#
# DIFFERS FROM scripts/score_candidates.R:
#   - Candidate list parameterized at top of script (not hardcoded inside
#     a CANDIDATES vector with hardcoded prefixes)
#   - Loader handles non-KReger directory naming (e.g., "delphi-epicast"
#     not "<team>-<model>") and slight column-format variability
#   - Outputs go to analysis/phase4/, not analysis/phase2/
#   - Critical: still validation-only (57 dates) — same data hygiene as
#     Phase 2; no test-set contamination at the selection stage
#
# WHY THIS RUNS ALONGSIDE Phase 2 (not replacing it):
#   The user has chosen a "parallel track" experiment design — keep the
#   Phase 2/3 KReger-only ensemble as the "primary" submission, and build
#   the wider-pool ensemble in Phase 4 as a secondary experiment. Both
#   ensembles are evaluated on test; whichever performs better will be
#   designated. See analysis/phase2/decisions.md for full rationale.
#
# INPUTS:
#   - target-data/time-series.csv             (ground truth)
#   - hub-config/tasks.json                   (origin date enumeration)
#   - model-output/<model_dir>/*.csv          (17 candidate models)
#
# OUTPUTS in analysis/phase4/:
#   - scores_long_expanded.csv                 (17 * 11 * 4 * 57 = 42,636 rows)
#   - per_model_summary_expanded.csv           (17 rows, sorted by WIS)
#   - correlation_point_error_expanded.csv     (17x17 matrix)
#   - correlation_wis_expanded.csv             (17x17 matrix)
#   - scores_metadata_expanded.txt             (run timestamp + versions)
#
# HOW TO RUN:
#   setwd("/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026")
#   source("scripts/score_candidates_expanded.R")
#
#   Runtime: ~3-4 minutes. Bottleneck is reading 17 * 57 = 969 CSVs.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md Phase 4
#   - analysis/phase2/decisions.md §2.5 (wider-pool design)
#   - scripts/score_candidates.R (template; same scoringutils workflow)
###############################################################################


# ============================================================================
# DEPENDENCIES
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(scoringutils)
})


# ============================================================================
# CONFIGURATION
# ============================================================================
HUB_PATH       <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"
TARGET_DATA    <- file.path(HUB_PATH, "target-data", "time-series.csv")
MODEL_OUTPUT   <- file.path(HUB_PATH, "model-output")
PHASE4_DIR     <- file.path(HUB_PATH, "analysis", "phase4")

# Validation/test cutoff — same constant used everywhere.
VALIDATION_END_DATE <- as.Date("2017-05-06")

# CANDIDATE POOL (17 models)
#
# Each entry is the directory name under model-output/. The "model_label"
# column is what we'll use in tables and correlation matrices — short
# enough to fit in printed output, distinct enough to identify the source.
#
# We separate KReger from non-KReger entries because (a) it makes the
# table easier to read and (b) the loader uses the directory name as-is
# (not prefixed). Both are loaded the same way.
KREGER_CANDIDATES <- c(
  "KReger-snaive_bc_bs",
  "KReger-ets_log",
  "KReger-arima_bc_bs",
  "KReger-nnetar_bc_bs",
  "KReger-tslm_fourier",
  "KReger-stl_arima_bc",
  "KReger-hist_week",
  "KReger-bsts_seasonal",
  "KReger-ets_bc",
  "KReger-arima_log",
  "KReger-nnetar_log"
)
# Note: KReger-theta_bc_bs deliberately excluded (per user decision —
# it was not part of the Phase 2 analysis pool, so including it now
# would be inconsistent with the documented methodology).

UPSTREAM_CANDIDATES <- c(
  "adk-rwdriftbs",       # Random walk + drift + bootstrap
  "ata266-snaivenonbs",  # Seasonal naive, no bootstrap
  "atsf-meanbs",         # Mean + bootstrap
  "atsf-meannonbs",      # Mean, no bootstrap
  "delphi-epicast",      # CMU Delphi research-grade benchmark
  "efm-NAIVEnobs",       # Naive, no bootstrap
  "hcr64-rw_driftNOBS",  # Random walk + drift, no bootstrap
  "slm854-naivebs",      # Naive + bootstrap
  "syed-snaivebs"        # Seasonal naive + bootstrap
)
# Note: hist-avg deliberately excluded (only 132/134 CSVs; missing dates
# would either complicate the analysis or require imputation. Per user
# decision: clean exclusion is preferred over partial-coverage inclusion).

ALL_CANDIDATES <- c(KREGER_CANDIDATES, UPSTREAM_CANDIDATES)
N_MODELS       <- length(ALL_CANDIDATES)

# Shorter labels for tables and plots. Maps directory names to compact
# identifiers. Keep these stable — downstream Phase 4 scripts will
# reference them.
MODEL_LABELS <- c(
  "KReger-snaive_bc_bs"  = "KR_snaive_bc",
  "KReger-ets_log"       = "KR_ets_log",
  "KReger-arima_bc_bs"   = "KR_arima_bc",
  "KReger-nnetar_bc_bs"  = "KR_nnetar_bc",
  "KReger-tslm_fourier"  = "KR_tslm",
  "KReger-stl_arima_bc"  = "KR_stl_arima",
  "KReger-hist_week"     = "KR_hist_wk",
  "KReger-bsts_seasonal" = "KR_bsts",
  "KReger-ets_bc"        = "KR_ets_bc",
  "KReger-arima_log"     = "KR_arima_log",
  "KReger-nnetar_log"    = "KR_nnetar_log",
  "adk-rwdriftbs"        = "adk_rwd",
  "ata266-snaivenonbs"   = "ata_snv",
  "atsf-meanbs"          = "atsf_mean",
  "atsf-meannonbs"       = "atsf_meanNB",
  "delphi-epicast"       = "delphi_epi",
  "efm-NAIVEnobs"        = "efm_naive",
  "hcr64-rw_driftNOBS"   = "hcr_rwd",
  "slm854-naivebs"       = "slm_naive",
  "syed-snaivebs"        = "syed_snv"
)

# Hub-spec quantile thresholds for coverage flags.
COV_50_LOWER <- 0.25
COV_50_UPPER <- 0.75
COV_95_LOWER <- 0.025
COV_95_UPPER <- 0.975

set.seed(20260428L)


# ============================================================================
# DIRECTORY SETUP
# ============================================================================
if (!dir.exists(PHASE4_DIR)) dir.create(PHASE4_DIR, recursive = TRUE)
cat(sprintf("Output directory: %s\n", PHASE4_DIR))


# ============================================================================
# LOAD GROUND TRUTH
# ============================================================================
cat("\n--- Loading ground truth ---\n")

truth <- readr::read_csv(TARGET_DATA, show_col_types = FALSE) |>
  dplyr::select(location, target_end_date, observed = observation)

cat(sprintf("  %d ground-truth rows\n", nrow(truth)))


# ============================================================================
# LOAD CANDIDATE FORECASTS
#
# Each candidate's CSVs live in model-output/<dir_name>/. Hub-format CSVs
# all have the same 8 columns:
#   origin_date, target, horizon, location, output_type, output_type_id,
#   target_end_date, value
#
# Non-KReger models may have minor format variations (date types, integer
# vs character horizon, etc.) — we coerce types defensively.
# ============================================================================
cat("\n--- Loading 17 candidate forecasts (validation-period only) ---\n")

#' Load all validation-period forecasts for one model directory.
#'
#' @param model_dir Character. The directory name under model-output/.
#' @return Tibble with columns matching scoringutils 2.x conventions:
#'   model (the dir name), location, origin_date, horizon,
#'   target_end_date, quantile_level, predicted.
load_model_forecasts <- function(model_dir) {
  dir_path <- file.path(MODEL_OUTPUT, model_dir)
  if (!dir.exists(dir_path)) {
    stop(sprintf("Model directory not found: %s", dir_path))
  }

  csvs <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(csvs) == 0L) {
    stop(sprintf("No CSVs in %s", dir_path))
  }

  # Read with permissive column types — non-KReger CSVs may have horizon
  # as character or origin_date with quotes around it. col_types = cols()
  # uses readr's heuristics for each file independently.
  all_csvs <- purrr::map_dfr(
    csvs,
    readr::read_csv,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_guess())
  )

  # Validation-window filter + type coercion to canonical form.
  all_csvs |>
    dplyr::filter(as.Date(origin_date) <= VALIDATION_END_DATE) |>
    dplyr::transmute(
      model           = model_dir,
      location        = as.character(location),
      origin_date     = as.Date(origin_date),
      horizon         = as.integer(horizon),
      target_end_date = as.Date(target_end_date),
      # output_type_id might come in as character ("0.025") or numeric
      # depending on how the upstream model wrote its CSVs. Force numeric.
      quantile_level  = as.numeric(output_type_id),
      predicted       = as.numeric(value)
    )
}

# Load each candidate and report row counts.
forecast_list <- list()
for (m in ALL_CANDIDATES) {
  fc <- load_model_forecasts(m)
  forecast_list[[m]] <- fc
  cat(sprintf("  %-25s  %5d rows  %d origin dates\n",
              m, nrow(fc), length(unique(fc$origin_date))))
}

forecasts <- dplyr::bind_rows(forecast_list)
cat(sprintf("\n  Total forecast rows: %d\n", nrow(forecasts)))


# ============================================================================
# COVERAGE SANITY: each model should have exactly 57 validation dates.
#
# Models with fewer dates indicate something wrong upstream — either the
# model didn't actually produce a forecast for some validation date, or
# the date filter caught something unexpected.
# ============================================================================
date_counts <- forecasts |>
  dplyr::group_by(model) |>
  dplyr::summarise(n_dates = dplyr::n_distinct(origin_date), .groups = "drop")

if (any(date_counts$n_dates != 57L)) {
  cat("\n  WARNING: not all candidates have 57 validation dates:\n")
  print(date_counts |> dplyr::filter(n_dates != 57L))
  stop("Validation-period origin date count mismatch — investigate before proceeding.")
}


# ============================================================================
# JOIN GROUND TRUTH
# ============================================================================
cat("\n--- Joining ground truth ---\n")

forecasts_scored <- forecasts |>
  dplyr::inner_join(truth, by = c("location", "target_end_date"))

dropped <- nrow(forecasts) - nrow(forecasts_scored)
if (dropped > 0L) {
  warning(sprintf("Dropped %d rows during truth join. Investigate.", dropped))
}
cat(sprintf("  %d rows after join (%d dropped)\n",
            nrow(forecasts_scored), dropped))


# ============================================================================
# SCORE WITH scoringutils
# ============================================================================
cat("\n--- Scoring via scoringutils ---\n")

forecast_obj <- forecasts_scored |>
  scoringutils::as_forecast_quantile(
    forecast_unit = c("model", "location", "origin_date",
                      "target_end_date", "horizon")
  )

scores <- scoringutils::score(forecast_obj) |>
  tibble::as_tibble()

cat(sprintf("  Scored %d forecast units\n", nrow(scores)))


# ============================================================================
# COMPUTE 95% COVERAGE FLAG MANUALLY
# ============================================================================
cov_95 <- forecasts_scored |>
  dplyr::filter(quantile_level %in% c(COV_95_LOWER, COV_95_UPPER)) |>
  tidyr::pivot_wider(
    id_cols      = c(model, location, origin_date, target_end_date,
                     horizon, observed),
    names_from   = quantile_level,
    values_from  = predicted,
    names_prefix = "q"
  ) |>
  dplyr::mutate(
    cov_95 = (observed >= .data[[paste0("q", COV_95_LOWER)]]) &
             (observed <= .data[[paste0("q", COV_95_UPPER)]])
  ) |>
  dplyr::select(model, location, origin_date, target_end_date,
                horizon, cov_95)


# ============================================================================
# BUILD SCORES_LONG_EXPANDED
# ============================================================================
cat("\n--- Building scores_long_expanded.csv ---\n")

medians <- forecasts_scored |>
  dplyr::filter(abs(quantile_level - 0.5) < 1e-9) |>
  dplyr::select(model, location, origin_date, target_end_date,
                horizon, point_forecast = predicted)

scores_long <- scores |>
  dplyr::transmute(
    model, location, origin_date, target_end_date, horizon,
    wis,
    wis_overprediction  = overprediction,
    wis_underprediction = underprediction,
    wis_dispersion      = dispersion,
    ae_median,
    cov_50 = interval_coverage_50
  ) |>
  dplyr::left_join(medians,
                   by = c("model", "location", "origin_date",
                          "target_end_date", "horizon")) |>
  dplyr::left_join(cov_95,
                   by = c("model", "location", "origin_date",
                          "target_end_date", "horizon")) |>
  dplyr::left_join(
    forecasts_scored |>
      dplyr::distinct(model, location, origin_date, target_end_date,
                      horizon, observed),
    by = c("model", "location", "origin_date", "target_end_date", "horizon")
  ) |>
  dplyr::select(model, location, origin_date, horizon, target_end_date,
                observed, point_forecast,
                wis, wis_overprediction, wis_underprediction, wis_dispersion,
                ae_median,
                cov_50, cov_95)

# Expected: 17 models * 11 locations * 4 horizons * 57 dates = 42,636 rows.
expected_rows <- N_MODELS * 11L * 4L * 57L
if (nrow(scores_long) != expected_rows) {
  warning(sprintf("scores_long has %d rows; expected %d", nrow(scores_long), expected_rows))
}
cat(sprintf("  %d rows (expected %d)\n", nrow(scores_long), expected_rows))

readr::write_csv(scores_long, file.path(PHASE4_DIR, "scores_long_expanded.csv"))


# ============================================================================
# PER-MODEL SUMMARY (sorted by WIS)
# ============================================================================
cat("\n--- Per-model summary ---\n")

per_model <- scores_long |>
  dplyr::group_by(model) |>
  dplyr::summarise(
    mae          = mean(ae_median, na.rm = TRUE),
    wis          = mean(wis, na.rm = TRUE),
    wis_over     = mean(wis_overprediction, na.rm = TRUE),
    wis_under    = mean(wis_underprediction, na.rm = TRUE),
    wis_disp     = mean(wis_dispersion, na.rm = TRUE),
    cov_50_rate  = mean(cov_50, na.rm = TRUE),
    cov_95_rate  = mean(cov_95, na.rm = TRUE),
    n_forecasts  = dplyr::n(),
    .groups      = "drop"
  ) |>
  dplyr::mutate(model_label = MODEL_LABELS[model]) |>
  dplyr::select(model_label, model, mae, wis, wis_over, wis_under,
                wis_disp, cov_50_rate, cov_95_rate, n_forecasts) |>
  dplyr::arrange(wis)

readr::write_csv(per_model, file.path(PHASE4_DIR, "per_model_summary_expanded.csv"))

cat("\n  Per-model WIS leaderboard (lower is better):\n")
print(per_model |>
        dplyr::select(model_label, wis, mae, cov_50_rate, cov_95_rate),
      n = Inf)


# ============================================================================
# CORRELATION MATRICES
#
# Same methodology as Phase 2: per-(location, horizon) Pearson correlation
# averaged across the 44 cells. Two matrices: point-error and WIS.
# ============================================================================
cat("\n--- Correlation matrices ---\n")

#' Compute average per-(location, horizon) correlation matrix for a metric.
compute_avg_correlation <- function(data, metric_col) {
  wide <- data |>
    dplyr::select(model, location, horizon, origin_date,
                  value = dplyr::all_of(metric_col)) |>
    tidyr::pivot_wider(names_from = model, values_from = value)

  groups <- wide |>
    dplyr::group_by(location, horizon) |>
    dplyr::group_split()

  cor_matrices <- purrr::map(groups, function(g) {
    mat <- g |> dplyr::select(dplyr::all_of(ALL_CANDIDATES)) |> as.matrix()
    stats::cor(mat, use = "complete.obs", method = "pearson")
  })

  Reduce(`+`, cor_matrices) / length(cor_matrices)
}

with_point_error <- scores_long |>
  dplyr::mutate(point_error = observed - point_forecast)

cor_point <- compute_avg_correlation(with_point_error, "point_error") |>
  round(3)
cor_wis <- compute_avg_correlation(scores_long, "wis") |>
  round(3)

# Replace row+column names with the shorter labels for readability.
new_names <- MODEL_LABELS[rownames(cor_point)]
rownames(cor_point) <- new_names; colnames(cor_point) <- new_names
rownames(cor_wis)   <- new_names; colnames(cor_wis)   <- new_names

# Save with the model labels as row+column identifiers.
write_cor_csv <- function(cor_mat, path) {
  out <- as.data.frame(cor_mat) |>
    tibble::rownames_to_column("model_label")
  readr::write_csv(out, path)
}

write_cor_csv(cor_point, file.path(PHASE4_DIR, "correlation_point_error_expanded.csv"))
write_cor_csv(cor_wis,   file.path(PHASE4_DIR, "correlation_wis_expanded.csv"))

# Print top extremes — these are what we care about for twin-pruning.
off_diag <- function(m) {
  m[upper.tri(m)] |>
    setNames(apply(which(upper.tri(m), arr.ind = TRUE), 1, function(idx) {
      paste(rownames(m)[idx[1]], colnames(m)[idx[2]], sep = " <-> ")
    }))
}

cat("\n  Top 10 highest point-error correlations (twin candidates):\n")
print(sort(off_diag(cor_point), decreasing = TRUE)[1:10])

cat("\n  Top 10 lowest point-error correlations (most diverse pairs):\n")
print(sort(off_diag(cor_point), decreasing = FALSE)[1:10])

# Empirical distribution summary — useful for picking the data-driven
# twin-pruning threshold in Phase 4B.
cat(sprintf(
  "\n  Off-diagonal correlation distribution:\n  N pairs = %d, min = %.3f, ",
  length(off_diag(cor_point)),
  min(off_diag(cor_point))
))
cat(sprintf("median = %.3f, mean = %.3f, max = %.3f\n",
            median(off_diag(cor_point)),
            mean(off_diag(cor_point)),
            max(off_diag(cor_point))))

cat("  Correlation quantiles for threshold selection:\n")
qs <- stats::quantile(off_diag(cor_point), c(0.50, 0.75, 0.90, 0.95, 0.99))
for (i in seq_along(qs)) {
  cat(sprintf("    %s percentile: %.3f\n", names(qs)[i], qs[i]))
}


# ============================================================================
# METADATA
# ============================================================================
metadata_path <- file.path(PHASE4_DIR, "scores_metadata_expanded.txt")
sink(metadata_path)
cat("Phase 4A: expanded candidate scoring run\n")
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("scoringutils: %s\n", packageVersion("scoringutils")))
cat(sprintf("\nCandidates (%d total):\n", N_MODELS))
cat(sprintf("  KReger (%d):\n", length(KREGER_CANDIDATES)))
cat(paste("    -", KREGER_CANDIDATES, collapse = "\n"))
cat(sprintf("\n  Upstream (%d):\n", length(UPSTREAM_CANDIDATES)))
cat(paste("    -", UPSTREAM_CANDIDATES, collapse = "\n"))
cat(sprintf("\n\nValidation cutoff: %s\n", VALIDATION_END_DATE))
cat(sprintf("Total forecast rows scored: %d\n", nrow(scores_long)))
cat("\nNote: KReger-theta_bc_bs and hist-avg deliberately excluded.\n")
cat("See analysis/phase2/decisions.md and Phase 4 design discussion.\n")
sink()


cat("\n========================================\n")
cat("  Phase 4A complete\n")
cat("========================================\n")
cat(sprintf("  Outputs in %s:\n", PHASE4_DIR))
cat("    - scores_long_expanded.csv          (42,636 rows)\n")
cat("    - per_model_summary_expanded.csv    (17 rows, sorted by WIS)\n")
cat("    - correlation_point_error_expanded.csv\n")
cat("    - correlation_wis_expanded.csv\n")
cat("    - scores_metadata_expanded.txt\n")
cat("========================================\n")
