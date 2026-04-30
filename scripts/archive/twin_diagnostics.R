###############################################################################
# FILE: scripts/twin_diagnostics.R
#
# PURPOSE:
#   Settle the three twin-pair selection decisions identified by Phase 2's
#   correlation analysis (FLUCAST_IMPLEMENTATION_PLAN.md Phase 2C). For
#   each twin pair, break out per-horizon and per-regime WIS to determine
#   whether one model dominates within specific subsets, or the pair is a
#   coin flip overall.
#
# TWIN PAIRS (from analysis/phase2/correlation_point_error.csv):
#   - ets_log <-> ets_bc           (point-error correlation 0.996)
#   - arima_bc_bs <-> arima_log    (0.984)
#   - nnetar_bc_bs <-> nnetar_log  (0.957)
#
# VERDICT RULE:
#   For each pair:
#     - Compute mean WIS in each (horizon, regime) cell, separately for
#       each model. There are 4 horizons * 4 regimes = 16 cells.
#     - In each cell, declare a winner if one model's WIS is >2% lower.
#       Otherwise call the cell a tie.
#     - If one model wins a strict majority of non-tie cells (more wins
#       than the other), it wins the pair.
#     - Otherwise (mostly ties, or wins evenly split): coin flip. Default
#       to the model with smaller mean wis_dispersion across the
#       validation set — sharper intervals when accuracy is tied.
#
# INPUTS:
#   - analysis/phase2/scores_long.csv  (from scripts/score_candidates.R)
#   - target-data/time-series.csv      (for regime classification)
#
# OUTPUTS:
#   - analysis/phase2/regimes.csv
#       One row per (origin_date, location, regime). 57 * 11 = 627 rows.
#       Saved as a side-effect of regime classification — useful for
#       Phase 3 ensemble work.
#   - analysis/phase2/twin_diagnostics.csv
#       Per-(pair, horizon, regime) WIS comparison. Long format with
#       columns: pair, horizon, regime, model_a, model_b, wis_a, wis_b,
#       diff, winner.
#   - Console output: per-pair verdict summary.
#
# HOW TO RUN:
#   setwd("/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026")
#   source("scripts/twin_diagnostics.R")
#
#   Runtime: ~1-2 min, dominated by regime calibration's slope cloud.
#
# REFERENCES:
#   - scripts/score_candidates.R (Phase 2 scoring; produces scores_long.csv)
#   - scripts/regime_helpers.R   (regime classifier)
###############################################################################


# ============================================================================
# DEPENDENCIES
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# Source the regime helpers (defines compute_regime_calibration,
# classify_regime, precompute_regimes).
source("scripts/regime_helpers.R")


# ============================================================================
# CONFIGURATION
# ============================================================================
HUB_PATH      <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"
TARGET_DATA   <- file.path(HUB_PATH, "target-data", "time-series.csv")
SCORES_LONG   <- file.path(HUB_PATH, "analysis", "phase2", "scores_long.csv")
ANALYSIS_DIR  <- file.path(HUB_PATH, "analysis", "phase2")

# Validation window (matches scripts/run_validation.R and score_candidates.R)
VALIDATION_START <- as.Date("2015-10-24")
VALIDATION_END   <- as.Date("2017-05-06")

# Twin pairs to evaluate. Order within each pair is (a, b); the verdict
# is reported as "<winner> over <loser>" using these names.
TWIN_PAIRS <- list(
  c("ets_log",      "ets_bc"),
  c("arima_bc_bs",  "arima_log"),
  c("nnetar_bc_bs", "nnetar_log")
)

# Cell-level WIS-margin threshold for declaring a winner (relative).
# 2% means: if model A's WIS is more than 2% lower than B's, A wins.
WIS_MARGIN <- 0.02


# ============================================================================
# LOAD DATA
# ============================================================================
cat("\n--- Loading scores and ground truth ---\n")

scores_long <- readr::read_csv(SCORES_LONG, show_col_types = FALSE) |>
  dplyr::mutate(
    origin_date     = as.Date(origin_date),
    target_end_date = as.Date(target_end_date)
  )

wili_data <- readr::read_csv(TARGET_DATA, show_col_types = FALSE) |>
  dplyr::transmute(
    location,
    target_end_date = as.Date(target_end_date),
    observation
  )

cat(sprintf("  scores_long: %d rows\n", nrow(scores_long)))
cat(sprintf("  wili_data:   %d rows\n", nrow(wili_data)))


# ============================================================================
# CALIBRATE REGIME CLASSIFIER
# ============================================================================
cat("\n--- Calibrating regime classifier ---\n")
cat("  (computing slope-z threshold and per-location calendar)\n")

calibration <- compute_regime_calibration(
  wili_data,
  validation_start_date = VALIDATION_START
)

cat(sprintf("  Slope-z threshold: %.3f (67th percentile of |z| in training data)\n",
            calibration$slope_z_threshold))
cat(sprintf("  Level thresholds calibrated for %d locations\n",
            nrow(calibration$level_thresholds)))
cat(sprintf("  Calendar fallback: %d (location, week_of_year) cells\n",
            nrow(calibration$location_calendar)))

# Sanity: regime distribution in the calendar fallback table.
calendar_dist <- calibration$location_calendar |>
  dplyr::count(regime) |>
  dplyr::arrange(dplyr::desc(n))
cat("\n  Calendar fallback regime distribution:\n")
print(calendar_dist)


# ============================================================================
# CLASSIFY VALIDATION ORIGIN DATES
# ============================================================================
cat("\n--- Classifying regimes for validation origin dates ---\n")

origin_dates <- sort(unique(scores_long$origin_date))
locations    <- sort(unique(scores_long$location))

cat(sprintf("  %d origin dates x %d locations = %d cells\n",
            length(origin_dates), length(locations),
            length(origin_dates) * length(locations)))

regime_lookup <- precompute_regimes(
  wili_data,
  origin_dates = origin_dates,
  locations    = locations,
  calibration  = calibration
)

# Save for Phase 3 reuse — no point recomputing.
readr::write_csv(regime_lookup, file.path(ANALYSIS_DIR, "regimes.csv"))
cat(sprintf("  Wrote: %s\n", file.path(ANALYSIS_DIR, "regimes.csv")))

# Distribution of regimes across the validation set.
cat("\n  Validation-period regime distribution:\n")
print(regime_lookup |> dplyr::count(regime) |> dplyr::arrange(dplyr::desc(n)))


# ============================================================================
# JOIN REGIMES INTO SCORES TABLE
# ============================================================================
scores_with_regime <- scores_long |>
  dplyr::inner_join(regime_lookup, by = c("origin_date", "location"))

# Sanity: row count must match — every scored forecast must have gotten
# a regime label (regimes are per-cell, scores have 4 horizons per cell,
# so a single regime label fans out to 4 score rows correctly via join).
if (nrow(scores_with_regime) != nrow(scores_long)) {
  stop(sprintf(
    "Regime join changed row count: %d -> %d. Investigate.",
    nrow(scores_long), nrow(scores_with_regime)
  ))
}


# ============================================================================
# PER-PAIR DIAGNOSTICS
# ============================================================================
cat("\n--- Computing per-pair (horizon, regime) WIS breakdown ---\n")

#' Compute WIS comparison table for one twin pair.
#' @param model_a Character. First model in the pair.
#' @param model_b Character. Second model in the pair.
#' @return Tibble with one row per (horizon, regime) cell.
pair_breakdown <- function(model_a, model_b) {
  scores_with_regime |>
    dplyr::filter(model %in% c(model_a, model_b)) |>
    dplyr::group_by(model, horizon, regime) |>
    dplyr::summarise(
      mean_wis  = mean(wis,            na.rm = TRUE),
      mean_disp = mean(wis_dispersion, na.rm = TRUE),
      n         = dplyr::n(),
      .groups   = "drop"
    ) |>
    tidyr::pivot_wider(
      id_cols     = c(horizon, regime),
      names_from  = model,
      values_from = c(mean_wis, mean_disp, n)
    ) |>
    # Build a stable wis_a / wis_b naming, independent of which model
    # is alphabetically first.
    dplyr::transmute(
      horizon, regime,
      wis_a   = .data[[paste0("mean_wis_",  model_a)]],
      wis_b   = .data[[paste0("mean_wis_",  model_b)]],
      disp_a  = .data[[paste0("mean_disp_", model_a)]],
      disp_b  = .data[[paste0("mean_disp_", model_b)]],
      n_a     = .data[[paste0("n_",         model_a)]],
      n_b     = .data[[paste0("n_",         model_b)]]
    ) |>
    dplyr::mutate(
      diff       = wis_a - wis_b,
      rel_diff   = diff / pmax(wis_a, wis_b),
      winner     = dplyr::case_when(
        rel_diff < -WIS_MARGIN  ~ "a_wins",
        rel_diff >  WIS_MARGIN  ~ "b_wins",
        TRUE                    ~ "tie"
      )
    )
}

# Build a single long table with all three pairs side-by-side.
all_pair_tables <- list()
for (pair in TWIN_PAIRS) {
  a <- pair[1]; b <- pair[2]
  tbl <- pair_breakdown(a, b) |>
    dplyr::mutate(
      pair    = sprintf("%s vs %s", a, b),
      model_a = a,
      model_b = b
    ) |>
    dplyr::relocate(pair, model_a, model_b)
  all_pair_tables[[sprintf("%s_vs_%s", a, b)]] <- tbl
}

twin_diagnostics <- dplyr::bind_rows(all_pair_tables)
readr::write_csv(twin_diagnostics, file.path(ANALYSIS_DIR, "twin_diagnostics.csv"))
cat(sprintf("  Wrote: %s\n", file.path(ANALYSIS_DIR, "twin_diagnostics.csv")))


# ============================================================================
# PER-PAIR VERDICTS
# ============================================================================

#' Apply the verdict rule for one pair.
#' @return Named character: c(verdict = ..., reason = ...).
verdict_for_pair <- function(model_a, model_b, pair_table) {
  pair_label <- sprintf("%s vs %s", model_a, model_b)
  # Use .data$pair to disambiguate from the calling scope's `pair` variable.
  cells <- pair_table |> dplyr::filter(.data$pair == pair_label)

  # Defensive: we should have exactly 12 cells (4 horizons * 3 regimes
  # observed in the validation period). If we don't, something upstream is
  # broken — fail loudly rather than computing a misleading verdict.
  if (nrow(cells) == 0L) {
    stop(sprintf("verdict_for_pair: no rows matched pair label '%s'. ",
                 pair_label),
         "Check that twin_diagnostics$pair contains the expected labels.")
  }
  if (nrow(cells) > 12L) {
    stop(sprintf(
      "verdict_for_pair: got %d rows for pair '%s' (expected up to 12). ",
      nrow(cells), pair_label
    ), "Filter is not restricting to one pair correctly.")
  }

  wins_a   <- sum(cells$winner == "a_wins")
  wins_b   <- sum(cells$winner == "b_wins")
  ties     <- sum(cells$winner == "tie")
  total    <- nrow(cells)

  # Get average dispersion across the validation set for the tiebreaker.
  disp_summary <- scores_with_regime |>
    dplyr::filter(model %in% c(model_a, model_b)) |>
    dplyr::group_by(model) |>
    dplyr::summarise(mean_disp = mean(wis_dispersion, na.rm = TRUE),
                     .groups = "drop")

  disp_a <- disp_summary$mean_disp[disp_summary$model == model_a]
  disp_b <- disp_summary$mean_disp[disp_summary$model == model_b]

  # Decision logic.
  if (wins_a > wins_b) {
    verdict <- model_a
    reason  <- sprintf("wins %d of %d cells (vs %d), %d ties",
                       wins_a, total, wins_b, ties)
  } else if (wins_b > wins_a) {
    verdict <- model_b
    reason  <- sprintf("wins %d of %d cells (vs %d), %d ties",
                       wins_b, total, wins_a, ties)
  } else {
    # Tie or coin flip. Default to smaller dispersion (sharper intervals).
    if (disp_a <= disp_b) {
      verdict <- model_a
      reason  <- sprintf(
        "coin flip (%d-%d-%d wins-losses-ties); %s has smaller dispersion (%.4f vs %.4f)",
        wins_a, wins_b, ties, model_a, disp_a, disp_b
      )
    } else {
      verdict <- model_b
      reason  <- sprintf(
        "coin flip (%d-%d-%d wins-losses-ties); %s has smaller dispersion (%.4f vs %.4f)",
        wins_a, wins_b, ties, model_b, disp_b, disp_a
      )
    }
  }

  c(verdict = verdict, reason = reason,
    wins_a = as.character(wins_a), wins_b = as.character(wins_b),
    ties = as.character(ties))
}


# ============================================================================
# PRINT VERDICTS
# ============================================================================
cat("\n========================================\n")
cat("  Twin-pair verdicts\n")
cat("========================================\n")

verdict_summary <- list()
for (pair in TWIN_PAIRS) {
  a <- pair[1]; b <- pair[2]
  v <- verdict_for_pair(a, b, twin_diagnostics)

  cat(sprintf("\n  %s vs %s:\n", a, b))
  cat(sprintf("    Verdict: KEEP %s\n", v[["verdict"]]))
  cat(sprintf("    Reason:  %s\n", v[["reason"]]))

  # Print the full cell breakdown for transparency.
  cells <- twin_diagnostics |>
    dplyr::filter(model_a == a, model_b == b) |>
    dplyr::transmute(
      horizon, regime,
      wis_a = round(wis_a, 4),
      wis_b = round(wis_b, 4),
      rel_diff_pct = sprintf("%+.1f%%", 100 * rel_diff),
      winner
    )
  cat("    Cell-by-cell breakdown:\n")
  print(cells, n = Inf)

  verdict_summary[[sprintf("%s_vs_%s", a, b)]] <- list(
    pair    = sprintf("%s vs %s", a, b),
    verdict = v[["verdict"]],
    reason  = v[["reason"]]
  )
}


# ============================================================================
# WRITE VERDICT SUMMARY
# ============================================================================
verdict_df <- tibble::tibble(
  pair    = vapply(verdict_summary, function(x) x$pair,    character(1)),
  verdict = vapply(verdict_summary, function(x) x$verdict, character(1)),
  reason  = vapply(verdict_summary, function(x) x$reason,  character(1))
)

readr::write_csv(verdict_df, file.path(ANALYSIS_DIR, "twin_verdicts.csv"))

cat("\n========================================\n")
cat("  Twin diagnostics complete\n")
cat("========================================\n")
cat(sprintf("  Outputs in %s:\n", ANALYSIS_DIR))
cat("    - regimes.csv               (627 rows; reusable for Phase 3)\n")
cat("    - twin_diagnostics.csv      (per-cell WIS breakdown)\n")
cat("    - twin_verdicts.csv         (per-pair decision)\n")
cat("========================================\n")
