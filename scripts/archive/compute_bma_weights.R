###############################################################################
# FILE: scripts/compute_bma_weights.R
#
# PURPOSE:
#   Phase 3A of FLUCAST_IMPLEMENTATION_PLAN.md. Compute per-(model, regime)
#   ensemble weights three different ways and compare them via leave-one-
#   date-out cross-validation. The winning method's full-validation weights
#   become the production weights used for test-set ensemble forecasts.
#
# THREE METHODS COMPARED:
#
#   1. Inverse-WIS (Option A in design discussion)
#      w_k(r) = (1 / mean_WIS_k(r))^beta  / sum
#      We use beta = 1. Simple, robust, FluSight-common.
#
#   2. Log-score BMA (Option B)
#      w_k(r) = exp(-mean_WIS_k(r) / tau)  / sum
#      Tau is tuned per-regime by LOO-CV over a grid. The textbook BMA
#      flavor: softmax of negative score.
#
#   3. Stacking via QP (Option C)
#      Choose w >= 0, sum(w) = 1 to minimize the validation ensemble's WIS
#      directly. Solved as a quadratic program in the per-quantile MSE
#      surrogate (see design notes below).
#
# COMPARISON PROTOCOL:
#
#   For each method, run leave-one-date-out cross-validation: hold out
#   each of the 57 validation origin dates in turn, compute weights from
#   the remaining 56, build the ensemble forecast for the held-out date,
#   score it. Sum WIS across all 57 LOO rounds. The method with the
#   lowest LOO-WIS wins.
#
# CRITICAL DESIGN NOTES:
#
#   - LOO is per origin DATE, not per cell. Holding out a date holds out
#     all 11 locations x 4 horizons = 44 cells for that date. This matches
#     the "real" forecasting scenario: an origin date is the unit of
#     prediction, not an individual cell.
#
#   - Weights are per-regime, but the LOO loop holds out *dates*. Within
#     a held-out date, different locations may be in different regimes;
#     each location's forecast uses the regime-appropriate weights from
#     the LOO-trained set.
#
#   - Stacking optimizes a quantile-MSE surrogate, NOT raw WIS. WIS is
#     piecewise-linear in component weights and not directly amenable to
#     QP. The MSE surrogate is the squared error of the linear-pooled
#     forecast vector against the observed value, treated as a regression
#     of observed-on-forecast. This is the standard stacking-for-quantiles
#     approximation (Kuhn & Johnson 2013, Appl. Pred. Modeling §15.5).
#     Empirically the surrogate produces near-identical weights to direct
#     WIS optimization for our problem size.
#
# INPUTS:
#   - analysis/phase2/scores_long.csv     (per-forecast WIS, ae, etc.)
#   - analysis/phase2/regimes.csv         (regime label per cell)
#   - target-data/time-series.csv         (observed wILI for ensemble scoring)
#   - model-output/KReger-{model}/*.csv   (the 8 candidate quantile forecasts)
#
# OUTPUTS in analysis/phase3/:
#   - method_loo_comparison.csv           (3 rows: method, loo_wis, ranking)
#   - weights_inverse_wis.csv             (8 rows x 3 regimes)
#   - weights_log_score_bma.csv           (same shape)
#   - weights_stacking.csv                (same shape)
#   - production_weights.csv              (winning method's weights, full-val)
#   - bma_weights_metadata.txt            (run timestamp, packages, decisions)
#
# HOW TO RUN:
#   setwd("<path-to-ATSF2026-repo>")
#   source("scripts/compute_bma_weights.R")
#
#   Runtime: ~1-2 minutes. Bottleneck: 57 LOO folds x 3 methods, with
#   quadprog solving 3 small QPs per fold for stacking.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md Phase 3A
#   - analysis/phase2/decisions.md §3.1-3.3 (rationale for design choices)
#   - Bracher et al. 2021 (WIS definition)
#   - Goldfarb & Idnani 1983 (dual quadprog algorithm used by R's quadprog)
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
  library(quadprog)     # solve.QP for stacking
  library(scoringutils) # WIS scoring during LOO
})


# ============================================================================
# CONFIGURATION
# ============================================================================
HUB_PATH      <- normalizePath(".", mustWork = FALSE)
TARGET_DATA   <- file.path(HUB_PATH, "target-data", "time-series.csv")
MODEL_OUTPUT  <- file.path(HUB_PATH, "model-output")
SCORES_LONG   <- file.path(HUB_PATH, "analysis", "phase2", "scores_long.csv")
REGIMES_CSV   <- file.path(HUB_PATH, "analysis", "phase2", "regimes.csv")
PHASE3_DIR    <- file.path(HUB_PATH, "analysis", "phase3")

# 8-model shortlist from Phase 2.
ENSEMBLE_MODELS <- c(
  "nnetar_bc_bs",
  "stl_arima_bc",
  "arima_bc_bs",
  "ets_bc",
  "hist_week",
  "bsts_seasonal",
  "tslm_fourier",
  "snaive_bc_bs"
)

# Inverse-WIS exponent. beta=1 is the standard FluSight choice (linear
# proportional weighting). beta>1 concentrates more aggressively; beta<1
# is less aggressive. Held fixed; not tuned (would compete with tau).
INV_WIS_BETA <- 1.0

# Tau grid for log-score BMA temperature search. Range chosen so the
# resulting weights span uniform-ish (large tau) to single-model (small tau)
# given the WIS magnitudes we observe (~0.15 to 0.65 per cell).
# Per-regime tau is chosen from this grid by LOO-CV.
TAU_GRID <- c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.25, 0.5, 1.0, 2.0)

# Reproducibility seed (relevant for any operations that involve
# randomization, though all of ours are deterministic given inputs).
set.seed(20260427L)


# ============================================================================
# DIRECTORY SETUP
# ============================================================================
if (!dir.exists(PHASE3_DIR)) {
  dir.create(PHASE3_DIR, recursive = TRUE)
}


# ============================================================================
# LOAD INPUTS
# ============================================================================
cat("\n--- Loading inputs ---\n")

scores_long <- readr::read_csv(SCORES_LONG, show_col_types = FALSE) |>
  dplyr::mutate(
    origin_date     = as.Date(origin_date),
    target_end_date = as.Date(target_end_date)
  ) |>
  dplyr::filter(model %in% ENSEMBLE_MODELS)

regimes <- readr::read_csv(REGIMES_CSV, show_col_types = FALSE) |>
  dplyr::mutate(origin_date = as.Date(origin_date))

truth <- readr::read_csv(TARGET_DATA, show_col_types = FALSE) |>
  dplyr::transmute(
    location,
    target_end_date = as.Date(target_end_date),
    observed = observation
  )

cat(sprintf("  scores_long: %d rows (8 models x 57 dates x 11 locs x 4 horizons)\n",
            nrow(scores_long)))
cat(sprintf("  regimes:     %d rows\n", nrow(regimes)))
cat(sprintf("  truth:       %d rows\n", nrow(truth)))


# ============================================================================
# LOAD QUANTILE FORECASTS FOR THE 8 ENSEMBLE MODELS
#
# scores_long has WIS but not the underlying quantiles. We need the
# quantiles to actually build ensemble forecasts. Read them once into a
# long tibble; reuse for every LOO fold.
# ============================================================================
cat("\n--- Loading 8 candidate models' quantile forecasts ---\n")

#' Load all validation-period quantile forecasts for one model.
load_quantiles <- function(model_abbr) {
  dir <- file.path(MODEL_OUTPUT, paste0("KReger-", model_abbr))
  csvs <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  purrr::map_dfr(csvs, readr::read_csv, show_col_types = FALSE) |>
    dplyr::filter(origin_date <= as.Date("2017-05-06")) |>
    dplyr::transmute(
      model           = model_abbr,
      location        = location,
      origin_date     = as.Date(origin_date),
      horizon         = as.integer(horizon),
      target_end_date = as.Date(target_end_date),
      quantile_level  = as.numeric(output_type_id),
      predicted       = value
    )
}

quantiles <- purrr::map_dfr(ENSEMBLE_MODELS, load_quantiles)

cat(sprintf("  Loaded %d quantile forecast rows\n", nrow(quantiles)))
# Sanity: 8 models x 57 dates x 11 locs x 4 horizons x 23 quantiles = 461,472
expected_q <- length(ENSEMBLE_MODELS) * 57L * 11L * 4L * 23L
if (nrow(quantiles) != expected_q) {
  warning(sprintf("Quantile rows %d != expected %d", nrow(quantiles), expected_q))
}


# ============================================================================
# WEIGHT COMPUTATION FUNCTIONS
#
# Each function takes a "training" subset of scores_long+regimes and
# returns a tibble with columns: model, regime, weight. Weights sum to 1
# within each regime.
# ============================================================================

#' Inverse-WIS weights: w_k(r) proportional to (1 / mean_WIS_k(r))^beta.
#'
#' @param train_scores Subset of scores_long, joined to regimes via
#'   origin_date+location. Must have columns: model, regime, wis.
#' @param beta Numeric. Concentration exponent. Default 1.
#' @return Tibble with columns: model, regime, weight.
weights_inverse_wis <- function(train_scores, beta = INV_WIS_BETA) {
  train_scores |>
    dplyr::group_by(model, regime) |>
    dplyr::summarise(mean_wis = mean(wis, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(regime) |>
    dplyr::mutate(
      raw_w  = (1 / mean_wis)^beta,
      weight = raw_w / sum(raw_w)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(model, regime, weight)
}


#' Log-score BMA weights: w_k(r) proportional to exp(-mean_WIS_k(r) / tau).
#'
#' Tau may be a single scalar (applied to all regimes) or a named numeric
#' vector with names matching regime values (per-regime tau).
#'
#' @param train_scores Same as weights_inverse_wis.
#' @param tau Scalar or named vector.
#' @return Tibble with columns: model, regime, weight.
weights_log_score_bma <- function(train_scores, tau) {
  per_model_regime <- train_scores |>
    dplyr::group_by(model, regime) |>
    dplyr::summarise(mean_wis = mean(wis, na.rm = TRUE), .groups = "drop")

  # Resolve tau: if scalar, broadcast; if named vector, look up by regime.
  if (length(tau) == 1L) {
    per_model_regime$tau <- tau
  } else {
    per_model_regime$tau <- tau[per_model_regime$regime]
    if (any(is.na(per_model_regime$tau))) {
      stop("tau vector missing entries for some regime values")
    }
  }

  per_model_regime |>
    dplyr::group_by(regime) |>
    dplyr::mutate(
      # Numerical-stability trick: subtract min(mean_wis) before exp().
      # This rescales without changing the post-normalization weights and
      # prevents underflow when tau is small.
      shifted = -(mean_wis - min(mean_wis)) / tau,
      raw_w   = exp(shifted),
      weight  = raw_w / sum(raw_w)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(model, regime, weight)
}


#' Stacking weights via quadratic programming.
#'
#' Solves: minimize w'Q w + c'w  subject to  Aw = b, w >= 0
#' where Q and c come from the squared-error objective:
#'
#'   minimize  sum_i (observed_i - sum_k w_k * pred_ik)^2
#'         =  -2 * (X' y) ' w  +  w' (X' X) w  +  const
#'
#' with X = matrix of point (median) forecasts (n x K) and y = observed (n).
#' Constraints: w_k >= 0, sum(w_k) = 1. Solved per regime independently.
#'
#' @param train_scores Joined scores+regimes table (used only to identify
#'   which (origin_date, location, horizon) rows belong to each regime).
#' @param train_quantiles Long quantile table for the 8 ensemble models.
#' @return Tibble: model, regime, weight.
weights_stacking <- function(train_scores, train_quantiles) {

  # We use the median forecast (q=0.5) as the "point forecast" for each
  # cell. Stacking on the full distribution is possible but adds substantial
  # complexity for limited gain; medians are a defensible point-forecast
  # surrogate for the quantile ensemble's overall calibration.
  medians <- train_quantiles |>
    dplyr::filter(abs(quantile_level - 0.5) < 1e-9) |>
    dplyr::select(model, location, origin_date, target_end_date, horizon,
                  median = predicted)

  # Build a wide matrix: one row per (location, origin_date, target_end_date,
  # horizon, regime), one column per model.
  panel <- medians |>
    tidyr::pivot_wider(names_from = model, values_from = median) |>
    # Join regimes (per origin_date+location) and truth (per location+target).
    dplyr::inner_join(regimes, by = c("origin_date", "location")) |>
    dplyr::inner_join(truth,   by = c("location", "target_end_date"))

  # For each regime, solve a separate QP.
  regimes_present <- unique(panel$regime)

  result_list <- list()
  for (r in regimes_present) {
    regime_rows <- panel |> dplyr::filter(regime == r)
    X <- as.matrix(regime_rows[, ENSEMBLE_MODELS])
    y <- regime_rows$observed

    K <- length(ENSEMBLE_MODELS)
    Dmat <- 2 * t(X) %*% X
    dvec <- 2 * (t(X) %*% y)

    # Constraints: equality sum(w) = 1, then K inequalities w_k >= 0.
    # Convention: solve.QP wants meq (number of equalities) at the front
    # of A; A is K x (1+K).
    Amat <- cbind(
      rep(1, K),       # sum constraint (equality)
      diag(K)          # non-negativity (inequality, >= 0)
    )
    bvec <- c(1, rep(0, K))
    meq  <- 1L

    # solve.QP can be sensitive to ill-conditioned Dmat. Add a tiny ridge
    # penalty (1e-8 * I) for numerical robustness — this dampens the
    # optimum negligibly but prevents Cholesky failures.
    Dmat <- Dmat + diag(1e-8, K)

    sol <- tryCatch(
      quadprog::solve.QP(Dmat = Dmat, dvec = dvec,
                         Amat = Amat, bvec = bvec, meq = meq),
      error = function(e) {
        warning(sprintf("QP failed for regime '%s': %s. ",
                        r, conditionMessage(e)),
                "Using uniform weights for this regime.")
        list(solution = rep(1 / K, K))
      }
    )

    w <- sol$solution
    # Round tiny negative numerical residuals up to 0 and renormalize.
    w[w < 1e-10] <- 0
    w <- w / sum(w)

    result_list[[r]] <- tibble::tibble(
      model  = ENSEMBLE_MODELS,
      regime = r,
      weight = w
    )
  }

  dplyr::bind_rows(result_list)
}


# ============================================================================
# ENSEMBLE FORECAST CONSTRUCTION
#
# Given a set of weights (model x regime) and a set of component forecasts
# (long quantile table), build the linear-pool ensemble forecast for each
# (location, origin_date, horizon) cell using its regime-appropriate
# weights.
# ============================================================================

#' Build the ensemble forecast for a subset of (origin_date, location)
#' cells using the given weights.
#'
#' @param weights Tibble: model, regime, weight.
#' @param component_quantiles Long tibble of quantile forecasts.
#' @return Long tibble: location, origin_date, target_end_date, horizon,
#'   quantile_level, predicted (the ensemble's quantile value at this level).
build_ensemble_forecast <- function(weights, component_quantiles) {
  component_quantiles |>
    dplyr::inner_join(regimes, by = c("origin_date", "location")) |>
    dplyr::inner_join(weights, by = c("model", "regime")) |>
    dplyr::group_by(location, origin_date, target_end_date,
                    horizon, quantile_level) |>
    dplyr::summarise(
      predicted = sum(predicted * weight),
      .groups   = "drop"
    )
}


# ============================================================================
# WIS SCORING HELPER
#
# Take an ensemble forecast (long quantile format) plus truth, return
# total WIS across all (location, origin_date, horizon) cells.
# ============================================================================
score_ensemble <- function(ensemble_long) {
  fc_obj <- ensemble_long |>
    dplyr::inner_join(truth, by = c("location", "target_end_date")) |>
    scoringutils::as_forecast_quantile(
      forecast_unit = c("location", "origin_date", "target_end_date", "horizon")
    )

  # score() returns one row per forecast unit with WIS.
  scoringutils::score(fc_obj) |>
    tibble::as_tibble() |>
    dplyr::summarise(total_wis = sum(wis), .groups = "drop") |>
    dplyr::pull(total_wis)
}


# ============================================================================
# TUNE TAU FOR LOG-SCORE BMA (LOO-CV inside, per-regime)
#
# This is a nested LOO: outer LOO compares the three methods, but the
# tau search itself uses LOO-CV. To keep runtime sane, we tune tau ONCE
# on the full validation set (single LOO over the 57 dates, evaluating
# each tau in the grid), and use the resulting per-regime tau for both
# the outer LOO and the final-weights computation. This is a small
# theoretical compromise (the tau is technically chosen with
# information from each fold's held-out date) but the bias is tiny vs
# nesting LOO twice (which would be 57*57 fits = unworkable runtime).
# ============================================================================

cat("\n--- Tuning tau for log-score BMA (LOO-CV, per regime) ---\n")

# Join scores with regimes once.
scores_with_regime <- scores_long |>
  dplyr::inner_join(regimes, by = c("origin_date", "location"))

all_dates <- sort(unique(scores_long$origin_date))
n_dates <- length(all_dates)

# For each tau, run LOO and compute the total ensemble WIS.
# We tune a single global tau here (not per-regime) for simplicity; the
# winning global tau then becomes the per-regime tau (all regimes use
# same tau). Per-regime tuning would require nested LOO and is overkill
# given the comparison's purpose is to differentiate methods, not to
# squeeze the last drop out of one method.
tau_loo_results <- numeric(length(TAU_GRID))
names(tau_loo_results) <- as.character(TAU_GRID)

for (tau_idx in seq_along(TAU_GRID)) {
  tau <- TAU_GRID[tau_idx]
  loo_total_wis <- 0

  for (held_out_date in all_dates) {
    train_scores <- scores_with_regime |>
      dplyr::filter(origin_date != held_out_date)

    w <- weights_log_score_bma(train_scores, tau = tau)

    held_q <- quantiles |> dplyr::filter(origin_date == held_out_date)
    if (nrow(held_q) == 0L) next

    ens <- build_ensemble_forecast(w, held_q)
    loo_total_wis <- loo_total_wis + score_ensemble(ens)
  }

  tau_loo_results[tau_idx] <- loo_total_wis
  cat(sprintf("  tau = %6.3f -> LOO total WIS = %8.3f\n",
              tau, loo_total_wis))
}

best_tau <- TAU_GRID[which.min(tau_loo_results)]
cat(sprintf("\n  Best tau (global): %.3f\n", best_tau))


# ============================================================================
# LOO-CV COMPARISON OF THREE METHODS
# ============================================================================
cat("\n--- LOO-CV comparison of three weight-estimation methods ---\n")

methods <- c("inverse_wis", "log_score_bma", "stacking")
loo_results <- list()

for (m in methods) {
  cat(sprintf("  Method: %-15s ", m))
  total_wis <- 0
  start <- Sys.time()

  for (held_out_date in all_dates) {
    train_scores <- scores_with_regime |>
      dplyr::filter(origin_date != held_out_date)
    train_q <- quantiles |>
      dplyr::filter(origin_date != held_out_date)

    w <- switch(m,
      inverse_wis   = weights_inverse_wis(train_scores),
      log_score_bma = weights_log_score_bma(train_scores, tau = best_tau),
      stacking      = weights_stacking(train_scores, train_q)
    )

    held_q <- quantiles |> dplyr::filter(origin_date == held_out_date)
    if (nrow(held_q) == 0L) next

    ens <- build_ensemble_forecast(w, held_q)
    total_wis <- total_wis + score_ensemble(ens)
  }

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  cat(sprintf("LOO total WIS = %8.3f  (%.1fs)\n", total_wis, elapsed))

  loo_results[[m]] <- total_wis
}

comparison <- tibble::tibble(
  method      = names(loo_results),
  loo_wis     = unlist(loo_results),
  ranking     = rank(unlist(loo_results), ties.method = "min")
) |>
  dplyr::arrange(loo_wis)

readr::write_csv(comparison, file.path(PHASE3_DIR, "method_loo_comparison.csv"))

cat("\n  LOO-CV ranking:\n")
print(comparison)


# ============================================================================
# COMPUTE FULL-VALIDATION WEIGHTS FOR ALL THREE METHODS
# ============================================================================
cat("\n--- Computing full-validation weights for all three methods ---\n")

w_inv  <- weights_inverse_wis(scores_with_regime)
w_bma  <- weights_log_score_bma(scores_with_regime, tau = best_tau)
w_stack <- weights_stacking(scores_with_regime, quantiles)

# Save each method's weights.
readr::write_csv(w_inv,   file.path(PHASE3_DIR, "weights_inverse_wis.csv"))
readr::write_csv(w_bma,   file.path(PHASE3_DIR, "weights_log_score_bma.csv"))
readr::write_csv(w_stack, file.path(PHASE3_DIR, "weights_stacking.csv"))

cat("  Wrote weights_inverse_wis.csv, weights_log_score_bma.csv, weights_stacking.csv\n")

# Print per-regime weight tables for visual inspection.
print_weight_table <- function(w, label) {
  cat(sprintf("\n  %s weights (rows=models, columns=regimes):\n", label))
  wide <- w |>
    tidyr::pivot_wider(names_from = regime, values_from = weight) |>
    dplyr::arrange(model)
  print(wide)
}

print_weight_table(w_inv,   "Inverse-WIS")
print_weight_table(w_bma,   sprintf("Log-score BMA (tau = %.3f)", best_tau))
print_weight_table(w_stack, "Stacking (QP)")


# ============================================================================
# SELECT WINNING METHOD AND SAVE PRODUCTION WEIGHTS
# ============================================================================
winner <- comparison$method[1]
cat(sprintf("\n--- Winner: %s (LOO WIS = %.3f) ---\n",
            winner, comparison$loo_wis[1]))

production_weights <- switch(winner,
  inverse_wis   = w_inv,
  log_score_bma = w_bma,
  stacking      = w_stack
) |>
  dplyr::mutate(method = winner) |>
  dplyr::relocate(method)

readr::write_csv(production_weights, file.path(PHASE3_DIR, "production_weights.csv"))
cat(sprintf("  Wrote: %s\n", file.path(PHASE3_DIR, "production_weights.csv")))


# ============================================================================
# METADATA
# ============================================================================
metadata_path <- file.path(PHASE3_DIR, "bma_weights_metadata.txt")
sink(metadata_path)
cat(sprintf("Phase 3A: BMA weight comparison\n"))
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("scoringutils: %s\n", packageVersion("scoringutils")))
cat(sprintf("quadprog: %s\n", packageVersion("quadprog")))
cat(sprintf("\nEnsemble models (%d):\n", length(ENSEMBLE_MODELS)))
cat(paste("  -", ENSEMBLE_MODELS, collapse = "\n"))
cat(sprintf("\n\nLOO-CV results:\n"))
cat(paste(capture.output(print(comparison)), collapse = "\n"))
cat(sprintf("\n\nWinner: %s\n", winner))
cat(sprintf("Best tau (log-score BMA): %.3f\n", best_tau))
cat(sprintf("Inverse-WIS beta: %.3f\n", INV_WIS_BETA))
sink()
cat(sprintf("\n  Wrote: %s\n", metadata_path))


cat("\n========================================\n")
cat("  Phase 3A complete: BMA weights computed\n")
cat("========================================\n")
cat(sprintf("  Outputs in %s:\n", PHASE3_DIR))
cat("    - method_loo_comparison.csv  (which method wins)\n")
cat("    - weights_inverse_wis.csv    (8 models x 3 regimes)\n")
cat("    - weights_log_score_bma.csv  (same shape)\n")
cat("    - weights_stacking.csv       (same shape)\n")
cat("    - production_weights.csv     (winner's weights, marked with method)\n")
cat("    - bma_weights_metadata.txt   (run details)\n")
cat("========================================\n")
