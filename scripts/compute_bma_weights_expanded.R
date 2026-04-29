###############################################################################
# FILE: scripts/compute_bma_weights_expanded.R
#
# PURPOSE:
#   Phase 4C of FLUCAST_IMPLEMENTATION_PLAN.md. Compute per-(model, regime)
#   ensemble weights for the EXPANDED 10-model shortlist from Phase 4B,
#   following the same three-way methodology used for the primary 8-model
#   ensemble (Phase 3A).
#
#   Three weight-estimation methods compared:
#     1. Inverse-WIS:    w_k(r) propto (1 / mean_WIS_k(r))^beta
#     2. Log-score BMA:  w_k(r) propto exp(-mean_WIS_k(r) / tau)
#     3. Stacking:       solve QP minimizing ensemble median MSE per regime
#
#   Comparison protocol: leave-one-date-out cross-validation. The
#   winning method's full-validation weights become the production
#   weights used to build the expanded ensemble in Phase 4D.
#
# DIFFERS FROM scripts/compute_bma_weights.R:
#   - ENSEMBLE_MODELS list comes from analysis/phase4/shortlist_expanded.csv
#     (10 models including delphi-epicast, slm_naive, atsf_meanNB)
#   - Outputs go to analysis/phase4/, prefixed with _expanded
#   - Otherwise identical methodology — important for Phase 4D's apples-
#     to-apples comparison between primary and expanded BMA ensembles
#
# CRITICAL: Same data hygiene as Phase 3A.
#   Weights estimated from validation-period (57 dates) only. Test seasons
#   remain unseen. Production weights are frozen — applied unchanged
#   when Phase 4D generates 134-date ensemble forecasts.
#
# INPUTS:
#   - analysis/phase4/shortlist_expanded.csv   (10 models, post-twin-pruning)
#   - analysis/phase4/scores_long_expanded.csv (per-forecast WIS for those 10)
#   - analysis/phase2/regimes.csv              (regime label per cell)
#   - target-data/time-series.csv              (truth for ensemble scoring)
#   - model-output/<model_dir>/*.csv           (10 candidates' quantiles)
#
# OUTPUTS in analysis/phase4/:
#   - method_loo_comparison_expanded.csv       (3 rows: method, loo_wis, ranking)
#   - weights_inverse_wis_expanded.csv         (10 rows x N_regimes)
#   - weights_log_score_bma_expanded.csv       (same shape)
#   - weights_stacking_expanded.csv            (same shape)
#   - production_weights_expanded.csv          (winning method's weights)
#   - bma_weights_expanded_metadata.txt        (run details)
#
# HOW TO RUN:
#   setwd("/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026")
#   source("scripts/compute_bma_weights_expanded.R")
#
#   Runtime: ~1-2 minutes (similar to Phase 3A).
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md Phase 4C
#   - scripts/compute_bma_weights.R (template)
#   - analysis/phase2/decisions.md §3.1-3.3 (rationale for design choices)
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
  library(quadprog)
  library(scoringutils)
})


# ============================================================================
# CONFIGURATION
# ============================================================================
HUB_PATH       <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"
TARGET_DATA    <- file.path(HUB_PATH, "target-data", "time-series.csv")
MODEL_OUTPUT   <- file.path(HUB_PATH, "model-output")

PHASE4_DIR     <- file.path(HUB_PATH, "analysis", "phase4")
PHASE2_DIR     <- file.path(HUB_PATH, "analysis", "phase2")

SHORTLIST_CSV  <- file.path(PHASE4_DIR, "shortlist_expanded.csv")
SCORES_LONG    <- file.path(PHASE4_DIR, "scores_long_expanded.csv")
REGIMES_CSV    <- file.path(PHASE2_DIR, "regimes.csv")

# Same hyperparameters as Phase 3A — keep methodology identical.
INV_WIS_BETA <- 1.0
TAU_GRID     <- c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.25, 0.5, 1.0, 2.0)
set.seed(20260428L)


# ============================================================================
# LOAD SHORTLIST AND DERIVE ENSEMBLE_MODELS
# ============================================================================
cat("\n--- Loading expanded shortlist ---\n")

shortlist <- readr::read_csv(SHORTLIST_CSV, show_col_types = FALSE)

# `model` column has the directory names; `model_label` has the short
# version we'll use in printed tables.
ENSEMBLE_MODELS <- shortlist$model
MODEL_LABELS    <- setNames(shortlist$model_label, shortlist$model)

cat(sprintf("  %d models in expanded ensemble:\n", length(ENSEMBLE_MODELS)))
for (m in ENSEMBLE_MODELS) {
  cat(sprintf("    %s  (label: %s)\n", m, MODEL_LABELS[m]))
}


# ============================================================================
# LOAD INPUTS
# ============================================================================
cat("\n--- Loading scores and regimes ---\n")

scores_long <- readr::read_csv(SCORES_LONG, show_col_types = FALSE) |>
  dplyr::mutate(
    origin_date     = as.Date(origin_date),
    target_end_date = as.Date(target_end_date)
  ) |>
  # Filter to just the shortlist's models (scores_long_expanded.csv has
  # all 20; we only need the 10 survivors for weight estimation).
  dplyr::filter(model %in% ENSEMBLE_MODELS)

regimes <- readr::read_csv(REGIMES_CSV, show_col_types = FALSE) |>
  dplyr::mutate(origin_date = as.Date(origin_date))

truth <- readr::read_csv(TARGET_DATA, show_col_types = FALSE) |>
  dplyr::transmute(
    location, target_end_date = as.Date(target_end_date),
    observed = observation
  )

cat(sprintf("  scores_long: %d rows (%d models x 57 dates x 11 locs x 4 horizons)\n",
            nrow(scores_long), length(ENSEMBLE_MODELS)))
cat(sprintf("  regimes:     %d rows\n", nrow(regimes)))
cat(sprintf("  truth:       %d rows\n", nrow(truth)))


# ============================================================================
# LOAD QUANTILE FORECASTS FOR THE 10 ENSEMBLE MODELS
#
# scores_long has WIS but not the underlying quantiles — we need the
# quantile values to actually build ensemble forecasts.
# ============================================================================
cat("\n--- Loading quantile forecasts for the 10 ensemble models ---\n")

#' Load all validation-period quantile forecasts for one model.
load_quantiles <- function(model_dir) {
  dir <- file.path(MODEL_OUTPUT, model_dir)
  csvs <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  loaded <- purrr::map_dfr(
    csvs,
    readr::read_csv,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_guess())
  )
  loaded |>
    dplyr::filter(as.Date(origin_date) <= as.Date("2017-05-06")) |>
    dplyr::transmute(
      model           = model_dir,
      location        = as.character(location),
      origin_date     = as.Date(origin_date),
      horizon         = as.integer(horizon),
      target_end_date = as.Date(target_end_date),
      quantile_level  = as.numeric(output_type_id),
      predicted       = as.numeric(value)
    )
}

quantiles <- purrr::map_dfr(ENSEMBLE_MODELS, load_quantiles)

cat(sprintf("  Loaded %d quantile rows\n", nrow(quantiles)))
expected_q <- length(ENSEMBLE_MODELS) * 57L * 11L * 4L * 23L
if (nrow(quantiles) != expected_q) {
  warning(sprintf("Quantile rows %d != expected %d", nrow(quantiles), expected_q))
}


# ============================================================================
# WEIGHT COMPUTATION FUNCTIONS
#
# Same three methods as Phase 3A. Code is essentially copied — methodology
# must be identical for Phase 4D's primary-vs-expanded comparison to be
# apples-to-apples. Comments here are abbreviated; full discussion is in
# scripts/compute_bma_weights.R.
# ============================================================================

#' Inverse-WIS weights: w_k(r) propto (1 / mean_WIS_k(r))^beta.
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

#' Log-score BMA weights: w_k(r) propto exp(-mean_WIS_k(r) / tau).
weights_log_score_bma <- function(train_scores, tau) {
  per_model_regime <- train_scores |>
    dplyr::group_by(model, regime) |>
    dplyr::summarise(mean_wis = mean(wis, na.rm = TRUE), .groups = "drop")

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
      shifted = -(mean_wis - min(mean_wis)) / tau,
      raw_w   = exp(shifted),
      weight  = raw_w / sum(raw_w)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(model, regime, weight)
}

#' Stacking weights via QP on median forecasts (per regime).
weights_stacking <- function(train_scores, train_quantiles) {

  medians <- train_quantiles |>
    dplyr::filter(abs(quantile_level - 0.5) < 1e-9) |>
    dplyr::select(model, location, origin_date, target_end_date, horizon,
                  median = predicted)

  panel <- medians |>
    tidyr::pivot_wider(names_from = model, values_from = median) |>
    dplyr::inner_join(regimes, by = c("origin_date", "location")) |>
    dplyr::inner_join(truth,   by = c("location", "target_end_date"))

  regimes_present <- unique(panel$regime)
  result_list <- list()

  for (r in regimes_present) {
    regime_rows <- panel |> dplyr::filter(regime == r)
    X <- as.matrix(regime_rows[, ENSEMBLE_MODELS])
    y <- regime_rows$observed

    K <- length(ENSEMBLE_MODELS)
    Dmat <- 2 * t(X) %*% X
    dvec <- 2 * (t(X) %*% y)

    Amat <- cbind(rep(1, K), diag(K))
    bvec <- c(1, rep(0, K))
    meq  <- 1L

    Dmat <- Dmat + diag(1e-8, K)

    sol <- tryCatch(
      quadprog::solve.QP(Dmat = Dmat, dvec = dvec,
                         Amat = Amat, bvec = bvec, meq = meq),
      error = function(e) {
        warning(sprintf("QP failed for regime '%s': %s", r, conditionMessage(e)))
        list(solution = rep(1 / K, K))
      }
    )

    w <- sol$solution
    w[w < 1e-10] <- 0
    w <- w / sum(w)

    result_list[[r]] <- tibble::tibble(
      model = ENSEMBLE_MODELS, regime = r, weight = w
    )
  }

  dplyr::bind_rows(result_list)
}


# ============================================================================
# ENSEMBLE FORECAST CONSTRUCTION
# ============================================================================

#' Build the ensemble forecast for a subset of (origin_date, location)
#' cells using the given weights.
build_ensemble_forecast <- function(weights, component_quantiles) {
  component_quantiles |>
    dplyr::inner_join(regimes, by = c("origin_date", "location")) |>
    dplyr::inner_join(weights, by = c("model", "regime")) |>
    dplyr::group_by(location, origin_date, target_end_date,
                    horizon, quantile_level) |>
    dplyr::summarise(predicted = sum(predicted * weight), .groups = "drop")
}


# ============================================================================
# WIS SCORING HELPER
# ============================================================================
score_ensemble <- function(ensemble_long) {
  fc_obj <- ensemble_long |>
    dplyr::inner_join(truth, by = c("location", "target_end_date")) |>
    scoringutils::as_forecast_quantile(
      forecast_unit = c("location", "origin_date", "target_end_date", "horizon")
    )

  scoringutils::score(fc_obj) |>
    tibble::as_tibble() |>
    dplyr::summarise(total_wis = sum(wis), .groups = "drop") |>
    dplyr::pull(total_wis)
}


# ============================================================================
# TAU TUNING (LOO over the full validation set)
# ============================================================================
cat("\n--- Tuning tau for log-score BMA (LOO-CV) ---\n")

scores_with_regime <- scores_long |>
  dplyr::inner_join(regimes, by = c("origin_date", "location"))

all_dates <- sort(unique(scores_long$origin_date))

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
# LOO-CV THREE-METHOD COMPARISON
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
  method  = names(loo_results),
  loo_wis = unlist(loo_results),
  ranking = rank(unlist(loo_results), ties.method = "min")
) |>
  dplyr::arrange(loo_wis)

readr::write_csv(comparison,
                 file.path(PHASE4_DIR, "method_loo_comparison_expanded.csv"))

cat("\n  LOO-CV ranking:\n")
print(comparison)

# Compare to Phase 3A's primary ensemble LOO results.
phase3_comparison_path <- file.path(HUB_PATH, "analysis", "phase3",
                                     "method_loo_comparison.csv")
if (file.exists(phase3_comparison_path)) {
  phase3_comp <- readr::read_csv(phase3_comparison_path, show_col_types = FALSE)
  cat("\n  Phase 3A (primary ensemble, 8 models) for comparison:\n")
  print(phase3_comp)

  # Compute primary-vs-expanded delta per method.
  cat("\n  Per-method delta (expanded - primary):\n")
  joined <- comparison |>
    dplyr::select(method, expanded_loo_wis = loo_wis) |>
    dplyr::inner_join(
      phase3_comp |> dplyr::select(method, primary_loo_wis = loo_wis),
      by = "method"
    ) |>
    dplyr::mutate(
      delta_wis     = expanded_loo_wis - primary_loo_wis,
      pct_change    = round(100 * delta_wis / primary_loo_wis, 2)
    ) |>
    dplyr::select(method, primary_loo_wis, expanded_loo_wis,
                  delta_wis, pct_change)
  print(joined)
}


# ============================================================================
# COMPUTE FULL-VALIDATION WEIGHTS
# ============================================================================
cat("\n--- Computing full-validation weights ---\n")

w_inv   <- weights_inverse_wis(scores_with_regime)
w_bma   <- weights_log_score_bma(scores_with_regime, tau = best_tau)
w_stack <- weights_stacking(scores_with_regime, quantiles)

readr::write_csv(w_inv,   file.path(PHASE4_DIR, "weights_inverse_wis_expanded.csv"))
readr::write_csv(w_bma,   file.path(PHASE4_DIR, "weights_log_score_bma_expanded.csv"))
readr::write_csv(w_stack, file.path(PHASE4_DIR, "weights_stacking_expanded.csv"))

# Print weight tables labelled by short model_label for readability.
print_weights <- function(w, label) {
  cat(sprintf("\n  %s weights:\n", label))
  wide <- w |>
    dplyr::mutate(model_label = MODEL_LABELS[model]) |>
    dplyr::select(model_label, regime, weight) |>
    tidyr::pivot_wider(names_from = regime, values_from = weight) |>
    dplyr::arrange(model_label)
  print(wide)
}

print_weights(w_inv,   "Inverse-WIS")
print_weights(w_bma,   sprintf("Log-score BMA (tau = %.3f)", best_tau))
print_weights(w_stack, "Stacking (QP)")


# ============================================================================
# DESIGNATE WINNER AND SAVE PRODUCTION WEIGHTS
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

readr::write_csv(production_weights,
                 file.path(PHASE4_DIR, "production_weights_expanded.csv"))
cat(sprintf("  Wrote: %s\n",
            file.path(PHASE4_DIR, "production_weights_expanded.csv")))


# ============================================================================
# METADATA
# ============================================================================
metadata_path <- file.path(PHASE4_DIR, "bma_weights_expanded_metadata.txt")
sink(metadata_path)
cat("Phase 4C: BMA weight comparison (expanded ensemble)\n")
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("scoringutils: %s\n", packageVersion("scoringutils")))
cat(sprintf("quadprog: %s\n", packageVersion("quadprog")))
cat(sprintf("\nEnsemble models (%d):\n", length(ENSEMBLE_MODELS)))
cat(paste("  -", ENSEMBLE_MODELS, " (", MODEL_LABELS, ")",
          sep = "", collapse = "\n"))
cat(sprintf("\n\nLOO-CV results:\n"))
cat(paste(capture.output(print(comparison)), collapse = "\n"))
cat(sprintf("\n\nWinner: %s\n", winner))
cat(sprintf("Best tau (log-score BMA): %.3f\n", best_tau))
cat(sprintf("Inverse-WIS beta: %.3f\n", INV_WIS_BETA))
sink()


cat("\n========================================\n")
cat("  Phase 4C complete: expanded BMA weights computed\n")
cat("========================================\n")
cat(sprintf("  Outputs in %s:\n", PHASE4_DIR))
cat("    - method_loo_comparison_expanded.csv\n")
cat("    - weights_inverse_wis_expanded.csv\n")
cat("    - weights_log_score_bma_expanded.csv\n")
cat("    - weights_stacking_expanded.csv\n")
cat("    - production_weights_expanded.csv\n")
cat("    - bma_weights_expanded_metadata.txt\n")
cat("========================================\n")
