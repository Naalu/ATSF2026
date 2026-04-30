# =============================================================================
# scripts/reproduce_full.R
#
# Driver — Path A: Full from-scratch reproduction
#
# Runs the entire FluCast pipeline from raw data to the final BMA-expanded
# ensemble submission. This script demonstrates that the project is fully
# reproducible from inputs alone, with no manual intervention.
#
# WARNING — long runtime expected:
#   * Phase 2 candidate scoring:        ~30 minutes (per candidate)
#   * Phase 2 bsts validation:           ~4.6 hours
#   * Phase 4 expanded candidate scoring: ~1 hour
#   * Phase 4 ensemble matrix generation: ~30 minutes
#   Total: ~6-12 hours depending on machine
#
# Use scripts/reproduce_smart.R if you want a checkpointed version that
# skips already-completed phases.
#
# Use scripts/reproduce_final.R if you only want to regenerate the final
# submission CSVs from existing intermediate artifacts (~30 seconds).
#
# REQUIREMENTS:
#   - R >= 4.3 with packages: tidyverse, fpp3, scoringutils, bsts,
#     hubValidations, fable, fabletools, tsibble
#   - Working directory: project root
#   - Inputs:
#       target-data/oracle-output.csv
#       hub-config/tasks.json
#       model-output/<other-team>/...     (classmate models for expanded pool)
#
# OUTPUTS at completion:
#   - model-output/KReger-bma_ensemble/  (134 CSVs — submission)
#   - model-metadata/KReger-bma_ensemble.yml
#   - analysis/phase[2-5]/...             (intermediate artifacts)
#
# Run from project root:
#   Rscript scripts/reproduce_full.R
# =============================================================================

# Track total runtime
.start_time <- Sys.time()

cat("\n",
    "===========================================================\n",
    " FluCast — Full Reproduction Pipeline                      \n",
    " Started: ", format(.start_time, "%Y-%m-%d %H:%M:%S"), "\n",
    "===========================================================\n",
    sep = "")

# Sanity check — must run from project root
if (!file.exists("hub-config/tasks.json")) {
  stop("This script must be run from the project root (where hub-config/ lives).")
}

# Helper: run a phase script with timing
run_phase <- function(phase_label, script_path) {
  cat(sprintf("\n>>> %s — %s <<<\n", phase_label, script_path))
  cat(sprintf("    started: %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  t0 <- Sys.time()
  tryCatch(
    source(script_path, echo = FALSE),
    error = function(e) {
      cat(sprintf("    FAILED: %s\n", e$message))
      stop(sprintf("Phase %s failed. See error above.", phase_label))
    }
  )
  t1 <- Sys.time()
  
  elapsed <- difftime(t1, t0, units = "mins")
  cat(sprintf("    elapsed: %.1f minutes\n", as.numeric(elapsed)))
}

# -----------------------------------------------------------------------------
# PHASE 1 — Baseline submission (KReger-snaive_bc_bs)
# -----------------------------------------------------------------------------
# Generates the original assigned-model baseline. This is what was submitted
# in the early proposal phase before ensemble work began.
# -----------------------------------------------------------------------------
run_phase("Phase 1A — Baseline validation forecasts", "scripts/run_validation.R")
run_phase("Phase 1B — Baseline test forecasts",       "scripts/run_test.R")

# -----------------------------------------------------------------------------
# PHASE 2 — Candidate model generation and scoring
# -----------------------------------------------------------------------------
# Generates 11 candidate models (KReger family), scores each on validation,
# computes correlation matrices and twin-pair WIS dominance for diversity
# pruning. This is the foundation of the primary 8-model pool.
# -----------------------------------------------------------------------------
run_phase("Phase 2A — Generate all candidate forecasts",
          "scripts/run_all_candidates.R")
run_phase("Phase 2B — bsts candidate (long-running)",
          "scripts/run_bsts_validation.R")
run_phase("Phase 2C — Score candidates on validation",
          "scripts/score_candidates.R")
run_phase("Phase 2D — Analyze candidates (correlations, regimes)",
          "scripts/analyze_candidates.R")
run_phase("Phase 2E — Twin-pair diagnostics (primary 8-model selection)",
          "scripts/twin_diagnostics.R")

# -----------------------------------------------------------------------------
# PHASE 3 — Primary pool BMA weights and ensemble
# -----------------------------------------------------------------------------
# Computes regime-conditional BMA weights on the 8-model primary pool and
# generates the corresponding ensemble. This was the original Phase 3
# submission scope.
# -----------------------------------------------------------------------------
run_phase("Phase 3A — Primary pool BMA weights",
          "scripts/compute_bma_weights.R")
run_phase("Phase 3B — Primary pool ensemble generation",
          "scripts/generate_ensemble.R")

# -----------------------------------------------------------------------------
# PHASE 4 — Expanded pool (with classmate submissions and delphi-epicast)
# -----------------------------------------------------------------------------
# Pulls 9 additional models from the FluSight Hub, re-scores, re-prunes for
# diversity, and computes new BMA weights. Generates all four ensembles
# (BMA × {primary, expanded} and Equal × {primary, expanded}) for the
# four-cell comparison matrix that drives the report's central finding.
# -----------------------------------------------------------------------------
run_phase("Phase 4A — Score expanded candidate pool (20 models)",
          "scripts/score_candidates_expanded.R")
run_phase("Phase 4B — Twin-pair diagnostics (expanded 10-model selection)",
          "scripts/twin_diagnostics_expanded.R")
run_phase("Phase 4C — Expanded pool BMA weights",
          "scripts/compute_bma_weights_expanded.R")
run_phase("Phase 4D — Generate all four ensembles + LOO comparison",
          "scripts/generate_ensembles_matrix.R")

# -----------------------------------------------------------------------------
# PHASE 5 — Promotion and submission
# -----------------------------------------------------------------------------
# Copies the winning BMA-expanded ensemble (LOO total WIS = 614.0) into
# model-output/ with submission-ready filenames and writes the metadata YAML.
# Validates each submission CSV via hubValidations.
# -----------------------------------------------------------------------------
run_phase("Phase 5A — Promote BMA-expanded to model-output/",
          "scripts/promote_ensemble.R")

# -----------------------------------------------------------------------------
# PHASE 5D-E — Diagnostics (validation coverage + test scoring)
# -----------------------------------------------------------------------------
# These do not affect the submission but produce supporting evidence:
# prediction interval calibration on validation, leaderboard against the
# rest of the class on test, and a four-cell test-set evaluation.
# -----------------------------------------------------------------------------
run_phase("Phase 5D — Test-set evaluation (4-cell)",
          "scripts/score_test_set.R")
run_phase("Phase 5D — Test-set leaderboard (vs all hub models)",
          "scripts/score_test_leaderboard.R")
run_phase("Phase 5E — Validation prediction interval coverage",
          "scripts/score_validation_coverage.R")

# -----------------------------------------------------------------------------
# Done
# -----------------------------------------------------------------------------
.end_time <- Sys.time()
.total_elapsed <- difftime(.end_time, .start_time, units = "hours")

cat("\n",
    "===========================================================\n",
    " FluCast pipeline complete                                 \n",
    " Finished: ", format(.end_time, "%Y-%m-%d %H:%M:%S"), "\n",
    " Total runtime: ", sprintf("%.2f hours", as.numeric(.total_elapsed)), "\n",
    "===========================================================\n",
    sep = "")

cat("\nSubmission ready at:\n")
cat("  model-output/KReger-bma_ensemble/   (134 CSVs)\n")
cat("  model-metadata/KReger-bma_ensemble.yml\n\n")
