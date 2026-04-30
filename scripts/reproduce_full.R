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
#   - Working directory: project root (drivers verify this)
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

# Sanity check — must run from project root (paths in sourced scripts
# resolve to getwd(), so this matters)
if (!file.exists("hub-config/tasks.json")) {
  stop("This script must be run from the project root (where hub-config/ lives).")
}

# Helper: run a phase script with timing
# config_path: optional config to pass via CONFIG_PATH variable (used by
# Phase 1 scripts run_validation.R and run_test.R, which require a config)
run_phase <- function(phase_label, script_path, config_path = NULL) {
  cat(sprintf("\n>>> %s — %s <<<\n", phase_label, script_path))
  cat(sprintf("    started: %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  # Set CONFIG_PATH in the global env if requested. Phase 1 scripts look
  # for this variable and error out if it's not set.
  if (!is.null(config_path)) {
    assign("CONFIG_PATH", config_path, envir = .GlobalEnv)
    cat(sprintf("    config: %s\n", config_path))
  }
  
  t0 <- Sys.time()
  tryCatch(
    source(script_path, echo = FALSE),
    error = function(e) {
      cat(sprintf("    FAILED: %s\n", e$message))
      stop(sprintf("Phase %s failed. See error above.", phase_label))
    },
    finally = {
      # Clear CONFIG_PATH after the phase to avoid contaminating later phases
      if (!is.null(config_path) && exists("CONFIG_PATH", envir = .GlobalEnv)) {
        rm("CONFIG_PATH", envir = .GlobalEnv)
      }
    }
  )
  t1 <- Sys.time()
  
  elapsed <- difftime(t1, t0, units = "mins")
  cat(sprintf("    elapsed: %.1f minutes\n", as.numeric(elapsed)))
}

# -----------------------------------------------------------------------------
# PHASE 1 — Baseline submission (KReger-snaive_bc_bs)
# -----------------------------------------------------------------------------
# Generates the original assigned-model baseline using the generic Phase 1
# drivers. Both run_validation.R and run_test.R need a CONFIG_PATH set
# (the engine is generic — same drivers, different configs run different models).
# -----------------------------------------------------------------------------
run_phase("Phase 1A — Baseline validation forecasts",
          "scripts/run_validation.R",
          config_path = "configs/snaive_bc_bs.R")
run_phase("Phase 1B — Baseline test forecasts",
          "scripts/run_test.R",
          config_path = "configs/snaive_bc_bs.R")

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
          "scripts/archive/score_candidates.R")
run_phase("Phase 2D — Analyze candidates (correlations, regimes)",
          "scripts/analyze_candidates.R")
run_phase("Phase 2E — Twin-pair diagnostics (primary 8-model selection)",
          "scripts/archive/twin_diagnostics.R")

# -----------------------------------------------------------------------------
# PHASE 3 — Primary pool BMA weights and ensemble
# -----------------------------------------------------------------------------
# Computes regime-conditional BMA weights on the 8-model primary pool and
# generates the corresponding ensemble. This was the original Phase 3
# submission scope.
# -----------------------------------------------------------------------------
run_phase("Phase 3A — Primary pool BMA weights",
          "scripts/archive/compute_bma_weights.R")
run_phase("Phase 3B — Primary pool ensemble generation",
          "scripts/archive/generate_ensemble.R")

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
# PHASE 5E — Validation prediction interval coverage
# -----------------------------------------------------------------------------
# Validation-only diagnostic: 50% / 95% empirical PI coverage for each of the
# four ensembles. Supports the report's PIC discussion in §3.
#
# Note: test-period evaluation is not run from this driver. Test-set scoring
# is reserved for the instructor.
# -----------------------------------------------------------------------------
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
