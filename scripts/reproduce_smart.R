# =============================================================================
# scripts/reproduce_smart.R
#
# Driver — Path B: Checkpointed reproduction
#
# Same pipeline as reproduce_full.R, but each phase first checks whether
# evidence of its successful completion exists. If it does, the phase is
# skipped with an informative message.
#
# Design note on sentinels:
#   A "sentinel" file is one whose presence means a phase has run
#   successfully at some point. For most phases, the natural sentinel is
#   the phase's primary output. But for some phases (e.g., Phase 2A and 2B
#   producing intermediate forecast CSVs), the natural output may have
#   been consumed and overwritten by a later phase, leaving no trace.
#
#   For those phases, this script uses a DOWNSTREAM sentinel: a file
#   produced by a later phase that wouldn't exist if the earlier phase
#   hadn't run. This avoids redundant re-runs of work that's clearly
#   already been completed.
#
#   Spencer can force a single phase to re-run by deleting its sentinel
#   AND ensuring the downstream phase's outputs are also gone (otherwise
#   the script will think the work is done).
#
# Use scripts/reproduce_full.R for unconditional from-scratch reproduction.
# Use scripts/reproduce_final.R to skip directly to the final assembly step.
#
# Run from project root:
#   Rscript scripts/reproduce_smart.R
# =============================================================================

.start_time <- Sys.time()

cat("\n",
    "===========================================================\n",
    " FluCast — Smart (Checkpointed) Reproduction Pipeline      \n",
    " Started: ", format(.start_time, "%Y-%m-%d %H:%M:%S"), "\n",
    "===========================================================\n",
    sep = "")

if (!file.exists("hub-config/tasks.json")) {
  stop("This script must be run from the project root.")
}

# Phase definitions: label, script, sentinel, optional config_path.
#
# Sentinels point to artifacts whose existence is sufficient evidence that
# the phase has been completed. For chained phases, the sentinel may be
# a downstream artifact rather than the phase's own immediate output.
phases <- list(
  list(label = "Phase 1A — Baseline validation",
       script = "scripts/run_validation.R",
       # Downstream sentinel: the primary pool ensemble inputs depend on this
       sentinel = "model-output/KReger-snaive_bc_bs/2017-05-06-KReger-snaive_bc_bs.csv",
       config_path = "configs/snaive_bc_bs.R"),
  list(label = "Phase 1B — Baseline test",
       script = "scripts/run_test.R",
       sentinel = "model-output/KReger-snaive_bc_bs/2020-02-29-KReger-snaive_bc_bs.csv",
       config_path = "configs/snaive_bc_bs.R"),
  list(label = "Phase 2A — All candidate forecasts",
       script = "scripts/run_all_candidates.R",
       # Downstream sentinel: the expanded shortlist requires all candidate
       # forecasts to have been generated and scored
       sentinel = "analysis/phase4/shortlist_expanded.csv"),
  list(label = "Phase 2B — bsts validation (long)",
       script = "scripts/run_bsts_validation.R",
       # bsts produces forecasts to model-output/KReger-bsts_seasonal/.
       # Use that location as the sentinel (validation date 2017-05-06 is
       # the last validation date — its presence confirms a complete run).
       sentinel = "model-output/KReger-bsts_seasonal/2017-05-06-KReger-bsts_seasonal.csv"),
  list(label = "Phase 2C — Score candidates",
       script = "scripts/archive/score_candidates.R",
       # Downstream sentinel
       sentinel = "analysis/phase2/regimes.csv"),
  list(label = "Phase 2D — Analyze candidates",
       script = "scripts/analyze_candidates.R",
       sentinel = "analysis/phase2/regimes.csv"),
  list(label = "Phase 2E — Twin diagnostics (primary)",
       script = "scripts/archive/twin_diagnostics.R",
       # Downstream: primary BMA weights require Phase 2E's shortlist
       sentinel = "analysis/phase3/weights_log_score_bma.csv"),
  list(label = "Phase 3A — Primary BMA weights",
       script = "scripts/archive/compute_bma_weights.R",
       sentinel = "analysis/phase3/weights_log_score_bma.csv"),
  list(label = "Phase 3B — Primary ensemble",
       script = "scripts/archive/generate_ensemble.R",
       # Downstream: Phase 4 matrix consumes primary ensemble outputs
       sentinel = "analysis/phase4/matrix_loo_comparison.csv"),
  list(label = "Phase 4A — Score expanded candidates",
       script = "scripts/score_candidates_expanded.R",
       # Downstream: shortlist_expanded requires Phase 4A's scoring
       sentinel = "analysis/phase4/shortlist_expanded.csv"),
  list(label = "Phase 4B — Twin diagnostics (expanded)",
       script = "scripts/twin_diagnostics_expanded.R",
       sentinel = "analysis/phase4/shortlist_expanded.csv"),
  list(label = "Phase 4C — Expanded BMA weights",
       script = "scripts/compute_bma_weights_expanded.R",
       sentinel = "analysis/phase4/weights_log_score_bma_expanded.csv"),
  list(label = "Phase 4D — Four-ensemble matrix",
       script = "scripts/generate_ensembles_matrix.R",
       sentinel = "analysis/phase4/matrix_loo_comparison.csv"),
  list(label = "Phase 5A — Promote to submission",
       script = "scripts/promote_ensemble.R",
       sentinel = "model-output/KReger-bma_ensemble/2020-02-29-KReger-bma_ensemble.csv"),
  list(label = "Phase 5E — Validation PI coverage",
       script = "scripts/score_validation_coverage.R",
       sentinel = "analysis/phase4/validation_coverage.csv")
)

n_run <- 0L
n_skip <- 0L

for (p in phases) {
  cat(sprintf("\n>>> %s\n", p$label))
  
  if (file.exists(p$sentinel)) {
    cat(sprintf("    SKIP (sentinel exists: %s)\n", p$sentinel))
    n_skip <- n_skip + 1L
    next
  }
  
  cat(sprintf("    RUN: %s\n", p$script))
  cat(sprintf("    started: %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  # Phase 1 scripts need CONFIG_PATH set in the global env
  if (!is.null(p$config_path)) {
    assign("CONFIG_PATH", p$config_path, envir = .GlobalEnv)
    cat(sprintf("    config: %s\n", p$config_path))
  }
  
  t0 <- Sys.time()
  tryCatch(
    source(p$script, echo = FALSE),
    error = function(e) {
      cat(sprintf("    FAILED: %s\n", e$message))
      stop(sprintf("Phase failed: %s", p$label))
    },
    finally = {
      if (!is.null(p$config_path) && exists("CONFIG_PATH", envir = .GlobalEnv)) {
        rm("CONFIG_PATH", envir = .GlobalEnv)
      }
    }
  )
  t1 <- Sys.time()
  
  cat(sprintf("    elapsed: %.1f minutes\n",
              as.numeric(difftime(t1, t0, units = "mins"))))
  n_run <- n_run + 1L
}

.end_time <- Sys.time()
cat("\n",
    "===========================================================\n",
    sprintf(" Pipeline complete                                         \n"),
    sprintf(" Phases run:  %d\n", n_run),
    sprintf(" Phases skipped (sentinels existed): %d\n", n_skip),
    sprintf(" Total runtime: %.2f minutes\n",
            as.numeric(difftime(.end_time, .start_time, units = "mins"))),
    "===========================================================\n",
    sep = "")
