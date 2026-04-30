# =============================================================================
# scripts/reproduce_smart.R
#
# Driver — Path B: Checkpointed reproduction
#
# Same pipeline as reproduce_full.R, but each phase first checks whether its
# expected output artifacts already exist. If they do, the phase is skipped
# with an informative message. Useful for:
#
#   - Re-running after a single phase has been modified
#   - Resuming after an interrupted full reproduction
#   - Spot-checking phases without paying the full pipeline cost
#
# A phase is considered "done" if its sentinel output file exists. The
# user can force a full re-run of any phase by deleting that file.
#
# Use scripts/reproduce_full.R for unconditional from-scratch reproduction.
# Use scripts/reproduce_final.R to skip directly to the final assembly step.
#
# Run from project root:
#   Rscript scripts/reproduce_smart.R
#
# To force a phase to re-run, delete its sentinel before running:
#   rm analysis/phase4/staging/ensemble_bma_expanded/2015-10-24-KReger-bma_expanded.csv
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

# Phase definitions: label, script, sentinel output (existence = done).
# When designing sentinels, prefer files that are produced LAST by the phase,
# so partial completion isn't mistaken for full completion.
phases <- list(
  list(label = "Phase 1A — Baseline validation",
       script = "scripts/run_validation.R",
       sentinel = "analysis/phase1/snaive_bc_bs/validation/2017-05-06.csv"),
  list(label = "Phase 1B — Baseline test",
       script = "scripts/run_test.R",
       sentinel = "model-output/KReger-snaive_bc_bs/2020-02-29-KReger-snaive_bc_bs.csv"),
  list(label = "Phase 2A — All candidate forecasts",
       script = "scripts/run_all_candidates.R",
       sentinel = "analysis/phase2/forecasts/nnetar_bc_bs/2017-05-06.csv"),
  list(label = "Phase 2B — bsts validation (long)",
       script = "scripts/run_bsts_validation.R",
       sentinel = "analysis/phase2/forecasts/bsts_seasonal/2017-05-06.csv"),
  list(label = "Phase 2C — Score candidates",
       script = "scripts/archive/score_candidates.R",
       sentinel = "analysis/phase2/wis_summary.csv"),
  list(label = "Phase 2D — Analyze candidates",
       script = "scripts/analyze_candidates.R",
       sentinel = "analysis/phase2/regimes.csv"),
  list(label = "Phase 2E — Twin diagnostics (primary)",
       script = "scripts/archive/twin_diagnostics.R",
       sentinel = "analysis/phase2/shortlist.csv"),
  list(label = "Phase 3A — Primary BMA weights",
       script = "scripts/archive/compute_bma_weights.R",
       sentinel = "analysis/phase3/weights_log_score_bma.csv"),
  list(label = "Phase 3B — Primary ensemble",
       script = "scripts/archive/generate_ensemble.R",
       sentinel = "analysis/phase3/staging/ensemble_bma_primary/2017-05-06-KReger-bma_primary.csv"),
  list(label = "Phase 4A — Score expanded candidates",
       script = "scripts/score_candidates_expanded.R",
       sentinel = "analysis/phase4/wis_summary_expanded.csv"),
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
  
  t0 <- Sys.time()
  tryCatch(
    source(p$script, echo = FALSE),
    error = function(e) {
      cat(sprintf("    FAILED: %s\n", e$message))
      stop(sprintf("Phase failed: %s", p$label))
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
