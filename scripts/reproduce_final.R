# =============================================================================
# scripts/reproduce_final.R
#
# Driver — Path C: Final-step regeneration
#
# Re-runs only the final assembly step: takes existing intermediate artifacts
# (component forecasts, BMA weights, ensemble configs) and produces the final
# submission CSVs. Useful for:
#
#   - Spot-checking the promotion step works
#   - Regenerating the submission after a small fix without re-scoring
#   - Sanity-checking that intermediate artifacts produce the expected output
#
# Runtime: ~30 seconds.
#
# REQUIRED ARTIFACTS (from earlier phases — script will fail clearly if any
# are missing):
#   - analysis/phase4/staging/ensemble_bma_expanded/  (134 staged CSVs)
#   - model-metadata/KReger-bma_ensemble.yml          (or it'll be regenerated)
#
# OUTPUTS:
#   - model-output/KReger-bma_ensemble/   (134 CSVs)
#   - model-metadata/KReger-bma_ensemble.yml (overwritten if exists)
#
# Run from project root:
#   Rscript scripts/reproduce_final.R
# =============================================================================

.start_time <- Sys.time()

cat("\n",
    "===========================================================\n",
    " FluCast — Final-Step Reproduction                         \n",
    " Started: ", format(.start_time, "%Y-%m-%d %H:%M:%S"), "\n",
    "===========================================================\n",
    sep = "")

if (!file.exists("hub-config/tasks.json")) {
  stop("This script must be run from the project root.")
}

# Check required input artifacts before starting
required <- list(
  "Phase 4D staged BMA-expanded forecasts" =
    "analysis/phase4/staging/ensemble_bma_expanded",
  "Phase 4 production weights" =
    "analysis/phase4/production_weights_expanded.csv"
)

missing <- character()
for (label in names(required)) {
  path <- required[[label]]
  if (!file.exists(path)) {
    missing <- c(missing, sprintf("%s (%s)", label, path))
  }
}

if (length(missing) > 0L) {
  cat("\nMISSING REQUIRED ARTIFACTS:\n")
  for (m in missing) cat("  - ", m, "\n", sep = "")
  cat("\nUse scripts/reproduce_full.R or scripts/reproduce_smart.R to ",
      "generate them first.\n", sep = "")
  stop("Required artifacts not present.")
}

# All required artifacts exist — proceed with promotion
cat("\nAll required artifacts present. Running promotion...\n\n")

t0 <- Sys.time()
source("scripts/promote_ensemble.R", echo = FALSE)
t1 <- Sys.time()

cat(sprintf("\nPromotion elapsed: %.1f seconds\n",
            as.numeric(difftime(t1, t0, units = "secs"))))

# Verify expected output count
n_csvs <- length(list.files("model-output/KReger-bma_ensemble",
                             pattern = "\\.csv$"))
cat(sprintf("\nFinal submission: %d CSVs in model-output/KReger-bma_ensemble/\n",
            n_csvs))

if (n_csvs == 134L) {
  cat("Expected count (134) matches.\n")
} else {
  cat(sprintf("WARNING: expected 134 CSVs, found %d.\n", n_csvs))
}

.end_time <- Sys.time()
cat(sprintf("\nTotal runtime: %.1f seconds\n",
            as.numeric(difftime(.end_time, .start_time, units = "secs"))))
cat("Done.\n\n")
