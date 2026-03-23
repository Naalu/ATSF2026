###############################################################################
# FILE: scripts/analyze_candidates.R
#
# Error Correlation Analysis for Ensemble Candidate Selection
#
# This script reads the forecast CSVs from all candidate models, computes
# forecast errors (predicted median minus actual), and produces:
#   1. A pairwise error correlation matrix with heatmap visualization
#   2. Individual model WIS scores
#   3. A recommended ensemble composition based on diversity + skill
#
# The goal is to select 4-5 models whose errors are UNCORRELATED -- when
# Model A gets it wrong, Model B gets it right (or at least gets it wrong
# in a different direction). This error decorrelation is what makes
# ensembling effective.
#
# USAGE:
#   1. Run all candidate models first (see run_all_candidates.R)
#   2. source("scripts/analyze_candidates.R")
#
# OUTPUTS:
#   - Console: correlation matrix, WIS table, recommended ensemble
#   - Plots: correlation heatmap, individual model WIS comparison
#   - Saved: analysis results as RDS for use in the presentation
#
# Dependencies: tidyverse, fpp3 (for tsibble), ggplot2, patchwork
#
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(fpp3)
})

cat("==========================================================\n")
cat("  CANDIDATE MODEL ANALYSIS: Error Correlation & WIS\n")
cat("==========================================================\n\n")

# ---------------------------------------------------------------------------
# CONFIGURATION -- adjust these paths to match your setup
# ---------------------------------------------------------------------------

# Path to the ATSF2026 repo root
hub_path <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"

# Team abbreviation
team_abbr <- "KReger"

# All candidate model abbreviations (must match configs)
candidate_models <- c(
  "snaive_bc_bs",
  "ets_log",
  "arima_bc_bs",
  "nnetar_bc_bs",
  "theta_bc_bs",
  "tslm_fourier",
  "stl_arima_bc",
  "hist_week"
)

# Short display names for plots (same order as candidate_models)
display_names <- c(
  "SNAIVE+BC",
  "ETS+log",
  "ARIMA+BC",
  "NNETAR+BC",
  "THETA+BC",
  "TSLM+Fourier",
  "STL+ARIMA",
  "HistWeek"
)

# ---------------------------------------------------------------------------
# STEP 1: LOAD ACTUALS
# ---------------------------------------------------------------------------

cat("--- Loading actuals ---\n")
actuals <- read_csv(
  file.path(hub_path, "target-data", "time-series.csv"),
  show_col_types = FALSE
) |>
  select(location, target_end_date, actual = observation)

cat(sprintf("  %d actual observations loaded\n\n", nrow(actuals)))

# ---------------------------------------------------------------------------
# STEP 2: LOAD ALL CANDIDATE FORECASTS
# ---------------------------------------------------------------------------

cat("--- Loading candidate forecasts ---\n")

# Read all CSVs for each model and extract the horizon-1 median (q0.50)
all_forecasts <- list()

for (i in seq_along(candidate_models)) {
  model_abbr <- candidate_models[i]
  model_id   <- paste0(team_abbr, "-", model_abbr)
  model_dir  <- file.path(hub_path, "model-output", model_id)

  if (!dir.exists(model_dir)) {
    cat(sprintf("  SKIPPED %s: directory not found\n", model_id))
    next
  }

  csv_files <- list.files(model_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    cat(sprintf("  SKIPPED %s: no CSV files\n", model_id))
    next
  }

  # Read all CSVs for this model, extracting the median (q0.50) forecast
  model_df <- map_dfr(csv_files, function(f) {
    read_csv(f, show_col_types = FALSE)
  }) |>
    filter(
      output_type_id == 0.5,    # Median (q0.50)
      horizon == 1               # 1-week-ahead only (most comparable)
    ) |>
    mutate(model = display_names[i]) |>
    select(model, origin_date, location, target_end_date, forecast = value)

  all_forecasts[[model_id]] <- model_df
  cat(sprintf("  %s: %d forecasts from %d files\n",
              display_names[i], nrow(model_df), length(csv_files)))
}

# Combine into one dataframe
forecasts_df <- bind_rows(all_forecasts)

if (nrow(forecasts_df) == 0) {
  stop("No forecast data found. Run the candidate models first.")
}

cat(sprintf("\n  Total: %d forecasts across %d models\n\n",
            nrow(forecasts_df), n_distinct(forecasts_df$model)))

# ---------------------------------------------------------------------------
# STEP 3: JOIN WITH ACTUALS AND COMPUTE ERRORS
# ---------------------------------------------------------------------------

cat("--- Computing forecast errors ---\n")

errors_df <- forecasts_df |>
  inner_join(actuals, by = c("location", "target_end_date")) |>
  mutate(error = forecast - actual)   # Positive = overprediction

n_models <- n_distinct(errors_df$model)
cat(sprintf("  %d forecast-actual pairs across %d models\n\n", nrow(errors_df), n_models))

# ---------------------------------------------------------------------------
# STEP 4: PAIRWISE ERROR CORRELATION MATRIX
# ---------------------------------------------------------------------------

cat("--- Computing pairwise error correlations ---\n")

# Pivot so each model's errors are in a separate column
errors_wide <- errors_df |>
  select(model, origin_date, location, error) |>
  # Create a unique row ID from origin_date + location
  unite("row_id", origin_date, location, sep = "_") |>
  pivot_wider(names_from = model, values_from = error)

# Compute Pearson correlation matrix (only on model columns)
model_cols <- intersect(display_names, names(errors_wide))
cor_matrix <- cor(
  errors_wide |> select(all_of(model_cols)),
  use = "pairwise.complete.obs"
)

# Print the correlation matrix
cat("\n  Pairwise Error Correlation Matrix:\n")
cat("  (Lower values = more diversity = better for ensembling)\n\n")
print(round(cor_matrix, 3))

# ---------------------------------------------------------------------------
# STEP 5: CORRELATION HEATMAP
# ---------------------------------------------------------------------------

cat("\n--- Generating correlation heatmap ---\n")

# Convert correlation matrix to long format for ggplot
cor_long <- as.data.frame(as.table(cor_matrix)) |>
  rename(Model1 = Var1, Model2 = Var2, Correlation = Freq)

p_heatmap <- ggplot(cor_long, aes(x = Model1, y = Model2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0.5, limits = c(0, 1),
    name = "Error\nCorrelation"
  ) +
  labs(
    title = "Pairwise Forecast Error Correlations (Horizon 1)",
    subtitle = "Lower correlation = greater diversity = better ensemble component",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

print(p_heatmap)

# ---------------------------------------------------------------------------
# STEP 6: INDIVIDUAL MODEL WIS (simplified)
# ---------------------------------------------------------------------------
#
# WIS = (1/K) * sum_{k=1}^{K} IS_alpha_k(F, y)
# where IS_alpha is the interval score at level alpha.
#
# For a simpler analysis, we compute MAE of the median (the "sharpness"
# component of WIS) and the RMSE. These are strong proxies for WIS when
# all models use the same quantile structure.
#
# For a full WIS calculation, use scoringutils::score() or hubEvals.
# ---------------------------------------------------------------------------

cat("\n--- Computing individual model accuracy ---\n")

accuracy_df <- errors_df |>
  group_by(model) |>
  summarise(
    n_forecasts = n(),
    MAE   = mean(abs(error), na.rm = TRUE),
    RMSE  = sqrt(mean(error^2, na.rm = TRUE)),
    Bias  = mean(error, na.rm = TRUE),        # Positive = systematic overprediction
    .groups = "drop"
  ) |>
  arrange(MAE)

cat("\n  Model Accuracy (Horizon 1, all locations, all dates):\n\n")
print(accuracy_df, n = 20)

# Bar chart of MAE by model
p_accuracy <- ggplot(accuracy_df, aes(x = reorder(model, MAE), y = MAE)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", MAE)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(
    title = "Median Absolute Error by Model (Horizon 1)",
    subtitle = "Lower = more accurate point forecasts",
    x = NULL, y = "MAE"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

print(p_accuracy)

# ---------------------------------------------------------------------------
# STEP 7: ENSEMBLE SELECTION ALGORITHM
# ---------------------------------------------------------------------------
#
# Greedy forward selection: start with the best individual model (lowest
# MAE), then iteratively add the model that has the lowest AVERAGE error
# correlation with models already selected. Stop when adding another
# model doesn't reduce the average pairwise correlation by at least 0.02.
# ---------------------------------------------------------------------------

cat("\n--- Ensemble selection (greedy diversity maximization) ---\n")

available <- model_cols
selected  <- character(0)

# Start with the best individual model by MAE
best_model <- accuracy_df$model[1]
if (best_model %in% available) {
  selected  <- c(selected, best_model)
  available <- setdiff(available, best_model)
}

cat(sprintf("\n  Seed model (best MAE): %s\n", best_model))

# Greedy selection: add model with lowest avg correlation to selected set
while (length(available) > 0 && length(selected) < 6) {

  # For each candidate, compute its average correlation with selected models
  avg_cors <- sapply(available, function(cand) {
    mean(cor_matrix[cand, selected], na.rm = TRUE)
  })

  # Pick the candidate with the lowest average correlation
  best_cand <- names(which.min(avg_cors))
  best_cor  <- min(avg_cors)

  # Stop if adding this model doesn't bring enough diversity
  if (best_cor > 0.80 && length(selected) >= 4) {
    cat(sprintf("  Stopping: next best candidate (%s) has avg correlation %.3f\n",
                best_cand, best_cor))
    break
  }

  selected  <- c(selected, best_cand)
  available <- setdiff(available, best_cand)
  cat(sprintf("  +%s (avg correlation with selected: %.3f)\n", best_cand, best_cor))
}

cat(sprintf("\n  RECOMMENDED ENSEMBLE (%d models):\n", length(selected)))
for (s in selected) {
  # Find its MAE rank
  rank <- which(accuracy_df$model == s)
  cat(sprintf("    - %s (MAE rank: %d/%d)\n", s, rank, nrow(accuracy_df)))
}

# ---------------------------------------------------------------------------
# STEP 8: SELECTED ENSEMBLE CORRELATION SUB-MATRIX
# ---------------------------------------------------------------------------

cat("\n  Correlation matrix for selected ensemble:\n\n")
sub_cor <- cor_matrix[selected, selected]
print(round(sub_cor, 3))

avg_off_diag <- mean(sub_cor[lower.tri(sub_cor)])
cat(sprintf("\n  Average off-diagonal correlation: %.3f\n", avg_off_diag))
cat("  (Target: below 0.70 for effective ensembling)\n")

# ---------------------------------------------------------------------------
# STEP 9: SAVE RESULTS
# ---------------------------------------------------------------------------

analysis_results <- list(
  errors_df      = errors_df,
  cor_matrix     = cor_matrix,
  accuracy_df    = accuracy_df,
  selected       = selected,
  selected_cor   = sub_cor,
  timestamp      = Sys.time()
)

output_path <- file.path(hub_path, "scripts", "candidate_analysis.rds")
saveRDS(analysis_results, output_path)
cat(sprintf("\n  Results saved to: %s\n", output_path))

# ---------------------------------------------------------------------------
# STEP 10: REGIME-SPECIFIC ANALYSIS (bonus for presentation)
# ---------------------------------------------------------------------------

cat("\n--- Regime-specific error correlation ---\n")

# Split into "peak" (ILI > location median) and "trough" regimes
# This reveals whether models that are correlated overall diverge at peaks

regime_df <- errors_df |>
  inner_join(actuals, by = c("location", "target_end_date")) |>
  group_by(location) |>
  mutate(
    median_ili = median(actual, na.rm = TRUE),
    regime = if_else(actual > median_ili, "Peak", "Trough")
  ) |>
  ungroup()

for (r in c("Peak", "Trough")) {
  cat(sprintf("\n  Error correlations during %s regime:\n", r))

  regime_wide <- regime_df |>
    filter(regime == r) |>
    select(model, origin_date, location, error) |>
    unite("row_id", origin_date, location, sep = "_") |>
    pivot_wider(names_from = model, values_from = error)

  regime_cols <- intersect(selected, names(regime_wide))
  if (length(regime_cols) >= 2) {
    regime_cor <- cor(
      regime_wide |> select(all_of(regime_cols)),
      use = "pairwise.complete.obs"
    )
    print(round(regime_cor, 3))
  }
}

cat("\n==========================================================\n")
cat("  ANALYSIS COMPLETE\n")
cat("==========================================================\n")
cat("\n  Key outputs for your presentation:\n")
cat("  1. Correlation heatmap (currently displayed)\n")
cat("  2. MAE comparison bar chart (currently displayed)\n")
cat("  3. Recommended ensemble composition (printed above)\n")
cat("  4. Regime-specific correlations (printed above)\n")
cat(sprintf("  5. Full results saved to: %s\n", output_path))
cat("\n  Use ggsave() to save plots for your slides.\n")
