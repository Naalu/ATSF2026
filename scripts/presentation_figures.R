###############################################################################
# FILE: scripts/presentation_figures.R
#
# Generates publication-quality figures and accurate metrics for the
# forecast proposal defense presentation. Saves all outputs as PNG files.
#
# OUTPUTS (saved to slides/ directory):
#   01_ili_timeseries.png        - ILI% time series with key features annotated
#   02_correlation_heatmap.png   - Pairwise error correlation matrix
#   03_model_accuracy.png        - MAE comparison across models
#   04_wis_comparison.png        - Full WIS across all quantiles
#   05_coverage_table.png        - 50% and 95% coverage by model
#   06_ensemble_selection.png    - Selected vs dropped models visualization
#   07_regime_correlation.png    - Peak vs trough error correlations
#   08_fan_chart_comparison.png  - Fan charts for 3 models on same date
#   09_error_timeseries.png      - Errors over time showing when models diverge
#
# USAGE:
#   source("scripts/presentation_figures.R")
#
# Dependencies: tidyverse, fpp3, patchwork, scales
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(fpp3)
  library(patchwork)
  library(scales)
})

source("scripts/forecast_helpers.R")

hub_path <- normalizePath(".", mustWork = FALSE)
team_abbr <- "KReger"

# Create output directory
slides_dir <- file.path(hub_path, "slides")
if (!dir.exists(slides_dir)) dir.create(slides_dir, recursive = TRUE)

# Model display names and colors (consistent across all plots)
model_info <- tribble(
  ~model_abbr,    ~display_name,   ~color,     ~include_ensemble,
  "snaive_bc_bs", "SNAIVE+BC",     "#E69F00",  TRUE,
  "ets_log",      "ETS+log",       "#56B4E9",  TRUE,
  "arima_bc_bs",  "ARIMA+BC",      "#009E73",  TRUE,
  "nnetar_bc_bs", "NNETAR+BC",     "#D55E00",  TRUE,
  "tslm_fourier", "TSLM+Fourier",  "#CC79A7",  FALSE,
  "stl_arima_bc", "STL+ARIMA",     "#0072B2",  FALSE,
  "hist_week",    "HistWeek",      "#F0E442",  TRUE
)

model_colors <- setNames(model_info$color, model_info$display_name)

cat("=== GENERATING PRESENTATION FIGURES ===\n\n")


# =========================================================================
# LOAD ALL DATA
# =========================================================================

cat("--- Loading data ---\n")

actuals <- read_csv(
  file.path(hub_path, "target-data", "time-series.csv"),
  show_col_types = FALSE
)

# Load ALL forecasts (all quantiles, all horizons) for WIS calculation
all_forecasts_full <- list()
for (i in seq_len(nrow(model_info))) {
  abbr <- model_info$model_abbr[i]
  dname <- model_info$display_name[i]
  model_id <- paste0(team_abbr, "-", abbr)
  model_dir <- file.path(hub_path, "model-output", model_id)

  if (!dir.exists(model_dir)) next
  csv_files <- list.files(model_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) next

  model_df <- map_dfr(csv_files, ~read_csv(.x, show_col_types = FALSE)) |>
    mutate(model = dname)
  all_forecasts_full[[dname]] <- model_df
  cat(sprintf("  %s: %d rows from %d files\n", dname, nrow(model_df), length(csv_files)))
}

forecasts_full <- bind_rows(all_forecasts_full)
cat(sprintf("  Total: %d rows across %d models\n\n", nrow(forecasts_full), n_distinct(forecasts_full$model)))


# =========================================================================
# FIGURE 1: ILI% Time Series with Key Features
# =========================================================================

cat("--- Figure 1: ILI time series ---\n")

p1 <- actuals |>
  filter(location == "US National") |>
  ggplot(aes(x = target_end_date, y = observation)) +
  geom_line(color = "#2C3E50", linewidth = 0.6) +
  geom_point(
    data = actuals |> filter(location == "US National", observation == 0),
    aes(x = target_end_date, y = observation),
    color = "#E74C3C", size = 1.5, alpha = 0.7
  ) +
  annotate("rect",
    xmin = as.Date("2017-12-01"), xmax = as.Date("2018-03-15"),
    ymin = -Inf, ymax = Inf,
    fill = "#E74C3C", alpha = 0.1
  ) +
  annotate("text",
    x = as.Date("2018-01-15"), y = 7,
    label = "2017-18\nSevere Season",
    size = 3, color = "#E74C3C", fontface = "bold"
  ) +
  annotate("text",
    x = as.Date("2016-07-01"), y = 0.5,
    label = "Zero-valued\nsummer weeks",
    size = 2.5, color = "#E74C3C"
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title = "US National Weighted ILI% (2003-2020)",
    subtitle = "Strong annual seasonality | Variance scales with level | Peak severity varies across seasons",
    x = NULL, y = "Weighted ILI %"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(slides_dir, "01_ili_timeseries.png"), p1,
       width = 10, height = 4.5, dpi = 200, bg = "white")
cat("  Saved 01_ili_timeseries.png\n")


# =========================================================================
# COMPUTE PROPER WIS
# =========================================================================

cat("\n--- Computing WIS ---\n")

# WIS = (1/K) * sum_{k=1}^{K} IS_{alpha_k}
# where IS_alpha = (upper - lower) + (2/alpha) * (lower - y) * I(y < lower)
#                                   + (2/alpha) * (y - upper) * I(y > upper)
#
# For quantile forecasts, WIS can also be computed as:
# WIS = (1/K) * sum_{k=1}^K 2 * pinball_loss(q_k, y, tau_k)
# where pinball_loss(q, y, tau) = (y - q) * tau if y >= q
#                                 (q - y) * (1 - tau) if y < q

wis_data <- forecasts_full |>
  filter(horizon == 1) |>
  select(model, origin_date, location, target_end_date,
         output_type_id, value) |>
  inner_join(
    actuals |> select(location, target_end_date, observation),
    by = c("location", "target_end_date")
  )

# Compute pinball loss for each quantile
wis_by_model <- wis_data |>
  mutate(
    tau = output_type_id,
    pinball = if_else(
      observation >= value,
      (observation - value) * tau,
      (value - observation) * (1 - tau)
    )
  ) |>
  # WIS = mean of 2 * pinball across all quantiles, then average across dates/locations
  group_by(model, origin_date, location) |>
  summarise(
    wis = mean(2 * pinball, na.rm = TRUE),
    .groups = "drop"
  )

# Overall WIS by model
wis_summary <- wis_by_model |>
  group_by(model) |>
  summarise(
    mean_wis = mean(wis, na.rm = TRUE),
    median_wis = median(wis, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(mean_wis)

cat("\n  Weighted Interval Score (Horizon 1):\n")
print(wis_summary, n = 10)

# Compute coverage
coverage_data <- forecasts_full |>
  filter(horizon == 1) |>
  select(model, origin_date, location, target_end_date, output_type_id, value) |>
  inner_join(
    actuals |> select(location, target_end_date, observation),
    by = c("location", "target_end_date")
  )

# 50% coverage: observation between q0.25 and q0.75
# 95% coverage: observation between q0.025 and q0.975
coverage_summary <- coverage_data |>
  pivot_wider(names_from = output_type_id, values_from = value, names_prefix = "q") |>
  mutate(
    in_50 = observation >= q0.25 & observation <= q0.75,
    in_95 = observation >= q0.025 & observation <= q0.975
  ) |>
  group_by(model) |>
  summarise(
    coverage_50 = mean(in_50, na.rm = TRUE),
    coverage_95 = mean(in_95, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

cat("\n  Coverage (Horizon 1):\n")
print(coverage_summary, n = 10)

# Combine all metrics
all_metrics <- wis_summary |>
  inner_join(coverage_summary, by = "model") |>
  # Also get MAE of median
  inner_join(
    coverage_data |>
      filter(output_type_id == 0.5) |>
      group_by(model) |>
      summarise(mae = mean(abs(value - observation), na.rm = TRUE), .groups = "drop"),
    by = "model"
  ) |>
  arrange(mean_wis)

cat("\n  Combined Metrics:\n")
print(all_metrics, n = 10)


# =========================================================================
# FIGURE 2: Correlation Heatmap (polished)
# =========================================================================

cat("\n--- Figure 2: Correlation heatmap ---\n")

# Recompute errors for correlation
errors_df <- forecasts_full |>
  filter(horizon == 1, output_type_id == 0.5) |>
  select(model, origin_date, location, target_end_date, forecast = value) |>
  inner_join(
    actuals |> select(location, target_end_date, actual = observation),
    by = c("location", "target_end_date")
  ) |>
  mutate(error = forecast - actual)

errors_wide <- errors_df |>
  select(model, origin_date, location, error) |>
  unite("row_id", origin_date, location, sep = "_") |>
  pivot_wider(names_from = model, values_from = error)

model_cols <- intersect(model_info$display_name, names(errors_wide))
cor_matrix <- cor(errors_wide |> select(all_of(model_cols)), use = "pairwise.complete.obs")

# Reorder by hierarchical clustering for visual grouping
hc <- hclust(as.dist(1 - cor_matrix))
ordered_names <- rownames(cor_matrix)[hc$order]

cor_long <- as.data.frame(as.table(cor_matrix)) |>
  rename(Model1 = Var1, Model2 = Var2, Correlation = Freq) |>
  mutate(
    Model1 = factor(Model1, levels = ordered_names),
    Model2 = factor(Model2, levels = rev(ordered_names))
  )

p2 <- ggplot(cor_long, aes(x = Model1, y = Model2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(
    label = sprintf("%.2f", Correlation),
    color = if_else(Correlation > 0.6, "white", "black")
  ), size = 3.5, fontface = "bold") +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1),
    name = "Error\nCorrelation"
  ) +
  labs(
    title = "Forecast Error Correlations",
    subtitle = "Low correlation = different mistakes = effective ensemble diversity",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(slides_dir, "02_correlation_heatmap.png"), p2,
       width = 8, height = 7, dpi = 200, bg = "white")
cat("  Saved 02_correlation_heatmap.png\n")


# =========================================================================
# FIGURE 3: WIS Comparison
# =========================================================================

cat("--- Figure 3: WIS comparison ---\n")

# Mark ensemble members
all_metrics <- all_metrics |>
  mutate(
    in_ensemble = model %in% c("NNETAR+BC", "SNAIVE+BC", "ETS+log", "ARIMA+BC", "HistWeek"),
    label = sprintf("%.3f", mean_wis)
  )

p3 <- ggplot(all_metrics, aes(x = reorder(model, mean_wis), y = mean_wis)) +
  geom_col(aes(fill = in_ensemble), alpha = 0.85, width = 0.7) +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5, fontface = "bold") +
  scale_fill_manual(
    values = c("TRUE" = "#2C7BB6", "FALSE" = "#BBBBBB"),
    labels = c("TRUE" = "Selected for ensemble", "FALSE" = "Dropped"),
    name = NULL
  ) +
  coord_flip() +
  labs(
    title = "Weighted Interval Score by Model (Horizon 1)",
    subtitle = "Lower = better calibrated probabilistic forecasts",
    x = NULL, y = "Mean WIS"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(slides_dir, "03_wis_comparison.png"), p3,
       width = 8, height = 5, dpi = 200, bg = "white")
cat("  Saved 03_wis_comparison.png\n")


# =========================================================================
# FIGURE 4: Coverage + WIS Combined Metrics Table
# =========================================================================

cat("--- Figure 4: Metrics table ---\n")

# Save as clean CSV for the presentation
metrics_for_slides <- all_metrics |>
  select(
    Model = model,
    `Mean WIS` = mean_wis,
    `MAE (Median)` = mae,
    `50% Coverage` = coverage_50,
    `95% Coverage` = coverage_95,
    `In Ensemble` = in_ensemble
  ) |>
  mutate(
    `Mean WIS` = round(`Mean WIS`, 3),
    `MAE (Median)` = round(`MAE (Median)`, 3),
    `50% Coverage` = sprintf("%.1f%%", `50% Coverage` * 100),
    `95% Coverage` = sprintf("%.1f%%", `95% Coverage` * 100)
  )

write_csv(metrics_for_slides, file.path(slides_dir, "04_metrics_table.csv"))
cat("  Saved 04_metrics_table.csv\n")
cat("\n")
print(metrics_for_slides)


# =========================================================================
# FIGURE 5: Regime-Specific Correlation (fixed bug)
# =========================================================================

cat("\n--- Figure 5: Regime-specific correlation ---\n")

# errors_df already has 'actual' from the join above
regime_df <- errors_df |>
  group_by(location) |>
  mutate(
    median_ili = median(actual, na.rm = TRUE),
    regime = if_else(actual > median_ili, "Peak (ILI > median)", "Trough (ILI <= median)")
  ) |>
  ungroup()

ensemble_models <- c("NNETAR+BC", "SNAIVE+BC", "ETS+log", "ARIMA+BC", "HistWeek")

regime_cors <- map_dfr(c("Peak (ILI > median)", "Trough (ILI <= median)"), function(r) {
  rwide <- regime_df |>
    filter(regime == r, model %in% ensemble_models) |>
    select(model, origin_date, location, error) |>
    unite("row_id", origin_date, location, sep = "_") |>
    pivot_wider(names_from = model, values_from = error)

  rcols <- intersect(ensemble_models, names(rwide))
  if (length(rcols) < 2) return(tibble())

  rcor <- cor(rwide |> select(all_of(rcols)), use = "pairwise.complete.obs")
  as.data.frame(as.table(rcor)) |>
    rename(Model1 = Var1, Model2 = Var2, Correlation = Freq) |>
    mutate(Regime = r)
})

p5 <- regime_cors |>
  filter(Model1 != Model2) |>
  ggplot(aes(x = Model1, y = Model2, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(
    label = sprintf("%.2f", Correlation),
    color = if_else(Correlation > 0.6, "white", "black")
  ), size = 3) +
  scale_color_identity() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0.5, limits = c(0, 1)
  ) +
  facet_wrap(~Regime) +
  labs(
    title = "Error Correlations by ILI Regime (Ensemble Models)",
    subtitle = "Models may diverge more during peaks - when diversity matters most",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave(file.path(slides_dir, "05_regime_correlation.png"), p5,
       width = 12, height = 5.5, dpi = 200, bg = "white")
cat("  Saved 05_regime_correlation.png\n")


# =========================================================================
# FIGURE 6: Fan Chart Comparison (3 models, same date)
# =========================================================================

cat("--- Figure 6: Fan chart comparison ---\n")

fan_date <- as.Date("2017-12-02")
fan_location <- "US National"

fan_data <- forecasts_full |>
  filter(
    origin_date == fan_date,
    location == fan_location,
    model %in% c("SNAIVE+BC", "ARIMA+BC", "NNETAR+BC")
  ) |>
  select(model, target_end_date, output_type_id, value) |>
  pivot_wider(names_from = output_type_id, values_from = value, names_prefix = "q")

fan_actuals <- actuals |>
  filter(
    location == fan_location,
    target_end_date >= fan_date - 56,
    target_end_date <= fan_date + 35
  )

p6 <- ggplot() +
  # 95% PI
  geom_ribbon(data = fan_data,
    aes(x = target_end_date, ymin = q0.025, ymax = q0.975, fill = model),
    alpha = 0.15
  ) +
  # 50% PI
  geom_ribbon(data = fan_data,
    aes(x = target_end_date, ymin = q0.25, ymax = q0.75, fill = model),
    alpha = 0.3
  ) +
  # Median
  geom_line(data = fan_data,
    aes(x = target_end_date, y = q0.5, color = model),
    linewidth = 1
  ) +
  # Actuals
  geom_line(data = fan_actuals,
    aes(x = target_end_date, y = observation),
    linewidth = 0.8, color = "black"
  ) +
  geom_point(data = fan_actuals,
    aes(x = target_end_date, y = observation),
    size = 1.5, color = "black"
  ) +
  geom_vline(xintercept = fan_date, linetype = "dashed", color = "gray50") +
  annotate("text", x = fan_date - 3, y = max(fan_actuals$observation) * 0.95,
           label = "Origin", hjust = 1, size = 3, color = "gray50") +
  scale_fill_manual(values = model_colors, name = "Model") +
  scale_color_manual(values = model_colors, name = "Model") +
  labs(
    title = sprintf("Fan Chart Comparison: %s (Origin: %s)", fan_location, fan_date),
    subtitle = "Shaded regions: 50% and 95% prediction intervals | Black: actual observed values",
    x = NULL, y = "Weighted ILI %"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40", size = 10),
    legend.position = "bottom"
  )

ggsave(file.path(slides_dir, "06_fan_chart_comparison.png"), p6,
       width = 10, height = 5.5, dpi = 200, bg = "white")
cat("  Saved 06_fan_chart_comparison.png\n")


# =========================================================================
# FIGURE 7: Error Time Series (when do models diverge?)
# =========================================================================

cat("--- Figure 7: Error time series ---\n")

error_ts <- errors_df |>
  filter(location == "US National", model %in% ensemble_models) |>
  select(model, origin_date, error, actual)

p7 <- ggplot(error_ts, aes(x = origin_date, y = error, color = model)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_line(alpha = 0.7, linewidth = 0.6) +
  scale_color_manual(values = model_colors, name = "Model") +
  labs(
    title = "Forecast Errors Over Time (US National, Horizon 1)",
    subtitle = "Models diverge most during peak transitions - where ensemble diversity pays off",
    x = NULL, y = "Error (Forecast - Actual)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40", size = 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(slides_dir, "07_error_timeseries.png"), p7,
       width = 10, height = 4.5, dpi = 200, bg = "white")
cat("  Saved 07_error_timeseries.png\n")


# =========================================================================
# FIGURE 8: Ensemble Architecture Diagram Data
# (Text-based for inclusion in slides as structured content)
# =========================================================================

cat("--- Figure 8: WIS by regime ---\n")

regime_wis <- wis_by_model |>
  inner_join(
    actuals |> select(location, target_end_date, observation),
    by = c("location"),
    relationship = "many-to-many"
  )

# Simpler approach: add actual to wis_by_model via origin_date + location
wis_with_actual <- wis_by_model |>
  inner_join(
    errors_df |> select(model, origin_date, location, actual) |> distinct(),
    by = c("model", "origin_date", "location")
  ) |>
  group_by(location) |>
  mutate(
    median_ili = median(actual, na.rm = TRUE),
    regime = if_else(actual > median_ili, "Peak", "Trough")
  ) |>
  ungroup()

regime_wis_summary <- wis_with_actual |>
  filter(model %in% ensemble_models) |>
  group_by(model, regime) |>
  summarise(mean_wis = mean(wis, na.rm = TRUE), .groups = "drop")

p8 <- ggplot(regime_wis_summary, aes(x = reorder(model, mean_wis), y = mean_wis, fill = regime)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = c("Peak" = "#D55E00", "Trough" = "#56B4E9"), name = "Regime") +
  coord_flip() +
  labs(
    title = "WIS by Regime (Ensemble Models)",
    subtitle = "All models struggle more during peaks - but by different amounts",
    x = NULL, y = "Mean WIS"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "gray40"),
    legend.position = "bottom"
  )

ggsave(file.path(slides_dir, "08_wis_by_regime.png"), p8,
       width = 8, height = 5, dpi = 200, bg = "white")
cat("  Saved 08_wis_by_regime.png\n")


# =========================================================================
# SAVE ALL METRICS AS RDS
# =========================================================================

presentation_data <- list(
  all_metrics     = all_metrics,
  wis_summary     = wis_summary,
  coverage_summary = coverage_summary,
  cor_matrix      = cor_matrix,
  ensemble_models = ensemble_models,
  regime_wis      = regime_wis_summary,
  timestamp       = Sys.time()
)

saveRDS(presentation_data, file.path(slides_dir, "presentation_data.rds"))

cat("\n=== ALL FIGURES SAVED ===\n")
cat(sprintf("  Output directory: %s\n", slides_dir))
cat("  Files:\n")
for (f in list.files(slides_dir)) {
  cat(sprintf("    %s\n", f))
}
