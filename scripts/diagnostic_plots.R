###############################################################################
# FILE: scripts/diagnostic_plots.R
#
# Visual QA for the SNAIVE + Box-Cox + Bootstrapped Residuals forecasts.
# Produces three types of plots:
#
#   (a) Fan charts for 5 representative origin dates — verify seasonal shape,
#       interval width, and actuals coverage.
#   (b) Sanity dashboard — median forecast (horizon 1) vs actuals over ALL
#       origin dates, plus 90% interval width over time.
#   (c) Lambda summary — bar chart of guerrero-estimated lambdas by location.
#
# Inputs:   139 CSV files in model-output/KReger-snaive_bc_bs/
# Outputs:  Plots displayed in the R graphics device.
#
# How to run:
#   source("scripts/KReger-snaive_bc_bs_forecast.R")  # loads data objects
#   source("scripts/diagnostic_plots.R")
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(patchwork)  # install.packages("patchwork") if needed


# =============================================================================
# PLOT (a): FAN CHARTS FOR 5 REPRESENTATIVE ORIGIN DATES
#
# Dates span different phases of the flu season:
#   - Early season (Nov): ILI just starting to rise
#   - Peak season (Jan): wide intervals expected
#   - Descending from peak (Mar): forecasts should decline
#   - Holiday period (Dec): potentially noisy reporting
#   - Late in dataset (Feb 2020): near end of available data
#
# What to look for:
#   - Red dots (actuals) fall within the fan most of the time
#   - Fan widens with increasing horizon (h1 narrower than h4)
#   - No negative y-axis values
#   - Upward fan spread is larger during peak season
# =============================================================================

cat("\n=== DIAGNOSTIC PLOT (a): Fan charts for 5 representative dates ===\n")

qa_dates <- as.Date(c(
  "2015-11-07",  # Early season
  "2017-01-14",  # Near a winter peak
  "2018-03-10",  # Descending from peak
  "2019-12-28",  # Holiday period
  "2020-02-15"   # Late in last season
))

qa_locations <- c("US National", "HHS Region 4", "HHS Region 10")

# Use index-based iteration to avoid R's for-loop Date-to-numeric coercion
for (i in seq_along(qa_dates)) {
  d <- qa_dates[i]
  date_str <- as.character(d)
  csv_file <- file.path(OUTPUT_DIR, paste0(date_str, "-", TEAM_MODEL, ".csv"))
  
  if (!file.exists(csv_file)) {
    cat("  Skipping", date_str, "— file not found\n")
    next
  }
  
  fc_df <- read_csv(csv_file, show_col_types = FALSE)
  
  # Fan chart for each location, combined with patchwork
  plots <- list()
  for (loc in qa_locations) {
    p <- plot_quantile_fan(fc_df, wili_tsibble, d, location = loc) +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
    plots <- c(plots, list(p))
  }
  
  combined <- wrap_plots(plots, ncol = 1) +
    plot_annotation(
      title = paste0("Fan charts for origin date: ", date_str),
      subtitle = "Blue = forecast fan (1-99%, 10-90%, 25-75%, median) | Red dots = actuals",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  print(combined)
  cat("  Plotted fan charts for", date_str, "\n")
}


# =============================================================================
# PLOT (b): SANITY DASHBOARD — MEDIAN vs ACTUALS OVER TIME
#
# For US National, reads all 139 CSVs and extracts horizon-1 data.
# Panel 1: Median forecast vs actual observed values
# Panel 2: Width of 90% prediction interval over time
#
# What to look for:
#   - Median tracks actuals (with ~1-year offset since SNAIVE)
#   - 90% PI width peaks during winter months
#   - No sudden jumps or flat spots suggesting data issues
# =============================================================================

cat("\n=== DIAGNOSTIC PLOT (b): Sanity dashboard — median vs actuals ===\n")

dashboard_loc <- "US National"

# Read all CSVs and extract horizon-1 data for US National
dashboard_data <- map_dfr(origin_dates, function(d) {
  csv_file <- file.path(OUTPUT_DIR, paste0(d, "-", TEAM_MODEL, ".csv"))
  if (!file.exists(csv_file)) return(NULL)
  
  read_csv(csv_file, show_col_types = FALSE) |>
    filter(location == dashboard_loc, horizon == 1) |>
    select(origin_date, target_end_date, output_type_id, value)
})

# Extract median (q0.5) and 90% interval bounds (q0.05, q0.95)
dashboard_summary <- dashboard_data |>
  filter(output_type_id %in% c(0.05, 0.5, 0.95)) |>
  pivot_wider(names_from = output_type_id, values_from = value,
              names_prefix = "q") |>
  mutate(interval_width_90 = q0.95 - q0.05)

actuals_for_dash <- wili_tsibble |>
  as_tibble() |>
  filter(location == dashboard_loc) |>
  select(target_end_date, observation)

dashboard_joined <- dashboard_summary |>
  left_join(actuals_for_dash, by = "target_end_date")

# Panel 1: Median forecast vs actuals
p_median <- ggplot(dashboard_joined, aes(x = target_end_date)) +
  geom_ribbon(aes(ymin = q0.05, ymax = q0.95),
              fill = "steelblue", alpha = 0.2) +
  geom_line(aes(y = observation), color = "black", linewidth = 0.4) +
  geom_line(aes(y = q0.5), color = "steelblue", linewidth = 0.6) +
  labs(
    title = paste0("Horizon-1 median forecast vs actuals: ", dashboard_loc),
    subtitle = "Blue = forecast median | Black = actual | Shaded = 90% PI",
    x = "Target end date", y = "wILI %"
  ) +
  theme_minimal()

# Panel 2: 90% interval width over time
p_width <- ggplot(dashboard_joined, aes(x = target_end_date, y = interval_width_90)) +
  geom_line(color = "darkorange", linewidth = 0.6) +
  labs(
    title = "90% prediction interval width (horizon 1)",
    subtitle = "Should be wider during peak flu season (Jan-Feb)",
    x = "Target end date", y = "Interval width (q0.95 - q0.05)"
  ) +
  theme_minimal()

dashboard_combined <- p_median / p_width +
  plot_annotation(
    title = paste0("Forecast Sanity Dashboard: ", dashboard_loc),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(dashboard_combined)
cat("  Dashboard plotted for", dashboard_loc, "\n")


# =============================================================================
# PLOT (c): LAMBDA SUMMARY — GUERRERO-ESTIMATED LAMBDAS BY LOCATION
#
# Verifies all lambdas are in a reasonable range:
#   λ ≈ 0 → log-like transform; λ ≈ 1 → no transform; λ < 0 → stronger
# =============================================================================

cat("\n=== DIAGNOSTIC PLOT (c): Lambda summary ===\n")

lambda_df <- tibble(
  location = names(lambdas),
  lambda   = unname(lambdas)
) |>
  arrange(lambda) |>
  mutate(location = factor(location, levels = location))

p_lambda <- ggplot(lambda_df, aes(x = location, y = lambda)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text(aes(label = sprintf("%.3f", lambda)),
            hjust = -0.1, size = 3) +
  coord_flip() +
  labs(
    title = "Guerrero-estimated Box-Cox lambda by location",
    subtitle = "λ ≈ 0 → log transform | λ ≈ 1 → no transform | λ < 0 → stronger stabilization",
    x = NULL, y = "Lambda (λ)"
  ) +
  scale_y_continuous(limits = c(min(lambda_df$lambda) - 0.05,
                                max(lambda_df$lambda) + 0.1)) +
  theme_minimal()

print(p_lambda)
cat("  Lambda summary plotted.\n")


# =============================================================================
# NUMERIC SUMMARY
# =============================================================================

cat("\n=== NUMERIC SUMMARY ===\n")
cat(sprintf("  Lambdas: min = %.3f, max = %.3f, median = %.3f\n",
            min(lambdas), max(lambdas), median(lambdas)))
cat(sprintf("  Dashboard location: %s\n", dashboard_loc))
cat(sprintf("  Median forecast range: %.2f to %.2f wILI%%\n",
            min(dashboard_joined$q0.5, na.rm = TRUE),
            max(dashboard_joined$q0.5, na.rm = TRUE)))
cat(sprintf("  Actual observation range: %.2f to %.2f wILI%%\n",
            min(dashboard_joined$observation, na.rm = TRUE),
            max(dashboard_joined$observation, na.rm = TRUE)))
cat(sprintf("  90%% PI width range: %.2f to %.2f\n",
            min(dashboard_joined$interval_width_90, na.rm = TRUE),
            max(dashboard_joined$interval_width_90, na.rm = TRUE)))

cat("\n=== Diagnostic plots complete ===\n")