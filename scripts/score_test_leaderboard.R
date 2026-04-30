# =============================================================================
# scripts/score_test_leaderboard.R
#
# Phase 5D-bonus — Cross-class test-set leaderboard
#
# Scores every model in model-output/ on the test window using the same
# methodology as the ATSF2026 dashboard, producing a leaderboard comparable
# to https://sjfox.github.io/ATSF2026-dashboard/eval.html.
#
# To match the dashboard's reported counts as closely as possible, the
# scoring uses each model's actual coverage of the test window — no
# pre-filtering to a common subset. This means models with partial test
# coverage are scored only on the cells they cover (matching the dashboard's
# per-model N column).
#
# Outputs (saved to analysis/phase4/):
#   - test_leaderboard.csv:        full leaderboard with mean WIS, MAE, coverage
#   - test_leaderboard_report.md:  readable markdown summary
#
# Run from project root:
#   Rscript scripts/score_test_leaderboard.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scoringutils)
})

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------

# Test window — first origin date after validation ends
TEST_START <- as.Date("2017-10-28")
TEST_END   <- as.Date("2020-08-29")  # extend through summer to capture all
                                      # test-period submissions; truth filter
                                      # will drop forecasts past truth coverage

OUT_DIR  <- "analysis/phase4"
OUT_CSV  <- file.path(OUT_DIR, "test_leaderboard.csv")
OUT_MD   <- file.path(OUT_DIR, "test_leaderboard_report.md")

# -----------------------------------------------------------------------------
# 1. Load truth
# -----------------------------------------------------------------------------

cat("=== Loading truth ===\n")

truth <- readr::read_csv("target-data/oracle-output.csv",
                         show_col_types = FALSE) |>
  dplyr::filter(target == "ili perc") |>
  dplyr::transmute(
    location,
    target_end_date = as.Date(target_end_date),
    observed        = oracle_value
  ) |>
  dplyr::distinct(location, target_end_date, .keep_all = TRUE)

cat(sprintf("  %d truth rows; date range %s to %s\n",
            nrow(truth),
            min(truth$target_end_date),
            max(truth$target_end_date)))

# -----------------------------------------------------------------------------
# 2. Discover models
# -----------------------------------------------------------------------------

# Each subdirectory of model-output/ is a (team-model) submission.
# Skip non-directories like README.md.
all_dirs <- list.files("model-output", full.names = FALSE)
model_names <- all_dirs[file.info(file.path("model-output", all_dirs))$isdir]
model_names <- model_names[!is.na(model_names)]

cat(sprintf("\n=== Found %d models ===\n", length(model_names)))
cat(paste(" -", model_names, collapse = "\n"), "\n")

# -----------------------------------------------------------------------------
# 3. Function: load + score one model
# -----------------------------------------------------------------------------

#' Load a model's test-period forecasts and score against truth.
#'
#' @param model_name Character. Name of the team-model directory under
#'   model-output/ (e.g., "delphi-epicast", "KReger-bma_ensemble").
#' @return One-row tibble with model name, mean WIS, mean MAE, coverage,
#'   and N cells. Returns NULL if the model has no scorable cells.
score_model <- function(model_name) {
  dir <- file.path("model-output", model_name)
  csv_files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Filter to test-period origin dates by filename
  origin_dates <- as.Date(stringr::str_extract(basename(csv_files),
                                                "^\\d{4}-\\d{2}-\\d{2}"))
  test_files <- csv_files[origin_dates >= TEST_START & origin_dates <= TEST_END]
  
  if (length(test_files) == 0L) {
    cat(sprintf("  %s: no test-period CSVs\n", model_name))
    return(NULL)
  }
  
  # Read and combine. Use col_types = cols(.default = "c") then re-cast
  # because some models may have NA-typed columns that confuse purrr::map_dfr.
  forecasts <- purrr::map_dfr(test_files, function(f) {
    od <- as.Date(stringr::str_extract(basename(f), "^\\d{4}-\\d{2}-\\d{2}"))
    suppressWarnings(
      readr::read_csv(f, show_col_types = FALSE,
                      col_types = readr::cols(.default = readr::col_character())) |>
        dplyr::mutate(origin_date = od)
    )
  }) |>
    dplyr::filter(output_type == "quantile") |>
    dplyr::transmute(
      origin_date,
      location        = as.character(location),
      horizon         = as.integer(horizon),
      target_end_date = as.Date(target_end_date),
      quantile_level  = as.numeric(output_type_id),
      predicted       = as.numeric(value)
    )
  
  # Join with truth, drop rows past truth coverage
  combined <- forecasts |>
    dplyr::left_join(truth, by = c("location", "target_end_date")) |>
    dplyr::filter(!is.na(observed), !is.na(predicted))
  
  if (nrow(combined) == 0L) {
    cat(sprintf("  %s: no rows after truth join\n", model_name))
    return(NULL)
  }
  
  # Score
  forecast_obj <- tryCatch(
    combined |>
      scoringutils::as_forecast_quantile(
        forecast_unit  = c("origin_date", "location", "horizon",
                           "target_end_date"),
        observed       = "observed",
        predicted      = "predicted",
        quantile_level = "quantile_level"
      ),
    error = function(e) {
      cat(sprintf("  %s: as_forecast_quantile failed: %s\n",
                  model_name, e$message))
      NULL
    }
  )
  
  if (is.null(forecast_obj)) return(NULL)
  
  scores <- scoringutils::score(forecast_obj)
  
  # Derive headline metrics. scoringutils 2.x exposes interval coverage via
  # add_coverage(); but for simplicity we compute coverage from quantile_level
  # crossings ourselves.
  # mean WIS, MAE = abs(median - observed)
  cov_50 <- combined |>
    dplyr::filter(quantile_level %in% c(0.25, 0.75)) |>
    tidyr::pivot_wider(names_from = quantile_level,
                       names_prefix = "q",
                       values_from = predicted) |>
    dplyr::filter(!is.na(q0.25), !is.na(q0.75)) |>
    dplyr::summarise(
      cov_50 = mean(observed >= q0.25 & observed <= q0.75) * 100
    ) |>
    dplyr::pull(cov_50)
  
  cov_95 <- combined |>
    dplyr::filter(quantile_level %in% c(0.025, 0.975)) |>
    tidyr::pivot_wider(names_from = quantile_level,
                       names_prefix = "q",
                       values_from = predicted) |>
    dplyr::filter(!is.na(q0.025), !is.na(q0.975)) |>
    dplyr::summarise(
      cov_95 = mean(observed >= q0.025 & observed <= q0.975) * 100
    ) |>
    dplyr::pull(cov_95)
  
  mae <- combined |>
    dplyr::filter(quantile_level == 0.5) |>
    dplyr::summarise(mae = mean(abs(predicted - observed))) |>
    dplyr::pull(mae)
  
  tibble::tibble(
    model     = model_name,
    mean_wis  = mean(scores$wis),
    mae       = mae,
    cov_50    = cov_50,
    cov_95    = cov_95,
    n_cells   = nrow(scores)
  )
}

# -----------------------------------------------------------------------------
# 4. Score all models
# -----------------------------------------------------------------------------

cat("\n=== Scoring ===\n")

leaderboard <- purrr::map_dfr(model_names, function(m) {
  cat(sprintf("\n%s\n", m))
  score_model(m)
})

leaderboard <- leaderboard |>
  dplyr::arrange(mean_wis) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::select(rank, model, mean_wis, mae, cov_50, cov_95, n_cells)

# -----------------------------------------------------------------------------
# 5. Display
# -----------------------------------------------------------------------------

cat("\n=== Leaderboard (test season) ===\n\n")

leaderboard |>
  dplyr::mutate(
    mean_wis = round(mean_wis, 3),
    mae      = round(mae,      3),
    cov_50   = round(cov_50,   1),
    cov_95   = round(cov_95,   1)
  ) |>
  print(n = Inf)

readr::write_csv(leaderboard, OUT_CSV)
cat(sprintf("\nSaved: %s\n", OUT_CSV))

# -----------------------------------------------------------------------------
# 6. Markdown report
# -----------------------------------------------------------------------------

md_lines <- c(
  "# Test-Set Leaderboard",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  sprintf("Test window: origin dates from %s onward", TEST_START),
  "",
  "Comparable to https://sjfox.github.io/ATSF2026-dashboard/eval.html",
  "Methodology: scoringutils 2.x; per-model coverage of the test window.",
  "",
  "| Rank | Model | Mean WIS | MAE | 50% Cov. | 95% Cov. | N |",
  "|---:|---|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(leaderboard))) {
  row <- leaderboard[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %d | %s | %.3f | %.3f | %.1f | %.1f | %d |",
    row$rank, row$model, row$mean_wis, row$mae,
    row$cov_50, row$cov_95, row$n_cells
  ))
}

# Highlight KReger-bma_ensemble's rank if present
if ("KReger-bma_ensemble" %in% leaderboard$model) {
  rk <- leaderboard$rank[leaderboard$model == "KReger-bma_ensemble"]
  md_lines <- c(md_lines, "",
                sprintf("**KReger-bma_ensemble rank: %d of %d**",
                        rk, nrow(leaderboard)))
}

writeLines(md_lines, OUT_MD)
cat(sprintf("Saved: %s\n", OUT_MD))

cat("\n=== Done ===\n")
