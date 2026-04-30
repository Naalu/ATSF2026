###############################################################################
# FILE: scripts/twin_diagnostics_expanded.R  (v2)
#
# PURPOSE:
#   Phase 4B of FLUCAST_IMPLEMENTATION_PLAN.md. Apply twin-pair pruning to
#   the 20-model expanded candidate pool from Phase 4A using hierarchical
#   clustering of the point-error correlation matrix.
#
# v2 CHANGES vs v1:
#   1. Compute THREE linkage methods (single, complete, average) in
#      parallel and save all three shortlists. Report compares them.
#      "complete" is designated the official choice consumed by Phase 4C.
#   2. Survivor selection uses Phase 2's per-(horizon, regime) WIS cell
#      methodology — not aggregate WIS — so the algorithm can identify
#      regime-conditional dominance (e.g., a model that wins all peaking
#      cells survives even if its aggregate WIS is marginally worse).
#
# CLUSTERING METHODOLOGY:
#   - Distance matrix: d_ij = 1 - Pearson correlation of point errors
#   - Three linkage methods computed:
#       single   = min distance between cluster members (chains aggressively)
#       complete = max distance (most conservative; requires ALL pairs > 0.93)
#       average  = mean distance (UPGMA; balanced — was v1 default)
#   - Each method's tree cut at height 0.07 (correlation >= 0.93)
#
# SURVIVOR SELECTION (PER-CELL, Phase 2 methodology):
#   For each cluster of size > 1:
#     1. Build a 4 (horizons) x N_regimes grid of cells.
#     2. For each cell, identify the cluster member with lowest mean WIS.
#        A cell is "won" by that member only if its WIS is more than 2%
#        lower than the second-best member in the cell. Otherwise tie.
#     3. The cluster's survivor is the member who wins the most cells.
#     4. Tiebreaker: lowest mean dispersion (sharper intervals when
#        accuracy is comparable).
#
# WHY PER-CELL (vs aggregate WIS):
#   The earlier WIS-only-aggregate tiebreaker handed nnetar_log (WIS 0.2834)
#   the win over nnetar_bc_bs (WIS 0.2838) by a 0.0004 margin. Phase 2's
#   per-cell analysis showed nnetar_bc_bs wins all 3 peaking-regime cells
#   (where forecast quality matters most for public-health decisions).
#   Per-cell selection restores Phase 2's deliberate choice and is the
#   methodologically richer signal.
#
# CRITICAL: Validation-data only.
#   All scoring decisions use validation-period data (57 dates, Seasons
#   1-2). Test seasons remain unseen. Same hygiene as Phase 2.
#
# INPUTS:
#   - analysis/phase4/scores_long_expanded.csv
#   - analysis/phase4/per_model_summary_expanded.csv
#   - analysis/phase4/correlation_point_error_expanded.csv
#   - analysis/phase2/regimes.csv
#
# OUTPUTS in analysis/phase4/:
#   - cluster_decisions_<method>.csv     (one per linkage method)
#   - shortlist_expanded_<method>.csv    (one per linkage method)
#   - shortlist_expanded.csv             (= complete-linkage; Phase 4C input)
#   - cluster_decisions.csv              (= complete-linkage)
#   - dendrograms.pdf                    (3-panel comparison plot)
#   - twin_diagnostics_expanded.txt      (report-ready summary)
#
# HOW TO RUN:
#   setwd("<path-to-ATSF2026-repo>")
#   source("scripts/twin_diagnostics_expanded.R")
#
#   Runtime: ~1 minute.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md Phase 4B
#   - scripts/twin_diagnostics.R (Phase 2 — per-cell verdict template)
###############################################################################


# ============================================================================
# DEPENDENCIES
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
})


# ============================================================================
# CONFIGURATION
# ============================================================================
HUB_PATH       <- normalizePath(".", mustWork = FALSE)
PHASE4_DIR     <- file.path(HUB_PATH, "analysis", "phase4")
PHASE2_DIR     <- file.path(HUB_PATH, "analysis", "phase2")

SCORES_LONG    <- file.path(PHASE4_DIR, "scores_long_expanded.csv")
PER_MODEL      <- file.path(PHASE4_DIR, "per_model_summary_expanded.csv")
COR_MATRIX     <- file.path(PHASE4_DIR, "correlation_point_error_expanded.csv")
REGIMES_CSV    <- file.path(PHASE2_DIR, "regimes.csv")

# Twin threshold and corresponding distance cut height.
TWIN_THRESHOLD <- 0.93
CUT_HEIGHT     <- 1 - TWIN_THRESHOLD

# Linkage methods to compare.
LINKAGE_METHODS  <- c("single", "complete", "average")
OFFICIAL_LINKAGE <- "complete"

# Per-cell margin (matches Phase 2's twin_diagnostics.R).
WIS_MARGIN <- 0.02


# ============================================================================
# LOAD INPUTS
# ============================================================================
cat("\n--- Loading Phase 4A scores and correlation matrix ---\n")

scores_long <- readr::read_csv(SCORES_LONG, show_col_types = FALSE) |>
  dplyr::mutate(
    origin_date     = as.Date(origin_date),
    target_end_date = as.Date(target_end_date)
  )

per_model <- readr::read_csv(PER_MODEL, show_col_types = FALSE)

cor_df <- readr::read_csv(COR_MATRIX, show_col_types = FALSE)
labels_vec <- cor_df$model_label
cor_matrix <- as.matrix(cor_df[, -1])
rownames(cor_matrix) <- labels_vec
colnames(cor_matrix) <- labels_vec

regimes <- readr::read_csv(REGIMES_CSV, show_col_types = FALSE) |>
  dplyr::mutate(origin_date = as.Date(origin_date))

cat(sprintf("  scores_long: %d rows\n", nrow(scores_long)))
cat(sprintf("  per_model:   %d rows\n", nrow(per_model)))
cat(sprintf("  correlation: %dx%d matrix\n", nrow(cor_matrix), ncol(cor_matrix)))

# Pre-join regimes to scores for the survivor function.
scores_with_regime <- scores_long |>
  dplyr::inner_join(regimes, by = c("origin_date", "location"))


# ============================================================================
# SURVIVOR SELECTION FUNCTION
#
# For a cluster (vector of `model` directory names), pick the surviving
# member using per-(horizon, regime) WIS cells.
# ============================================================================

#' Pick a cluster's survivor using per-cell WIS comparison.
#'
#' @param cluster_models Character vector of `model` directory names.
#'   Must have length >= 2. Singletons handled outside.
#' @param scoring_table Long table with: model, horizon, regime, wis,
#'   wis_dispersion (regimes already joined).
#' @return Character: the surviving model's directory name.
pick_cluster_survivor <- function(cluster_models, scoring_table) {

  if (length(cluster_models) < 2L) {
    stop(sprintf("pick_cluster_survivor needs >= 2 members, got %d",
                 length(cluster_models)))
  }

  # Filter to cluster members and compute per-cell mean WIS.
  cluster_scores <- scoring_table |>
    dplyr::filter(model %in% cluster_models) |>
    dplyr::group_by(model, horizon, regime) |>
    dplyr::summarise(
      cell_wis  = mean(wis, na.rm = TRUE),
      cell_disp = mean(wis_dispersion, na.rm = TRUE),
      .groups   = "drop"
    )

  # For each (horizon, regime) cell: rank members by WIS, identify the
  # winner if margin to second-best exceeds WIS_MARGIN.
  cell_winners <- cluster_scores |>
    dplyr::group_by(horizon, regime) |>
    dplyr::arrange(cell_wis, .by_group = TRUE) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::summarise(
      best_model  = model[rank == 1L][1],
      best_wis    = cell_wis[rank == 1L][1],
      second_wis  = cell_wis[rank == 2L][1],
      rel_diff    = (second_wis - best_wis) / pmax(best_wis, second_wis),
      cell_winner = ifelse(rel_diff > WIS_MARGIN, best_model, "tie"),
      .groups     = "drop"
    )

  # Count cells won by each member, including 0-win members.
  win_counts <- cell_winners |>
    dplyr::filter(cell_winner != "tie") |>
    dplyr::count(cell_winner, name = "n_wins") |>
    dplyr::rename(model = cell_winner)

  all_members <- tibble::tibble(model = cluster_models)
  win_counts <- all_members |>
    dplyr::left_join(win_counts, by = "model") |>
    dplyr::mutate(n_wins = dplyr::coalesce(n_wins, 0L))

  # Tiebreaker: aggregate dispersion across all cluster cells.
  disp_per_member <- cluster_scores |>
    dplyr::group_by(model) |>
    dplyr::summarise(mean_disp = mean(cell_disp, na.rm = TRUE),
                     .groups = "drop")

  win_counts <- win_counts |>
    dplyr::left_join(disp_per_member, by = "model") |>
    dplyr::arrange(dplyr::desc(n_wins), mean_disp)

  win_counts$model[1]
}


# ============================================================================
# LINKAGE METHOD RUNNER
# ============================================================================

#' Run hierarchical clustering with one linkage method, then survivor pick.
run_linkage <- function(method) {
  cat(sprintf("\n--- Linkage: %s ---\n", method))

  dist_obj <- as.dist(1 - cor_matrix)
  hc <- stats::hclust(dist_obj, method = method)
  cluster_assignments <- stats::cutree(hc, h = CUT_HEIGHT)

  clusters_df <- tibble::tibble(
    model_label = names(cluster_assignments),
    cluster     = cluster_assignments
  )

  # Map labels back to directory `model` names.
  label_to_model <- per_model |> dplyr::select(model_label, model)
  clusters_with_model <- clusters_df |>
    dplyr::inner_join(label_to_model, by = "model_label") |>
    dplyr::inner_join(per_model |> dplyr::select(model, wis, wis_disp),
                      by = "model")

  cluster_ids <- sort(unique(clusters_with_model$cluster))

  decisions_list <- list()
  for (cl in cluster_ids) {
    members <- clusters_with_model |> dplyr::filter(cluster == cl)
    cluster_size <- nrow(members)

    if (cluster_size == 1L) {
      survivor <- members$model[1]
      decision_method <- "singleton"
    } else {
      survivor <- pick_cluster_survivor(members$model, scores_with_regime)
      decision_method <- "per_cell_wis"
    }

    decisions_list[[as.character(cl)]] <- members |>
      dplyr::mutate(
        is_survivor      = (model == survivor),
        cluster_size     = cluster_size,
        decision_method  = decision_method,
        cluster_members  = paste(sort(model_label), collapse = ", ")
      )
  }

  cluster_decisions <- dplyr::bind_rows(decisions_list) |>
    dplyr::arrange(cluster, dplyr::desc(is_survivor), wis)

  shortlist <- cluster_decisions |>
    dplyr::filter(is_survivor) |>
    dplyr::transmute(cluster, model_label, model, wis, wis_disp,
                     cluster_size, cluster_members) |>
    dplyr::arrange(wis)

  n_clusters   <- length(cluster_ids)
  n_singletons <- sum(table(clusters_df$cluster) == 1)
  n_with_twins <- n_clusters - n_singletons

  cat(sprintf("  %d clusters (%d singletons + %d with twins) -> %d survivors\n",
              n_clusters, n_singletons, n_with_twins, nrow(shortlist)))

  list(
    cluster_decisions = cluster_decisions,
    shortlist         = shortlist,
    hclust_obj        = hc
  )
}


# ============================================================================
# RUN ALL THREE LINKAGE METHODS
# ============================================================================
cat("\n=== Running all three linkage methods ===\n")

results <- list()
for (m in LINKAGE_METHODS) {
  results[[m]] <- run_linkage(m)

  readr::write_csv(
    results[[m]]$cluster_decisions,
    file.path(PHASE4_DIR, sprintf("cluster_decisions_%s.csv", m))
  )
  readr::write_csv(
    results[[m]]$shortlist,
    file.path(PHASE4_DIR, sprintf("shortlist_expanded_%s.csv", m))
  )
}


# ============================================================================
# OFFICIAL DESIGNATION
# ============================================================================
cat(sprintf("\n--- Designating '%s' as official ---\n", OFFICIAL_LINKAGE))

official_shortlist <- results[[OFFICIAL_LINKAGE]]$shortlist
official_decisions <- results[[OFFICIAL_LINKAGE]]$cluster_decisions

readr::write_csv(official_shortlist,
                 file.path(PHASE4_DIR, "shortlist_expanded.csv"))
readr::write_csv(official_decisions,
                 file.path(PHASE4_DIR, "cluster_decisions.csv"))

cat(sprintf("  shortlist_expanded.csv = %s linkage's shortlist\n",
            OFFICIAL_LINKAGE))


# ============================================================================
# DISPLAY OFFICIAL DECISIONS
# ============================================================================
cat(sprintf("\n=== Official shortlist (%s linkage, per-cell survivor) ===\n",
            OFFICIAL_LINKAGE))

cat("\n  Cluster decisions:\n\n")
for (cl in sort(unique(official_decisions$cluster))) {
  members <- official_decisions |> dplyr::filter(cluster == cl)
  size    <- nrow(members)
  surv_idx <- which(members$is_survivor)
  surv    <- members$model_label[surv_idx]
  surv_wis <- members$wis[surv_idx]

  if (size == 1L) {
    cat(sprintf("  Cluster %2d  singleton:        KEEP %-15s (WIS = %.4f)\n",
                cl, surv, surv_wis))
  } else {
    cat(sprintf("  Cluster %2d  size %d, per-cell:  KEEP %-15s (WIS = %.4f)\n",
                cl, size, surv, surv_wis))
    losers <- members[!members$is_survivor, ]
    for (i in seq_len(nrow(losers))) {
      cat(sprintf("                                DROP %-15s (WIS = %.4f, +%.4f)\n",
                  losers$model_label[i], losers$wis[i],
                  losers$wis[i] - surv_wis))
    }
  }
}

cat(sprintf("\n  Final shortlist (%d models, sorted by WIS):\n",
            nrow(official_shortlist)))
print(
  official_shortlist |>
    dplyr::transmute(
      model_label,
      wis = round(wis, 4),
      wis_disp = round(wis_disp, 4),
      cluster_size,
      "kr/up" = ifelse(grepl("^KR_", model_label), "kr", "up")
    ),
  n = Inf
)


# ============================================================================
# CROSS-METHOD COMPARISON
# ============================================================================
cat("\n=== Cross-method shortlist comparison ===\n")

cross_method <- per_model |>
  dplyr::select(model_label, wis) |>
  dplyr::mutate(
    in_single   = model_label %in% results[["single"]]$shortlist$model_label,
    in_complete = model_label %in% results[["complete"]]$shortlist$model_label,
    in_average  = model_label %in% results[["average"]]$shortlist$model_label
  ) |>
  dplyr::arrange(wis)

cat("\n  Per-method inclusion (1 = in shortlist):\n\n")
cross_method |>
  dplyr::transmute(
    model_label, wis = round(wis, 4),
    single = as.integer(in_single),
    complete = as.integer(in_complete),
    average = as.integer(in_average)
  ) |>
  print(n = Inf)

cat(sprintf("\n  Single   linkage: %d survivors\n",
            nrow(results[["single"]]$shortlist)))
cat(sprintf("  Complete linkage: %d survivors  <- OFFICIAL\n",
            nrow(results[["complete"]]$shortlist)))
cat(sprintf("  Average  linkage: %d survivors\n",
            nrow(results[["average"]]$shortlist)))


# ============================================================================
# COMPARISON DENDROGRAM
# ============================================================================
cat("\n--- Saving 3-panel dendrogram comparison ---\n")
dendrogram_path <- file.path(PHASE4_DIR, "dendrograms.pdf")

grDevices::pdf(dendrogram_path, width = 14, height = 6)
graphics::par(mfrow = c(1, 3))
for (m in LINKAGE_METHODS) {
  is_official <- (m == OFFICIAL_LINKAGE)
  marker <- if (is_official) " (OFFICIAL)" else ""
  plot(
    results[[m]]$hclust_obj,
    main = sprintf("Linkage: %s%s", m, marker),
    sub  = sprintf("Cut h = %.2f (corr >= %.2f)", CUT_HEIGHT, TWIN_THRESHOLD),
    xlab = "",
    cex  = 0.65
  )
  graphics::abline(h = CUT_HEIGHT, col = "red", lty = 2, lwd = 1.5)
}
graphics::par(mfrow = c(1, 1))
grDevices::dev.off()
cat(sprintf("  Wrote: %s\n", dendrogram_path))


# ============================================================================
# REPORT-READY SUMMARY
# ============================================================================
summary_path <- file.path(PHASE4_DIR, "twin_diagnostics_expanded.txt")
sink(summary_path)
cat("Phase 4B v2: twin-pair pruning of expanded candidate pool\n")
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat(sprintf("\nMethodology:\n"))
cat(sprintf("  Distance:          1 - Pearson cor of point-forecast errors\n"))
cat(sprintf("  Twin threshold:    %.2f (cut at distance = %.2f)\n",
            TWIN_THRESHOLD, CUT_HEIGHT))
cat(sprintf("  Survivor:          per-(horizon, regime) cell winner count\n"))
cat(sprintf("  Tiebreak:          mean dispersion (Phase 2 methodology)\n"))
cat(sprintf("  Cell win margin:   %.0f%% relative WIS\n", 100 * WIS_MARGIN))

cat(sprintf("\nLinkage method comparison:\n"))
n_kr_in <- sum(grepl("^KR_", per_model$model_label))
n_up_in <- nrow(per_model) - n_kr_in
for (m in LINKAGE_METHODS) {
  marker <- if (m == OFFICIAL_LINKAGE) " <- OFFICIAL" else ""
  sl <- results[[m]]$shortlist
  n_kr <- sum(grepl("^KR_", sl$model_label))
  n_up <- nrow(sl) - n_kr
  cat(sprintf("  %-9s: %2d total (KReger %d/%d, Upstream %d/%d)%s\n",
              m, nrow(sl), n_kr, n_kr_in, n_up, n_up_in, marker))
}

cat(sprintf("\nOFFICIAL shortlist (%s linkage):\n", OFFICIAL_LINKAGE))
print(official_shortlist |>
        dplyr::transmute(model_label, model, wis = round(wis, 4),
                         cluster_size))
sink()
cat(sprintf("\n  Wrote: %s\n", summary_path))


# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n========================================\n")
cat("  Phase 4B v2 complete\n")
cat("========================================\n")
cat(sprintf("  Three linkage methods evaluated:\n"))
for (m in LINKAGE_METHODS) {
  marker <- if (m == OFFICIAL_LINKAGE) " <- OFFICIAL" else ""
  cat(sprintf("    %-9s -> %2d survivors%s\n",
              m, nrow(results[[m]]$shortlist), marker))
}
cat(sprintf("\n  Outputs in %s:\n", PHASE4_DIR))
cat("    - cluster_decisions_<method>.csv  (3 files)\n")
cat("    - shortlist_expanded_<method>.csv (3 files)\n")
cat("    - shortlist_expanded.csv          (= complete; for Phase 4C)\n")
cat("    - cluster_decisions.csv           (= complete)\n")
cat("    - dendrograms.pdf                 (3-panel comparison)\n")
cat("    - twin_diagnostics_expanded.txt   (report summary)\n")
cat("========================================\n")
