###############################################################################
# FILE: scripts/regime_helpers.R
#
# PURPOSE:
#   Classify (origin_date, location) cells into one of four flu-season
#   regimes for use by the regime-aware BMA ensemble in Phase 3 and by the
#   twin-pair diagnostics in Phase 2.
#
#   Regimes:
#     "growing"    — wILI is increasing; epidemic ramping up
#     "peaking"    — wILI is at high level; near or at season peak
#     "declining"  — wILI is decreasing; epidemic winding down
#     "off_season" — wILI is low and flat; summer/early-fall lull
#
# DESIGN — 4-step priority cascade (per location):
#   1. Slope test. Fit a linear regression to the last 4 weeks of wILI
#      ending at as_of_date. Standardize the slope by dividing by the
#      rolling SD of those weeks. If |z-slope| > slope_z_threshold,
#      classify as "growing" (positive) or "declining" (negative).
#      Otherwise fall through.
#
#   2. Level test. Compare the most-recent 4-week mean wILI against this
#      location's historical 90th and 10th percentiles (computed from
#      pre-as_of_date data only). Above 90th -> "peaking", below 10th ->
#      "off_season". Otherwise fall through.
#
#   3. National tiebreaker. If steps 1 and 2 both fell through, run the
#      same slope+level checks on US National wILI. If those classify it,
#      apply that label to the local location. National usually leads
#      regional, so this catches cases where local signal is ambiguous
#      but a regime shift is happening at the national level.
#
#   4. Calendar fallback. Use the location's empirical week-of-year wILI
#      profile (computed from pre-validation data) to assign one of the
#      four regimes based on what the average season looks like at this
#      week-of-year for this location.
#
# DATA-DRIVEN CALIBRATION:
#   Both the slope threshold (step 1) and the calendar boundaries (step 4)
#   are derived from training data, NOT hardcoded. This means the classifier
#   is calibrated to the actual signal-to-noise of the wILI series and to
#   the actual seasonal pattern, rather than to a reviewer-questionable
#   "0.05" or "Oct-Dec".
#
# CRITICAL: As-of constraint.
#   Every threshold and boundary is computed from data strictly before
#   as_of_date. No future information leaks into the regime label assigned
#   at any origin date. This is the "honest as-of forecasting" requirement
#   from the user's design choices.
#
# PUBLIC API:
#   compute_regime_calibration(wili_data, validation_start_date)
#     -> list with $slope_z_threshold and $location_calendar (lookup table)
#
#   classify_regime(wili_data, as_of_date, location, calibration)
#     -> character: one of c("growing", "peaking", "declining", "off_season")
#
#   precompute_regimes(wili_data, origin_dates, locations, calibration)
#     -> tibble with columns: origin_date, location, regime
#
# INPUTS:
#   wili_data: tibble/tsibble with columns location (chr),
#     target_end_date (Date), observation (dbl). The hub's
#     target-data/time-series.csv loaded into this format.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md §1.6 (regime classifier was designed
#     but not previously coded)
#   - User memory: 4-step priority cascade with smoothed slope + location-
#     specific thresholds + level percentiles + seasonal default
###############################################################################

# We don't library() anything at the top — this file is sourced by other
# scripts that already have dplyr/tibble loaded. Calling functions are
# qualified with `dplyr::` etc. to be explicit about provenance.


# ============================================================================
# PRIVATE HELPERS
# ============================================================================

#' Compute the 4-week mean wILI ending at as_of_date for one location.
#'
#' Returns NA if fewer than 4 observations are available in the window —
#' the caller should treat NA as "fall through to the next cascade step".
#'
#' @param wili_data Long tibble with location, target_end_date, observation.
#' @param as_of_date Date. The "as-of" cutoff.
#' @param loc Character. Location to filter to.
#' @return Numeric scalar, the 4-week mean, or NA.
.recent_mean <- function(wili_data, as_of_date, loc) {
  recent <- wili_data |>
    dplyr::filter(location == loc,
                  target_end_date <= as_of_date,
                  target_end_date > as_of_date - 28L) |>
    dplyr::pull(observation)
  if (length(recent) < 4L) return(NA_real_)
  mean(recent, na.rm = TRUE)
}


#' Compute the standardized slope of recent wILI for one location.
#'
#' Fits y ~ t for the last 4 weeks of data ending at as_of_date, then
#' divides the slope by the rolling SD of the same window. The result is
#' interpretable as "how many standard deviations of typical week-to-week
#' variation does this slope represent per week?".
#'
#' Returns a list with $z (standardized slope) and $sign (raw sign), or
#' NULL if data is insufficient.
.standardized_slope <- function(wili_data, as_of_date, loc) {
  window <- wili_data |>
    dplyr::filter(location == loc,
                  target_end_date <= as_of_date,
                  target_end_date > as_of_date - 28L) |>
    dplyr::arrange(target_end_date)

  if (nrow(window) < 4L) return(NULL)

  # Time variable in weeks (0, 1, 2, 3) — slope is then "wILI per week".
  t <- as.numeric(window$target_end_date - min(window$target_end_date)) / 7
  y <- window$observation

  # Defensive: if all values are identical, slope is 0 and SD is 0 — would
  # produce NaN. Treat as "no signal" and return z = 0.
  if (sd(y, na.rm = TRUE) < 1e-9) {
    return(list(z = 0, sign = 0))
  }

  # OLS slope. lm() is overkill for 4 points but reads more clearly than
  # cov(t,y)/var(t) and produces identical results.
  fit <- stats::lm(y ~ t)
  slope <- coef(fit)[["t"]]

  # Standardize: z = slope / sd(y). Interpretable as "this slope's
  # magnitude in units of typical wILI variability per week".
  z <- slope / sd(y, na.rm = TRUE)

  list(z = z, sign = sign(slope))
}


#' Look up location-specific level thresholds from the calibration object.
#'
#' Returns the 10th and 90th percentiles of historical wILI for this
#' location, computed from data strictly before validation_start_date.
.level_thresholds <- function(calibration, loc) {
  th <- calibration$level_thresholds
  matched <- th[th$location == loc, ]
  if (nrow(matched) == 0L) {
    stop(sprintf("No level thresholds calibrated for location '%s'", loc))
  }
  list(p10 = matched$p10[1], p90 = matched$p90[1])
}


#' Look up the calendar-fallback regime for one location and one date.
.calendar_regime <- function(calibration, loc, as_of_date) {
  cal <- calibration$location_calendar
  woy <- .week_of_year(as_of_date)
  matched <- cal[cal$location == loc & cal$week_of_year == woy, ]
  if (nrow(matched) == 0L) {
    # Should not happen — every (loc, woy) is calibrated. Defensive default.
    return("off_season")
  }
  matched$regime[1]
}


#' Compute MMWR-ish week of year (1-53) from a Date.
#'
#' We use base R's %V format string (ISO 8601 week), which is close enough
#' to MMWR for our calendar-fallback purposes — the boundaries we derive
#' are smooth enough that ISO vs MMWR week numbering doesn't matter at
#' the regime-label level.
.week_of_year <- function(d) {
  as.integer(format(d, "%V"))
}


# ============================================================================
# PUBLIC: CALIBRATION
# ============================================================================

#' Compute data-driven thresholds from pre-validation wILI history.
#'
#' Calibration consists of:
#'   - slope_z_threshold: derived from the empirical distribution of
#'     standardized slopes in the training data. We use the 67th percentile
#'     of |z|. Slopes above this are "clearly trending" — the top third of
#'     historical slope magnitudes.
#'   - level_thresholds: for each location, the 10th and 90th percentiles
#'     of historical weekly wILI. Used by step 2 of the cascade.
#'   - location_calendar: for each (location, week_of_year), the empirical
#'     average wILI percentile and a derived regime label.
#'
#' All quantities are computed from data with target_end_date strictly
#' less than validation_start_date — no validation/test contamination.
#'
#' @param wili_data Long tibble: location, target_end_date, observation.
#' @param validation_start_date Date. First validation origin date
#'   (typically 2015-10-24 for this project).
#' @return List with three named elements: slope_z_threshold (scalar),
#'   level_thresholds (tibble), location_calendar (tibble).
compute_regime_calibration <- function(wili_data, validation_start_date) {

  # ----- Subset to pre-validation history -----
  train <- wili_data |>
    dplyr::filter(target_end_date < validation_start_date)

  if (nrow(train) == 0L) {
    stop("No pre-validation training data available for regime calibration.")
  }

  # ----- Slope threshold -----
  # Compute the standardized slope at every Saturday in the training data
  # for which we have at least 4 weeks of prior history. This yields a
  # cloud of |z| values whose 67th percentile defines "clearly trending".
  saturdays <- train |>
    dplyr::filter(format(target_end_date, "%w") == "6") |>
    dplyr::distinct(target_end_date) |>
    dplyr::pull(target_end_date) |>
    sort()

  # Restrict to dates with enough lead-up history (need 4 prior weeks).
  saturdays <- saturdays[saturdays >= min(train$target_end_date) + 28L]

  locations <- unique(train$location)

  # Brute-force loop. ~600 dates × 11 locations = ~6,600 standardizations.
  # Each takes <1 ms; total runtime <10 sec. Optimization not needed.
  z_values <- numeric()
  for (d in saturdays) {
    for (loc in locations) {
      r <- .standardized_slope(train, as.Date(d, origin = "1970-01-01"), loc)
      if (!is.null(r) && is.finite(r$z)) {
        z_values <- c(z_values, abs(r$z))
      }
    }
  }

  # 67th percentile of |z|. Slopes whose |z| exceeds this are in the
  # "clearly trending" top third — empirically calibrated, not assumed.
  slope_z_threshold <- as.numeric(stats::quantile(z_values, 0.67, na.rm = TRUE))

  # ----- Level thresholds (per location) -----
  level_thresholds <- train |>
    dplyr::group_by(location) |>
    dplyr::summarise(
      p10 = stats::quantile(observation, 0.10, na.rm = TRUE),
      p90 = stats::quantile(observation, 0.90, na.rm = TRUE),
      .groups = "drop"
    )

  # ----- Calendar fallback (per location, per week of year) -----
  # For each (location, week_of_year), compute mean wILI and rank that
  # mean within the location's distribution. Use the rank to label.
  calendar_means <- train |>
    dplyr::mutate(woy = .week_of_year(target_end_date)) |>
    dplyr::group_by(location, woy) |>
    dplyr::summarise(mean_wili = mean(observation, na.rm = TRUE),
                     .groups = "drop")

  # Derive regime label per (location, woy) from where mean_wili sits in
  # the location's woy distribution AND the slope of that mean across
  # adjacent weeks (which captures whether the season is rising/falling
  # at this point in the calendar).
  location_calendar <- calendar_means |>
    dplyr::group_by(location) |>
    dplyr::arrange(woy, .by_group = TRUE) |>
    dplyr::mutate(
      # Lead/lag week means for slope calculation.
      mean_wili_prev = dplyr::lag(mean_wili, default = mean_wili[1]),
      mean_wili_next = dplyr::lead(mean_wili, default = mean_wili[dplyr::n()]),
      # Calendar slope in units of wILI per week.
      cal_slope = (mean_wili_next - mean_wili_prev) / 2,
      # Quantile of this woy's mean within all woy means for the location.
      level_pct = stats::ecdf(mean_wili)(mean_wili),
      # Regime label: peaking if level is high, off_season if low,
      # otherwise growing/declining based on calendar slope sign.
      regime = dplyr::case_when(
        level_pct >= 0.90              ~ "peaking",
        level_pct <= 0.20              ~ "off_season",
        cal_slope > 0                  ~ "growing",
        cal_slope < 0                  ~ "declining",
        TRUE                           ~ "off_season"
      )
    ) |>
    dplyr::select(location, week_of_year = woy, regime, mean_wili) |>
    dplyr::ungroup()

  list(
    slope_z_threshold = slope_z_threshold,
    level_thresholds  = level_thresholds,
    location_calendar = location_calendar
  )
}


# ============================================================================
# PUBLIC: CLASSIFICATION
# ============================================================================

#' Classify the regime of one (location, as_of_date) cell.
#'
#' Runs the 4-step cascade. First step that produces a non-NA classification
#' wins. The cascade is designed so steps 1-2 are conservative (fall
#' through often) and steps 3-4 are guaranteed-non-null backstops.
#'
#' @param wili_data Long tibble: location, target_end_date, observation.
#' @param as_of_date Date. The forecasting origin date.
#' @param location Character. Location to classify.
#' @param calibration List returned by compute_regime_calibration().
#' @return Character: "growing", "peaking", "declining", or "off_season".
classify_regime <- function(wili_data, as_of_date, location, calibration) {

  # ----- Step 1: Slope test on this location -----
  loc_slope <- .standardized_slope(wili_data, as_of_date, location)
  if (!is.null(loc_slope) &&
      is.finite(loc_slope$z) &&
      abs(loc_slope$z) > calibration$slope_z_threshold) {
    return(if (loc_slope$z > 0) "growing" else "declining")
  }

  # ----- Step 2: Level test on this location -----
  recent <- .recent_mean(wili_data, as_of_date, location)
  if (!is.na(recent)) {
    th <- .level_thresholds(calibration, location)
    if (recent >= th$p90) return("peaking")
    if (recent <= th$p10) return("off_season")
  }

  # ----- Step 3: National tiebreaker -----
  # Run steps 1-2 on US National. If they classify, apply that label here.
  nat_slope <- .standardized_slope(wili_data, as_of_date, "US National")
  if (!is.null(nat_slope) &&
      is.finite(nat_slope$z) &&
      abs(nat_slope$z) > calibration$slope_z_threshold) {
    return(if (nat_slope$z > 0) "growing" else "declining")
  }

  nat_recent <- .recent_mean(wili_data, as_of_date, "US National")
  if (!is.na(nat_recent)) {
    nat_th <- .level_thresholds(calibration, "US National")
    if (nat_recent >= nat_th$p90) return("peaking")
    if (nat_recent <= nat_th$p10) return("off_season")
  }

  # ----- Step 4: Calendar fallback -----
  # Always returns a label — last-resort backstop never falls through.
  .calendar_regime(calibration, location, as_of_date)
}


# ============================================================================
# PUBLIC: BATCH PRECOMPUTATION
# ============================================================================

#' Precompute regime labels for many (origin_date, location) cells at once.
#'
#' Builds a lookup table the rest of the pipeline can join against. Saves
#' repeated calls to classify_regime() inside Phase 3's BMA loop.
#'
#' @param wili_data Long tibble: location, target_end_date, observation.
#' @param origin_dates Date vector of origin dates to classify.
#' @param locations Character vector of locations.
#' @param calibration List from compute_regime_calibration().
#' @return Tibble with columns origin_date, location, regime.
precompute_regimes <- function(wili_data, origin_dates, locations,
                                calibration) {

  # Build all (date, location) combinations.
  grid <- tidyr::expand_grid(
    origin_date = origin_dates,
    location    = locations
  )

  # Classify each row.
  grid$regime <- mapply(
    function(d, loc) classify_regime(wili_data, d, loc, calibration),
    grid$origin_date, grid$location,
    USE.NAMES = FALSE
  )

  tibble::as_tibble(grid)
}
