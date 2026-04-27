# FluCast — Phase 2–3 Design Decisions

This document records the methodological decisions made between Phase 1
completion and the start of Phase 3 ensemble implementation. Each
decision lists the alternatives considered, the choice made, and the
reasoning. Intended audience: future-self at report-writing time;
graders or reviewers reproducing the work.

**Project context**: INF 599 final submission for the FluSight ILI Sandbox
Hub (ATSF2026). Validation = 57 origin dates from 2015-10-24 through
2017-05-06 (Seasons 1–2). Test = 77 dates from 2017-10-28 through
2020-02-29 (Seasons 3–5).

**Implementation plan reference**: `FLUCAST_IMPLEMENTATION_PLAN.md`
(project root). All phase numbers below refer to that plan.

---

## Phase 1: Candidate Pool Expansion

### 1.1 Skipping Prophet

**Considered**: Add `fable.prophet` or raw `prophet` package as an
ensemble candidate.

**Decision**: Skip Prophet entirely.

**Reasoning**:
- `fable.prophet` 0.x does not provide a `generate.fbl_prophet`
  S3 method, so the engine's bootstrap simulation path doesn't dispatch.
  Verified by inspecting `methods("generate")` after install.
- The raw `prophet` package would have required a custom `simulate_fn`
  similar to the bsts integration (~80 lines).
- bsts already covers the "MCMC-based Bayesian uncertainty" niche that
  Prophet would have filled.
- Recent FluSight literature reports Prophet underperforming simpler
  seasonal baselines on weekly ILI data due to its automatic
  changepoint detection over-fitting short series.
- The implementation plan explicitly authorized this exit
  ("fallback: skip if cmdstanr fails"). Skipping was always on the table.

### 1.2 Off-by-one in season length count

**Found**: Implementation plan documented "56 validation dates (29 + 27)";
actual count from `tasks.json` is 57.

**Cause**: Counting Saturday-to-Saturday inclusive in Season 2016-17:
2016-10-29 to 2017-05-06 spans 28 weeks, not 27. Inclusive endpoint
issue in the plan.

**Resolution**: Updated `scripts/run_validation.R` to derive the count
from `tasks.json` rather than hardcode it. The validation set is 57
dates; test set is correspondingly 77.

---

## Phase 2: Scoring & Selection

### 2.1 Scoring library

**Considered**: Hand-implement WIS vs. use `scoringutils`.

**Decision**: `scoringutils` 2.2.0.

**Reasoning**: It's the FluSight community standard, battle-tested,
ensures reproducibility against other teams' analyses. Compatible with
the v2 `as_forecast_quantile() → score()` workflow.

### 2.2 95% coverage computed manually

**Found**: `scoringutils::score()` returns 50% and 90% interval coverage
by default, not 95%.

**Resolution**: 95% coverage flag computed manually from q0.025 and
q0.975 quantiles. FluSight reports use 95%, so the manual computation
keeps results comparable to the literature.

### 2.3 Correlation methodology

**Considered**: 
- Flatten across (location, horizon, date) — single 2508-element
  correlation per pair
- Per-(location, horizon) correlations averaged across 44 cells

**Decision**: Per-(location, horizon) average, matching the proposal-era
methodology.

**Reasoning**: Locations have different scales (Region 6 baseline ~3%,
Region 10 baseline ~0.5%). Horizons have different error magnitudes (h=1
WIS ~0.15, h=4 WIS ~0.50). Computing correlations within homogeneous
cells before averaging avoids these scale differences masking the actual
inter-model relationships.

### 2.4 Both point-error and WIS correlations reported

**Decision**: Compute and report both.

**Reasoning**: Point-error correlation is the standard diversity metric
in the literature. WIS-based correlation is more methodologically
appropriate for our use case (we're averaging quantile distributions,
not point forecasts). They typically agree; where they diverge, that's
informative. Reporting both costs ~5 lines of code and provides a
robustness check.

---

## Phase 2 Twin-Pair Verdicts

Three model pairs had point-error correlation > 0.93, triggering the
twin-pair rule from the proposal-era methodology. For each, we ran a
per-(horizon, regime) WIS comparison with a 2% relative-WIS margin to
declare cell winners, breaking ties via mean dispersion.

| Pair | Verdict | Margin |
|---|---|---|
| `ets_log` vs `ets_bc` | **KEEP `ets_bc`** | 9-0 with 3 ties |
| `arima_bc_bs` vs `arima_log` | **KEEP `arima_bc_bs`** | 10-0 with 2 ties |
| `nnetar_bc_bs` vs `nnetar_log` | **KEEP `nnetar_bc_bs`** | 3-1 with 8 ties; wins all peaking-regime cells |

**Notable finding for the report**: `nnetar_bc_bs` is meaningfully better
than `nnetar_log` in peaking-regime cells (h=1, h=2, h=3 peaking all go
to `nnetar_bc_bs`). The aggregate verdict is a coin flip, but the
regime-conditional pattern shows `nnetar_bc_bs` is sharper exactly when
forecast quality matters most for public-health decisions.

---

## Phase 2 Selection Philosophy

### 2.5 Ensemble pool composition

**Considered**:
- Continuity-first: keep the proposal's 5-model ensemble (NNETAR + ARIMA
  + ETS + SNAIVE + HistWeek), explain swaps in report
- Data-driven: use whatever WIS+correlation analysis says is best
- Hybrid: data-driven core + diversity anchors

**Decision**: Data-driven core + all 4 diversity anchors at low default
weight (BMA decides).

**Reasoning**: The aggregate WIS leaderboard (top 4: nnetar_bc_bs,
stl_arima_bc, arima_bc_bs, ets_bc) hides the possibility that
"badly-scoring" models are actually good in specific subsets (a regime,
a location). Including the diversity anchors and trusting the
regime-conditional BMA weighting to assign them low weights in normal
regimes—but potentially higher weights in unusual regimes—maximizes the
chance of catching that signal. If the BMA weighting works as
designed, the anchors cost nothing; if a regime exists where one anchor
genuinely outperforms, including them captures that gain.

**Final 8-model shortlist**: nnetar_bc_bs, stl_arima_bc, arima_bc_bs,
ets_bc, hist_week, bsts_seasonal, tslm_fourier, snaive_bc_bs.

---

## Phase 2 Regime Classifier

### 2.6 Regime scope

**Decision**: Per-location regimes plus US National signal as
tiebreaker.

**Reasoning**: Locations have meaningfully different seasonal dynamics
(timing, amplitude). National signal usually leads regional, so when
local signal is ambiguous, national context is informative.

### 2.7 As-of historical data

**Decision**: Use only data with `target_end_date < origin_date` (or
`< validation_start_date` for calibration).

**Reasoning**: "Honest as-of forecasting" — no future information leaks
into regime labels. Defensible to any reviewer; reproducible by anyone
with the same target-data file.

### 2.8 Threshold derivation

**Considered**: Hardcoded thresholds (e.g., "slope > 0.05") vs.
data-driven calibration.

**Decision**: Both step-1 slope threshold and step-4 calendar boundaries
derived from training data (pre-2015).

**Reasoning**: Hardcoded thresholds are the first thing a reviewer would
flag. Data-driven thresholds are calibrated to the actual signal-to-noise
of this dataset. Specifically:
- **Slope threshold**: 67th percentile of the distribution of
  standardized slopes (|z| where z = slope / rolling SD) computed across
  all (location, week) cells in the pre-2015 training history.
  Empirically calibrated value: ~0.702.
- **Calendar boundaries**: per-(location, week-of-year) regime labels
  derived from each location's empirical week-of-year wILI percentile
  and slope.

### 2.9 Regime distribution

Validation period (627 cells = 57 dates × 11 locations):

| Regime | Cells | Percentage |
|---|---|---|
| growing | 346 | 55% |
| declining | 201 | 32% |
| peaking | 80 | 13% |
| off_season | 0 | 0% |

Despite the calendar-fallback distribution including off_season cells
(110 in pre-validation training), no validation-period cells fall
through to off_season — the slope and level checks always find signal.
This is good news for BMA: only 3 active regimes to estimate weights
for, with at least 80 cells in each.

---

## Phase 3: Ensemble Implementation Decisions

### 3.1 BMA weight granularity

**Considered**:
- Pure regime-conditional: weights depend only on regime (3 regimes)
- Regime + horizon: weights depend on both (3 × 4 = 12 combinations)
- Hybrid: regime-primary with horizon smoothing

**Decision**: Pure regime-conditional.

**Reasoning**:
- **Sample size adequacy**: regime-only gives 27–115 cells per (model,
  regime); regime+horizon would split that across 4 horizons → 7–29
  cells per (model, regime, horizon). 7 cells is well below the rule-of-
  thumb 30 observations per parameter for stable estimation.
- **The forecasts already encode horizon information correctly**: each
  model's h=4 quantiles are appropriately wider than h=1 quantiles.
  The ensemble doesn't need separate per-horizon weights to handle
  horizon — it needs the components' horizon-specific outputs to be
  combined faithfully.
- **Defensibility**: weights depend on regime because regime materially
  changes which model is best (e.g., `nnetar_bc_bs` dominates peaking).
  Horizon-conditional refinement was rejected on sample-size grounds.
- **Occam's razor**: simpler model wins absent strong evidence to use
  the more complex one.

### 3.2 Quantile aggregation method

**Considered**:
- Linear pool of quantile values (FluSight standard)
- Vincentization (mathematically identical to linear pool for our case)
- Sample mixing: draw N samples from each component proportional to
  weight, recompute quantiles

**Decision**: Linear pool primary; sample mixing as a *post-hoc*
validation check.

**Reasoning**:
- Linear pool is the FluSight community standard. Using it ties our
  methodology to community practice and makes our results comparable.
- For unimodal, similarly-shaped components on the same support
  (which our 8 components are), linear pooling and sample mixing should
  agree to within ~1% WIS. Implementing both as primary would be
  redundant.
- Sample mixing has the theoretical advantage of producing a coherent
  posterior predictive distribution, but our task only requires
  quantile outputs, so the advantage doesn't matter operationally.
- Post-hoc check: implement sample mixing in Phase 3C/4. If results
  agree to within ~1% WIS, report linear pool. If they meaningfully
  diverge, that's a finding worth investigating.

### 3.3 Weight estimation approach

**Considered**:
- Validation-only frozen weights (matches plan)
- Expanding window during test season
- Hybrid: frozen with periodic regime-specific updates

**Decision**: Validation-only frozen weights.

**Reasoning**:
- **No leakage**: weights computed from data strictly before any test
  date.
- **Maximum learning from validation**: every validation cell
  contributes to weight estimation (627 cells × 8 models = 5016
  observations).
- **Reproducible**: a third party could verify the exact weights from
  validation data.
- **Honest in the report**: "Weights were learned from 627 cells across
  2 validation seasons and applied unchanged to test seasons" is a
  clean methodological statement.
- Expanding-window methods are theoretically better but add complexity
  that's hard to defend without empirical evidence the simpler method
  is insufficient. We can't generate that evidence without violating
  the "test-set-unseen" principle.

---

## File Map

| Path | Purpose |
|---|---|
| `analysis/phase2/scores_long.csv` | Master score table (27,588 rows) |
| `analysis/phase2/per_model_summary.csv` | Aggregate WIS/MAE/coverage per model |
| `analysis/phase2/correlation_point_error.csv` | 11×11 correlation matrix |
| `analysis/phase2/correlation_wis.csv` | 11×11 WIS-based correlation matrix |
| `analysis/phase2/regimes.csv` | Regime label per (origin_date, location) |
| `analysis/phase2/twin_diagnostics.csv` | Per-(pair, horizon, regime) WIS comparison |
| `analysis/phase2/twin_verdicts.csv` | Final twin-pair decisions |
| `scripts/score_candidates.R` | Phase 2A+2B scoring script |
| `scripts/regime_helpers.R` | Regime classifier |
| `scripts/twin_diagnostics.R` | Twin-pair diagnostic |

---

*Last updated: end of Phase 2, before Phase 3 ensemble implementation.*
