# New Candidate Model Configs + Analysis Scripts

## Files to copy into your ATSF2026 repo

### Config files → `configs/`

| File | model_abbr | Model | Transform | Runtime/date |
|------|-----------|-------|-----------|-------------|
| `arima_bc_bs.R` | arima_bc_bs | Auto-ARIMA(p,d,q)(P,D,Q)[52] | Box-Cox guerrero | ~15-25 sec |
| `nnetar_bc_bs.R` | nnetar_bc_bs | NNETAR neural network | Box-Cox guerrero | ~20-40 sec |
| `theta_bc_bs.R` | theta_bc_bs | Theta method (M3 winner) | Box-Cox guerrero | ~5-10 sec |
| `tslm_fourier.R` | tslm_fourier | Regression + 5 Fourier pairs | log(x+0.001) | ~2-5 sec |
| `stl_arima_bc.R` | stl_arima_bc | STL decomposition + ARIMA | Box-Cox guerrero | ~15-25 sec |

### Script files → `scripts/`

| File | Purpose |
|------|---------|
| `run_all_candidates.R` | Batch runner: sources each config and runs full loop |
| `analyze_candidates.R` | Error correlation analysis + ensemble selection |

## Quick Start

### Step 1: Copy files

```bash
cd /path/to/ATSF2026
cp <downloaded>/configs/*.R  configs/
cp <downloaded>/scripts/*.R  scripts/
```

### Step 2: Quick test (5 dates per model, ~10 minutes total)

```r
# Edit run_all_candidates.R: set QUICK_TEST <- TRUE
source("scripts/run_all_candidates.R")
```

### Step 3: Full run (all dates, multiple hours)

Recommended: run in parallel across 3-4 R sessions.

**Session 1 (fast models):**
```r
source("scripts/forecast_engine.R")
for (cfg in c("configs/hist_week.R", "configs/snaive_bc_bs.R",
              "configs/tslm_fourier.R", "configs/theta_bc_bs.R")) {
  source(cfg); config$run_loop <- TRUE; run_forecast(config)
}
```

**Session 2:**
```r
source("scripts/forecast_engine.R")
source("configs/arima_bc_bs.R"); config$run_loop <- TRUE; run_forecast(config)
```

**Session 3:**
```r
source("scripts/forecast_engine.R")
source("configs/nnetar_bc_bs.R"); config$run_loop <- TRUE; run_forecast(config)
```

**Session 4:**
```r
source("scripts/forecast_engine.R")
source("configs/stl_arima_bc.R"); config$run_loop <- TRUE; run_forecast(config)
```

### Step 4: Analyze

```r
source("scripts/analyze_candidates.R")

# Save plots for your presentation slides
ggsave("slides/correlation_heatmap.png", p_heatmap, width = 8, height = 7)
ggsave("slides/model_accuracy.png", p_accuracy, width = 8, height = 5)
```

## Troubleshooting

### ARIMA is very slow
The auto-ARIMA search with period=52 can be slow. If a single date
takes >60 seconds, you can constrain the search in the config:
```r
model_fn = function(lambda) {
  ARIMA(box_cox(observation, lambda) ~ pdq(0:2, 0:1, 0:2) + PDQ(0:1, 0:1, 0:1, period = 52))
}
```

### STL + ARIMA fails with generate()
decomposition_model() with bootstrap generation can have edge cases.
If it errors, switch to parametric intervals:
```r
simulate_method = "parametric"
```
Or convert to the custom_simulate_fn path (see hist_week.R for the pattern).

### NNETAR produces warnings about convergence
This is normal — some of the 20 random-start networks may not converge.
The averaging across networks handles this gracefully. Warnings can be
suppressed with suppressWarnings() (the engine already does this in the
forecast loop).

### A model fails on early origin dates
Models requiring 104+ training weeks will skip origin dates that don't
have enough history. This is handled by the engine's min_train_weeks
check — it prints a warning and returns NULL. The analyze_candidates.R
script uses inner_join, so it only analyzes dates where all models
produced forecasts.
