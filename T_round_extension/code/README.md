# code/ — simulation source (snapshot)

Snapshot of the code that produced `../data/`, matching the repo-root live copies.

## Files & dependency order
Source in this order (each depends on the ones above):

| File | What it is | Depends on |
|------|-----------|------------|
| `model.R` | The **model engine** (pure base R): all primitives (production `λ`, knowledge growth, Pareto/copula population draw, importance-sampling posteriors, the fixed deterministic `ce_reweight_posterior`), the 2-round `run_simulation()` reference, and the **T-round** `run_simulation_T()` + receding-horizon forward planner. | base R only |
| `app.R` | The **Shiny app** — presentation layer: UI, server, plots; `source("model.R")` + a `T` rounds slider. | `model.R`, shiny/ggplot2 |
| `sweep.R` | Base sweep infrastructure: `summarize_sweep`, ggplot helpers. | `model.R` |
| `sweep_T.R` | The **T-round manifest**: `SWEEP_BASE_PARAMS_T` (base defaults), `SWEEP_CONFIGS_T` (the 16 sweeps), `run_sweep_T`, `sweep_one_T`, `main_sweep_T`. | `model.R`, `sweep.R` |

> Note: the former `simulate_T.R` (a separate T-round file that sourced `app.R`'s primitives) has
> been absorbed into `model.R`; the model is now one self-contained base-R engine.

## How to load
```r
source("model.R")                       # the model engine (no Shiny/ggplot needed)
source("sweep.R"); source("sweep_T.R")  # sweep infrastructure (needs ggplot2/dplyr/tidyr)
```

## Key entry points
- `run_simulation_T(seed, T_rounds, n, epsilon, b, tau_k, ..., strategies = 1:9, detail = FALSE)` —
  one trial; returns per-strategy `total_expected`, `g_rounds` (schedule), `alpha`, `b_idx`,
  `gini_g1`, and `rho_s`. With `detail = TRUE` it also returns per-round K/R/λ/p trajectories (used
  by the Shiny app). This is the model.
- `run_simulation(seed, ...)` — the 2-round v5 reference; `run_simulation_T(T_rounds = 2)` reproduces
  it bit-identically (the regression anchor).
- `sweep_one_T("<name>", seeds = 1:200)` — run one sweep, write `{summary, raw, rawlong}.rds` + plots.
- `main_sweep_T(seeds = 1:200)` — run the whole 16-sweep manifest (checkpointed).

## Notes for re-running
- Scripts `setwd()` to the **project root** — run from there.
- `horizon_long` pins `n_steps=400` in its config; all other sweeps use the base `n_steps=50`
  (granularity-stable for T≤5; see `../DATA_DICTIONARY.md`).
- Parallelism: `parallel::mclapply` (fork), reproducible because each run seeds its own RNG.
