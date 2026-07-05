# code/ — simulation source (reference snapshot)

Snapshot of the code that produced `../data/`. The canonical *live* copies are in the project
root; these match them at the time of this package.

## Files & dependency order
Source them in this order (each depends on the ones above):

| File | What it is | Depends on |
|------|-----------|------------|
| `app.R` | The **v5 model** — all primitives (production `λ`, knowledge growth, Pareto/copula population draw, importance-sampling posteriors, the 2-round planners) plus a Shiny UI. **The `ce_reweight_posterior` fix lives here** (deterministic reweight; see STATE_OF_PLAY §5). | base R only |
| `simulate_T.R` | The **T-round core**: generalizes the planners to T rounds (`run_simulation_T`, the receding-horizon forward planner, per-round myopic/naive/seed, schedule metrics). | `app.R` primitives |
| `sweep.R` | Base sweep infrastructure: `summarize_sweep`, ggplot helpers, the original 2-round configs. | `app.R` |
| `sweep_T.R` | The **T-round manifest**: `SWEEP_BASE_PARAMS_T` (base defaults), `SWEEP_CONFIGS_T` (the 13 sweeps), `run_sweep_T`, `sweep_one_T`, `main_sweep_T`. | `app.R`, `simulate_T.R`, `sweep.R` |

## How to load (without launching Shiny)
`app.R` ends with a `shinyApp(...)` call; strip it before sourcing programmatically:

```r
src <- readLines("app.R", warn = FALSE)
eval(parse(text = paste(src[1:(grep("^shinyApp\\(", src)[1] - 1L)], collapse = "\n")),
     envir = globalenv())
source("simulate_T.R"); source("sweep.R"); source("sweep_T.R")
```

## Key entry points
- `run_simulation_T(seed, T_rounds, n, epsilon, b, tau_k, ..., strategies = 1:9)` — one trial,
  returns per-strategy `total_expected`, `g_rounds` (schedule), `alpha`, `b_idx`, `gini_g1`, and
  `rho_s`. This is the model.
- `sweep_one_T("<name>", seeds = 1:200)` — run one sweep, write `{summary, raw, rawlong}.rds` + plots.
- `main_sweep_T(seeds = 1:200)` — run the whole manifest (checkpointed: skips sweeps already saved).

## Notes for re-running
- Scripts `setwd()` to the **project root** and read `app.R` there — run from the project root,
  not from inside this package, or adjust the paths.
- Parallelism: `parallel::mclapply` (fork), cores = `detectCores()-1`. Reproducible because each
  run seeds its own RNG from the `seed` argument.
- `n_steps` is the greedy granularity (`δ = B/n_steps`). Results are now granularity-stable, so 50
  is fine; the diagnostics use up to 800 to confirm convergence.
