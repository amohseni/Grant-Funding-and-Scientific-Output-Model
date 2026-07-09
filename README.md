# Grant Funding & Scientific Output Model

A simulation model of how a **science funder should allocate grants** across researchers and over
time to maximize scientific output — and, in particular, whether a **forward-looking Bayesian
planner** (which learns from publications and anticipates how today's grants grow tomorrow's
capacity) beats simpler **myopic** or **naive** allocation rules.

Each researcher has latent **knowledge** K and **resources** R; output is `λ = γ·K·R/(K+R)`
(you need both). Funding raises resources this round and **compounds knowledge** for future rounds.
The funder infers each researcher's (K,R) from noisy publication counts and two signals, then
allocates a fixed budget under one of **9 strategies** (no-funding → naive → myopic → forward).

---

## The default model: the T-round extension

**The current, canonical model is the T-round model** (`model.R`): the funder allocates over
**T ∈ {1..5+} rounds**. Everything else is a **special case** of it:

| Model | File | Relationship |
|-------|------|--------------|
| **T-round** (default) | `model.R` | the model **engine** — all simulation logic (`run_simulation_T` etc.), pure base R, sourced by everything |
| Interactive app | `app.R` | Shiny front-end for the T-round model (`source("model.R")` + UI/server/plots); a `T` slider selects the number of rounds |
| v5 (2-round) | `model.R` → `run_simulation()` | the **T=2 special case**, kept in the engine as the reference the T-round path reproduces *bit-identically* at T=2 |
| v3 (2-round, older planner) | `archive/models/simulation_v3.R` | superseded earlier version, kept for history |

> The model engine and the Shiny app are factored apart: **`model.R`** holds the model (no Shiny),
> **`app.R`** is the presentation layer. The app exposes rounds `T ∈ 1..5` (default 2). The former
> `simulate_T.R` has been absorbed into `model.R`.

**Headline result (from the corrected data):** forward planning beats myopic across the board, and
by more as knowledge compounds faster (higher ε) and the horizon lengthens — via *back-loading*
grants (lean early, heavy late) so resources arrive once knowledge has compounded up. The value of the
peer-review signal is set by *inequality × precision*.
Full findings: [`T_round_extension/README.md`](T_round_extension/README.md).

---

## Start here

| I want to… | Go to |
|-----------|-------|
| **Analyze the sweep data** | [`T_round_extension/`](T_round_extension/) — self-contained package (data + dictionary + docs). Read its `README.md` and `DATA_DICTIONARY.md`. |
| **Read the results** | [`RESULTS.md`](RESULTS.md) — statistical digest by claim; figures in [`figures/`](figures/) (regenerate with `Rscript figures/make_figures.R`). |
| **See the plan to publication** | [`ROADMAP.md`](ROADMAP.md) — steps to finish the study; [`docs/OSF_UPLOAD.md`](docs/OSF_UPLOAD.md) — the OSF archiving plan. |
| **Understand project status / history** | [`docs/PROGRESS.md`](docs/PROGRESS.md) (running log) and [`T_round_extension/STATE_OF_PLAY.md`](T_round_extension/STATE_OF_PLAY.md) (clean status, the bug story, loose ends). |
| **See what to simulate next** | [`T_round_extension/NEXT_SIMULATIONS.md`](T_round_extension/NEXT_SIMULATIONS.md). |
| **Run the model / a sweep** | `model.R` + `sweep_T.R` (see "Running" below). |
| **Run the interactive app** | `app.R` (Shiny; `source("model.R")` + UI). Rounds `T` is a slider (1–5). |

---

## Directory map

```
Grant-Funding-and-Scientific-Output-Model/
│
├── README.md              ← this file (project overview + organization)
│
├── model.R                ← THE DEFAULT MODEL ENGINE — run_simulation_T + planners (pure base R)
├── app.R                  ← Shiny app (T-round; sources model.R). Deployment entry point.
├── sweep.R                ← shared sweep infrastructure (aggregation, plotting)
├── sweep_T.R              ← the T-round sweep manifest (16 sweeps) + runner
│                            (app.R, sweep*.R and tests all `source("model.R")`)
│
├── T_round_extension/     ← ★ CURRENT WORK: self-contained analysis package
│   ├── README.md              model in one page, findings, quickstart
│   ├── STATE_OF_PLAY.md       status, the bug story, loose ends, next steps
│   ├── DATA_DICTIONARY.md     every data file & column
│   ├── NEXT_SIMULATIONS.md    proposed additional sweeps
│   ├── code/  data/  validation/  diagnostics/  report/
│                            (its own snapshot of code + the corrected dataset + docs;
│                             portable — meant to be analyzed on its own)
│
├── sweep_results/         ← raw parameter-sweep outputs
│   ├── T_run_fixed/           the CURRENT corrected T-round manifest (16 sweeps, 200 trials, M=400)
│   └── legacy/                superseded runs (see sweep_results/README.md):
│       ├── T_round_buggy_M200/    pre-bugfix T-round run — DO NOT USE
│       ├── sweep_2round_11-05-2026/  old 2-round (v5) sweep
│       ├── micro_run_12-04-2026/     small early iteration
│       ├── plots_DHA_*/, old_csv/, logs/, pg_focused_summary.rds
│
├── docs/                  ← process & metadata docs
│   ├── PROGRESS.md            running development log
│   ├── CHAT_HANDOFF_PROMPT.md older analysis-handoff prompt (historical; its CSV paths now under sweep_results/legacy/old_csv/)
│   ├── assertions_log.txt     validation output (current model — all checks pass)
│   └── benchmark_report.txt   runtime benchmark
│
├── archive/               ← older model versions & scripts (kept for history)
│   ├── models/simulation_v3.R    the v3 Shiny app
│   └── sweep_scripts/            older sweep scripts (parameter_sweep.R, sweep_11-05-2026.R, sweep_12-04-2026.R)
│
├── tests/                 ← current test & diagnostic scripts (anchor test, assertions, the bug-hunt)
├── rsconnect/             ← shinyapps.io deployment config
└── (.git, .Rhistory, .RData, .claude — housekeeping)
```

**Which data is canonical?** `sweep_results/T_run_fixed/` (raw) — mirrored, documented, and
analysis-ready inside `T_round_extension/data/`. Everything under `sweep_results/legacy/` is
superseded; the `T_round_buggy_M200` run in particular must not be used (see the bug story in
`STATE_OF_PLAY.md`).

---

## Running the model (from this directory)

```r
source("model.R")                       # the model engine (pure base R)
source("sweep.R"); source("sweep_T.R")  # sweep infrastructure (need ggplot2/dplyr/tidyr)

run_simulation_T(seed = 1, T_rounds = 3, n = 50, epsilon = 0.3, strategies = 1:9)  # one trial
sweep_one_T("horizon_growth", seeds = 1:200)                                        # one sweep
main_sweep_T(seeds = 1:200)                                                         # full manifest (~30 min)
```

Launch the interactive app with `shiny::runApp("app.R")`.

Reproduction details and the exact production driver are in `T_round_extension/validation/`.

---

## Deployment

The interactive app is `app.R`, deployed to shinyapps.io (`rsconnect/`). It sources `model.R`, so
the deploy bundle must include both files. **Note:** the `ce_reweight_posterior` bug fix (see
`STATE_OF_PLAY.md` §5) changed `app.R`'s forward planner — the live app should be **redeployed** so
it matches these results.

---

## License & citation

- **Code** (all `.R` files) — MIT License (`LICENSE`).
- **Data & generated results** (`sweep_results/`, `T_round_extension/data/`, figures) — CC-BY-4.0
  (`LICENSE-DATA.md`).
- To cite, see `CITATION.cff`.

## Reproducibility

`docs/ENVIRONMENT.md` records the R and package versions. `reproduce.R` regenerates the entire
dataset from seeds and re-runs the validation suite (`Rscript reproduce.R`, ~45 min). The path from
here to publication is in `ROADMAP.md`.
