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

**The current, canonical model is the T-round model** (`simulate_T.R`): the funder allocates over
**T ∈ {1..5+} rounds**. Everything else is a **special case** of it:

| Model | File | Relationship |
|-------|------|--------------|
| **T-round** (default) | `simulate_T.R` | the general model — funds over T rounds |
| v5 (2-round) | `app.R` | the **T=2 special case**; the current interactive Shiny app; reproduced *bit-identically* by the T-round model at T=2 |
| v3 (2-round, older planner) | `archive/models/simulation_v3.R` | superseded earlier version, kept for history |

> **In progress:** `app.R` (the Shiny app) is currently the v5 2-round model. The next step is to
> upgrade it to the T-round model so the interactive app matches the default. Until then, the
> T-round model lives in `simulate_T.R` and is exercised through the sweep code and tests.

**Headline result (from the corrected data):** forward planning beats myopic across the board, and
by more as knowledge compounds faster (higher ε) and the horizon lengthens — via *front-loading*
grants to compound knowledge. The value of the peer-review signal is set by *inequality × precision*.
Full findings: [`T_round_extension/README.md`](T_round_extension/README.md).

---

## Start here

| I want to… | Go to |
|-----------|-------|
| **Analyze the sweep data** | [`T_round_extension/`](T_round_extension/) — self-contained package (data + dictionary + docs). Read its `README.md` and `DATA_DICTIONARY.md`. |
| **Understand project status / history** | [`docs/PROGRESS.md`](docs/PROGRESS.md) (running log) and [`T_round_extension/STATE_OF_PLAY.md`](T_round_extension/STATE_OF_PLAY.md) (clean status, the bug story, loose ends). |
| **See what to simulate next** | [`T_round_extension/NEXT_SIMULATIONS.md`](T_round_extension/NEXT_SIMULATIONS.md). |
| **Run the model / a sweep** | `simulate_T.R` + `sweep_T.R` (see "Running" below). |
| **Run the interactive app** | `app.R` (Shiny). |

---

## Directory map

```
Grant-Funding-and-Scientific-Output-Model/
│
├── README.md              ← this file (project overview + organization)
│
├── app.R                  ← Shiny app (v5 2-round; being upgraded to T-round). Deployment entry point.
├── simulate_T.R           ← THE DEFAULT MODEL — the T-round simulator & planners
├── sweep.R                ← shared sweep infrastructure (aggregation, plotting)
├── sweep_T.R              ← the T-round sweep manifest (13 sweeps) + runner
│                            (these four are the live source; app.R sources the others' primitives)
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
│   ├── T_run_fixed/           the CURRENT corrected T-round manifest (13 sweeps, 200 trials, M=400)
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
# Load the model without launching Shiny (app.R holds the shared primitives):
src <- readLines("app.R", warn = FALSE)
eval(parse(text = paste(src[1:(grep("^shinyApp\\(", src)[1] - 1L)], collapse = "\n")),
     envir = globalenv())
source("simulate_T.R"); source("sweep.R"); source("sweep_T.R")

run_simulation_T(seed = 1, T_rounds = 3, n = 50, epsilon = 0.3, strategies = 1:9)  # one trial
sweep_one_T("horizon_growth", seeds = 1:200)                                        # one sweep
main_sweep_T(seeds = 1:200)                                                         # full manifest (~30 min)
```

Reproduction details and the exact production driver are in `T_round_extension/validation/`.

---

## Deployment

The interactive app is `app.R`, deployed to shinyapps.io (`rsconnect/`). **Note:** the
`ce_reweight_posterior` bug fix (see `STATE_OF_PLAY.md` §5) changed `app.R`'s forward planner — the
live app should be **redeployed** so it matches these results.
