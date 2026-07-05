# T-Round Grant-Allocation Model — Analysis Package

**What this is.** A self-contained snapshot of the T-round extension of the grant-funding /
scientific-output model: the simulation code, the full corrected parameter-sweep dataset
(13 sweeps, 172 cells, 200 trials each), the validation suite, the bug-hunt diagnostics, and
a visual results report. It is organized so someone **new to the project can understand the
model and start analyzing the data in ~15 minutes**.

> **If you only want to analyze the data:** everything you need is in **`data/`** (R `.rds`
> files) and **`DATA_DICTIONARY.md`**. The `.rds` files load with base R `readRDS()` — no
> project setup required. Start with `data/<sweep>_summary.rds` (one row per parameter cell).

---

## 1. The model in one page

Each trial simulates **n researchers** over **T funding rounds**. Researcher *i* has a latent
**knowledge** Kᵢ and **resources** Rᵢ (both Pareto-distributed). A funder allocates a fixed
budget across researchers and rounds under one of **9 strategies**, trying to maximize total
scientific output.

**Dynamics (per round t):**
- **Production (output rate):**  `λ(K,R) = γ·K·R / (K + R)`  — a researcher needs *both*
  knowledge and resources; whichever is scarcer bottlenecks output.
- **Grants are non-persistent:**  round-t resources are `R₀ + gₜ` (the grant `gₜ` applies to
  that round only; it does not accumulate).
- **Knowledge compounds:**  `K ← K + ε·K·R/(K+R)` after each round — funding a researcher this
  round permanently raises their knowledge (and thus future output). `ε` is the growth rate.
- **Observation:**  publications `pₜ ~ Poisson(λ)` — the funder sees noisy pub counts, plus two
  one-time signals: a **resource signal** σ_R and a **grant/peer-review signal** σ_K (noisy
  reads of R and K).

**The funder's problem** is Bayesian: infer each researcher's (K,R) from pubs + signals, then
allocate. The reported output metric is **expected** output `Σ λ` (noise-free given the state),
so results are low-variance.

### The 9 strategies
| # | Name | Uses grant signal? | Plans ahead? | Seed floor? |
|---|------|--------------------|--------------|-------------|
| S1 | No funding | – | – | – |
| S2 | Uniform | – | – | flat every round |
| S3 | Naive (∝ pubs) | – | – | – |
| S4 | Myopic (pubs) | no | no (this round only) | – |
| S5 | Myopic (pubs + grant) | **yes** | no | – |
| S6 | Myopic (pubs + seed) | no | no | round-1 floor |
| S7 | Forward (pubs) | no | **yes (full horizon)** | – |
| S8 | Forward (pubs + grant) | **yes** | **yes** | – |
| S9 | Forward (pubs + seed) | no | **yes** | round-1 floor |

**Myopic** planners allocate each round's budget tranche to maximize *that round's* output.
**Forward** planners are a receding-horizon **certainty-equivalent (CE)** Bayesian planner: they
plan the whole remaining horizon (anticipating knowledge compounding and, one step ahead, the
information their grants will reveal), execute the current round, then re-plan next round with
the newly observed pubs.

### The key contrasts (differences between strategies)
These are the analysis workhorses — each isolates the value of one ingredient:
| Contrast | Definition | Interpretation |
|----------|------------|----------------|
| `signal_fwd` | S8 − S7 | **value of the grant/peer-review signal** (under forward planning) |
| `signal_myo` | S5 − S4 | value of the grant signal (under myopic) |
| `fwd_vs_myo_PG` | **S8 − S5** | **value of planning ahead** (pubs+grant) — the horizon headline |
| `fwd_vs_myo_P` | S7 − S4 | value of planning ahead (pubs only) |
| `fwd_vs_myo_PS` | S9 − S6 | value of planning ahead (pubs+seed) |
| `seed_fwd` | S9 − S7 | value of a uniform seed floor (under forward) |
| `seed_myo` | S6 − S4 | value of a seed floor (under myopic) |
| `optimal_vs_naive` | S8 − S3 | best Bayesian planner over the naive baseline |

---

## 2. Headline findings (from the corrected data)

1. **Forward planning beats myopic across the board** (`fwd_vs_myo_PG > 0` in every cell), and
   the advantage **grows with knowledge-growth ε** — up to +1.66 at ε=0.85, T=5 (z≈24). The
   mechanism is **front-loading**: `b_idx_S8` (schedule center-of-mass) rises with ε (0.51→0.62).
   It decomposes into **compounding value** (∝ ε) + **information value** (small, survives ε→0).
2. **Grant-signal value is governed by inequality × precision.** `signal_fwd` reaches +21 under
   a heavy knowledge tail (α_K=1.3) and a sharp signal (τ_K=0.3), and falls to ~0 under a light
   tail and a noisy signal. Half-value at τ_K≈1.5.
3. **Signal dominance survives K–R correlation.** As ρ_c rises to +0.8 (realized ρ_s=0.78) the
   signal's value drops but stays large under heavy tails (22.9 → 15.6).
4. **Uniform seed floors never help** (`seed_fwd < 0` everywhere; worse at larger budgets).
5. **Robustness nulls:** resource-signal noise doesn't affect grant-signal value; no finite-n
   artifact; budget scale doesn't drive the forward gain (horizon does).

See `report/results_report.html` for the visual version, and `STATE_OF_PLAY.md` for the full
account (including the bug that was found and fixed — read that before trusting any *earlier*
numbers you may have seen).

---

## 3. Folder map

```
T_round_extension/
├── README.md              ← you are here (model + findings + how to use)
├── STATE_OF_PLAY.md       ← full status: what's done, the bug story, loose ends, next steps
├── DATA_DICTIONARY.md     ← every data file & column explained (READ THIS to analyze data)
├── NEXT_SIMULATIONS.md    ← proposed additional sweeps (grids + rationale, ready to run)
├── code/                  ← the simulation code (reference snapshot) + how it fits together
├── data/                  ← THE DATASET: 13 sweeps × {summary, raw, rawlong} .rds + plots/
├── validation/            ← tests + logs proving the code is correct (anchor, assertions, bench)
├── diagnostics/           ← the scripts that found & confirmed the bug (provenance)
└── report/                ← results_report.html — self-contained visual report
```

---

## 4. Quickstart — load and explore the data (R)

```r
# One row per parameter cell, with _mean and _se for every metric:
s <- readRDS("data/horizon_growth_summary.rds")

# The horizon headline: forward-vs-myopic gain across horizon T and growth eps
library(dplyr)
s %>%
  select(T_rounds, epsilon, fwd_vs_myo_PG_mean, fwd_vs_myo_PG_se) %>%
  mutate(z = fwd_vs_myo_PG_mean / fwd_vs_myo_PG_se) %>%
  arrange(epsilon, T_rounds)

# Trial-level data (200 rows per cell) for your own aggregation / stats:
raw <- readRDS("data/horizon_growth_raw.rds")   # 5000 rows = 25 cells x 200 trials

# Per-round funding schedule (long format), e.g. how forward front-loads:
sched <- readRDS("data/horizon_growth_rawlong.rds")
```

**"Statistically distinguishable from zero"** convention used throughout: `abs(mean) > 2*se`.
**Reproducibility:** trial `seed` = the trial index, shared across cells (common random numbers),
so cross-cell contrasts are low-variance. `set.seed()` is applied per-run inside the simulator.

---

## 5. Parameter glossary (symbol ↔ code name)

| Symbol | Code | Meaning | Base default |
|--------|------|---------|--------------|
| T | `T_rounds` | number of funding rounds (horizon) | 2 |
| n | `n` | number of researchers | 50 |
| ε | `epsilon` | knowledge-growth rate per round | 0.1 |
| α_K | `k_shape` | knowledge Pareto tail exponent (**smaller = heavier tail / more inequality**) | 2 |
| α_R | `r_shape` | resource Pareto tail exponent | 2 |
| τ_K | `tau_k` | grant/peer-review signal noise SD (**smaller = sharper**) | 1 |
| τ_R | `tau_r` | resource signal noise SD | 1 |
| ρ_c | `rho_kr` | K–R Gaussian-copula correlation (report realized rank corr `rho_s`) | 0 |
| b | `b` | budget fraction; total budget `B_total = 2·b·n·E[R]`, split across T rounds | 0.5 |
| x_seed | `x_seed` | seed-floor fraction (S6, S9) | 0.25 |
| n_pre | `n_pre_rounds` | pre-round baseline pub observations | 0 |
| M | `M` | posterior importance samples | 400 |
| — | `n_steps` | greedy allocation granularity (`δ = B/n_steps`) | 50 |

Each sweep **varies a subset** of these and holds the rest at the base default. The exact
varied axes per sweep are in `DATA_DICTIONARY.md` §3.

---

*Generated as a housekeeping snapshot. Canonical live code lives in the project root
(`app.R`, `simulate_T.R`, `sweep_T.R`); this package's `code/` is a matching copy.*
