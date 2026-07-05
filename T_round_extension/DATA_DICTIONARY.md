# Data Dictionary

Everything in `data/`. All files are R serialized objects — load with `readRDS("data/<file>.rds")`,
which returns a `data.frame` (technically a tibble for summaries). No project setup needed.

Each of the 13 **sweeps** produces three files:

| File | Grain (one row = …) | Rows | Use it for |
|------|---------------------|------|------------|
| `<sweep>_summary.rds` | one **parameter cell**, aggregated over 200 trials | = # cells | main analysis, plotting, effect sizes |
| `<sweep>_raw.rds` | one **trial** (cell × seed) | # cells × 200 | your own aggregation, distributions, custom stats |
| `<sweep>_rawlong.rds` | one **(trial × strategy × round)** | # cells × 200 × (#strat × T) | funding **schedule** analysis (per-round allocation) |

Plus `data/plots/*.png` (quick-look ggplots the runner auto-generated) and `data/run.log`
(the console log of the production run: per-sweep timing).

---

## 1. `<sweep>_summary.rds` — cell-level (the main file)

**Columns, in order:**

- **Varied parameter(s)** — the axes this sweep varied (see §3). e.g. `horizon_growth` has
  `T_rounds`, `epsilon`.
- `n_trials` — trials aggregated into this cell (always **200**).
- For every metric `X` below, two columns: **`X_mean`** (mean over the 200 trials) and
  **`X_se`** (standard error = sd/√200).

**Metrics (the `X` in `X_mean`/`X_se`):**

*Per-strategy (j = 1..9):*
| Metric | Meaning |
|--------|---------|
| `out_Sj` | total **expected** output of strategy j (`Σ_t Σ_i λ`, summed over rounds & researchers) — the primary outcome |
| `alpha_Sj` | round-1 share of that strategy's **actual spend** (`sum(g₁)/sum(all g)`); NA for S1 (no spend) |

*Forward-planner schedule diagnostics (j ∈ {7,8,9} only):*
| Metric | Meaning |
|--------|---------|
| `b_idx_Sj` | schedule **center-of-mass** = `Σ_t t·αₜ / (T+1)`, where `αₜ = spendₜ/B_total`. **0.5 = uniform**; **>0.5 = front-loaded** (spends earlier); <0.5 = back-loaded |
| `gini_g1_Sj` | Gini coefficient of round-1 grants (concentration of funding; 0 = equal, →1 = concentrated). Also for S5. |

*Population diagnostic:*
| Metric | Meaning |
|--------|---------|
| `rho_s` | **realized** Spearman (rank) correlation between K and R in the drawn population. The knob is `rho_kr` (ρ_c); `rho_s` is the interpretable outcome. ≈0 for all sweeps except `correlation`. |

*Contrasts (differences that isolate one effect — see README §1 for interpretations):*
`fwd_vs_myo_P` (S7−S4), `fwd_vs_myo_PG` (S8−S5), `fwd_vs_myo_PS` (S9−S6),
`signal_myo` (S5−S4), `signal_fwd` (S8−S7), `seed_myo` (S6−S4), `seed_fwd` (S9−S7),
`optimal_vs_naive` (S8−S3).

> **Note:** the summary stores only the *varied* params, not the *fixed* ones. The fixed
> context for each sweep is the base default (§2), except for the axes it varies. This is
> documented per-sweep in §3 — keep it handy when comparing across sweeps.

---

## 2. `<sweep>_raw.rds` and `<sweep>_rawlong.rds`

**`_raw.rds`** — one row per trial. Columns: the varied param(s), `cell_id` (integer index of
the grid cell, 1..#cells), `seed` (trial index 1..200), then the **same metrics as the summary
but un-aggregated** (`out_Sj`, `alpha_Sj`, `b_idx_Sj`, `gini_g1_Sj`, `rho_s`, all contrasts —
no `_mean`/`_se`). Use for custom aggregation, paired tests (same `seed` across cells = common
random numbers), or output distributions.

**`_rawlong.rds`** — the funding **schedule**, one row per (trial × strategy × round):
| Column | Meaning |
|--------|---------|
| `cell_id` | grid cell index (join to `_raw` on `cell_id`; or map to params via §3) |
| `trial` | seed / trial index (1..200) |
| `strategy` | 1..9 |
| `round` | 1..T |
| `alpha_t` | fraction of **total budget** spent this round by this strategy (`spendₜ/B_total`) |
| `mean_g_t` | mean grant per researcher this round |
| `gini_g_t` | Gini of the round-t grant vector (funding concentration that round) |

Reconstruct any strategy's full spend schedule (α₁..α_T) by filtering `strategy` and pivoting on
`round`. This is where the front-loading story lives.

---

## 3. The 13 sweeps — varied axes, grid, and purpose

All at **200 trials/cell**. Unvaried parameters take the **base defaults** (§4). Symbol ↔ code
map is in README §5.

| Sweep | Cells | Varied axes (values) | Purpose / headline metric |
|-------|-------|----------------------|---------------------------|
| `horizon_growth` | 25 | `T_rounds`={1,2,3,4,5} × `epsilon`={0.05,0.15,0.3,0.55,0.85} | **The horizon headline.** `fwd_vs_myo_PG` across horizon × growth. |
| `horizon_noise` | 25 | `T_rounds`={1..5} × `tau_k`={0.1,0.3,1,3,10} | Horizon × signal noise. `signal_fwd`, `fwd_vs_myo_PG`. |
| `horizon_scale` | 20 | `T_rounds`={1..5} × `b`={0.1,0.3,0.5,1} | Horizon × budget. |
| `signal_precision` | 13 | `tau_k`={0.05,0.1,0.2,0.3,0.5,0.7,1,1.5,2,3,5,10,20} | Fine `signal_fwd` decay vs signal noise (T=2). |
| `signal_value` | 16 | `k_shape`={1.3,1.8,2.5,3.5} × `tau_k`={0.3,1,3,10} | **Inequality × precision.** `signal_fwd` (T=2). |
| `funder_scale` | 6 | `b`={0.1,0.2,0.3,0.5,0.75,1} | Budget scale (T=2). *(b-axis capped at 1.)* |
| `seed_value` | 16 | `b`={0.1,0.3,0.5,1} × `x_seed`={0.1,0.25,0.5,0.75} | `seed_fwd` — is seeding ever good? (T=2). |
| `alpha_regime` | 16 | `tau_k`={0.3,1,3,10} × `epsilon`={0.05,0.2,0.5,0.8} | Where forward spends: `b_idx_S8` (T=2). |
| `pre_rounds` | 5 | `n_pre_rounds`={0,1,5,10,20} | Effect of baseline pub history (T=2). |
| `regime_map` | 27 | `b`={0.1,0.3,1} × `tau_k`={0.3,1,3} × `k_shape`={1.5,2,3} | 3-D regime map of `signal_fwd` (T=2). |
| `correlation` | 12 | `rho_kr`={−0.5,0,0.5,0.8} × `k_shape`={1.3,2,3.5}, at `tau_k`=0.3 | **Does signal value survive K–R correlation?** (T=2). |
| `pop_size` | 4 | `n`={20,50,100,200} | Finite-n check (T=2). |
| `resource_noise` | 3 | `tau_r`={0.3,1,3} | Resource-signal noise effect (T=2). |
| `tail_map` | 9 | `k_shape`={1.3,2,3.5} × `r_shape`={1.3,2,3.5}, at `tau_k`=0.3 | **Resource-side inequality.** `signal_fwd` by knowledge vs resource tail. ⚠ budget `B∝E[R]` co-varies with `r_shape` — interpret that axis as inequality *and* budget. |
| `info_value` | 12 | `T_rounds`={2,3,4,5} × `tau_k`={0.3,1,3}, at `epsilon`=1e-4 | **Information channel isolated** (ε→0, no compounding): forward's residual edge is pure information value (`fwd_vs_myo_PG`). |
| `horizon_long` | 10 | `T_rounds`={5,6,7,8,10} × `epsilon`={0.3,0.85} | Does the forward advantage saturate past T=5? **Uses `n_steps`=400** (not the base 50 — too coarse beyond T=5; see below). |

**Sweeps that vary T:** `horizon_growth`, `horizon_noise`, `horizon_scale` (T=1..5), `info_value`
(T=2..5), and `horizon_long` (T=5..10). All others are at **T=2**.

> **⚠ Granularity note for `horizon_long`.** The greedy step is `δ = B/n_steps`; at fixed
> `n_steps=50` the per-round granularity thins as T grows, and beyond T=5 it inflates the
> forward advantage severalfold (an artifact). `horizon_long` is therefore run at **`n_steps=400`**
> (converged; verified in `tests/horizon_long_convergence.R`). The T=1..5 sweeps are unaffected —
> `n_steps=50` is granularity-stable there. If you extend any horizon sweep past T=5, raise `n_steps`.

---

## 4. Base parameter defaults (applied wherever a sweep doesn't vary the axis)

From `SWEEP_BASE_PARAMS_T` in `code/sweep_T.R`:

```
T_rounds = 2      n = 50           M = 400          n_steps = 50
epsilon  = 0.1    b = 0.5          x_seed = 0.25    n_pre_rounds = 0
k_shape  = 2      r_shape = 2      k_min = 1        r_min = 1
tau_k    = 1      tau_r = 1        rho_kr = 0       gamma = 1
use_resource_signal = TRUE
```

Budget: `E[R] = r_min·r_shape/(r_shape−1)` (Pareto mean = 2 at defaults); `B = b·n·E[R]`;
**total budget `B_total = 2·B`, split evenly across the T rounds** for non-forward strategies;
forward strategies plan the whole `B_total` freely across rounds.

> **Watch-outs when analyzing:**
> - The T=2 cross-sections run at **`epsilon = 0.1`** (base), *not* 0.3. Growth is only swept in
>   `horizon_growth` and `alpha_regime`.
> - `k_shape` **smaller = heavier tail = more inequality** (counter-intuitive direction).
> - `tau_k` **smaller = sharper signal.**
> - Output (`out_Sj`) is **expected** output (`Σλ`), not realized pub counts — deterministic
>   given the drawn state, hence low variance. Poisson pub draws only drive the planners'
>   learning, not the reported metric.
