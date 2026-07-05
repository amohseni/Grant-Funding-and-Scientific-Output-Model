# Next Simulations to Run

Proposed additional parameter sweeps for the T-round model, beyond the 13 in `data/`. Each is
motivated by a specific analysis question, sized to slot into `SWEEP_CONFIGS_T` (`code/sweep_T.R`)
as a new entry, and cheap to run (all ~200 seeds; the full set below is ~90 cells ≈ 15–20 min on
7 cores). Symbol↔code names and base defaults are in `DATA_DICTIONARY.md`.

**How to run these:** add each `grid_fn` block below to `SWEEP_CONFIGS_T`, then
`sweep_one_T("<name>", seeds = 1:200, out_dir = "sweep_results/T_run_supplement")`. Same output
schema as the existing data, so they drop straight into the same analysis.

**Priority:** ⭐⭐⭐ = do this; ⭐⭐ = strong value; ⭐ = nice-to-have / robustness.

---

## A. Fills a genuine gap

### A1 · `tail_map` — resource-side inequality ⭐⭐⭐
**Question.** The whole "signal value depends on inequality" story varies only the **knowledge**
tail (`k_shape`, α_K); the **resource** tail (`r_shape`, α_R) is fixed at 2 in *every* existing
sweep. Since production `KR/(K+R)` is roughly symmetric in K and R, resource inequality may drive
results too. Does signal value track knowledge inequality *specifically*, or *any* inequality?

**Grid.** `k_shape` {1.3, 2, 3.5} × `r_shape` {1.3, 2, 3.5}, at `tau_k`=0.3, T=2. **9 cells.**
Watch: `signal_fwd`, `out_S*`, and how the K-tail vs R-tail effects compare.

```r
tail_map = list(
  name = "tail_map", tier = 2,
  description = "Knowledge tail k_shape vs resource tail r_shape at T=2, sharp signal (tau_k=0.3).",
  grid_fn = function() cbind(expand.grid(k_shape = c(1.3, 2, 3.5), r_shape = c(1.3, 2, 3.5)), tau_k = 0.3),
  varied_params = c("k_shape", "r_shape"),
  primary_plot = list(type = "heatmap", x_var = "r_shape", y_var = "k_shape",
                      fill_var = "signal_fwd_mean", text_var = "signal_fwd_mean",
                      title = "Signal value vs knowledge and resource tails", fill_label = "S8-S7"),
  secondary_plot = NULL)
```

> **⚠ Budget confound to document when analyzing.** The budget uses `E[R] = r_min·α_R/(α_R−1)`, so
> a heavier resource tail (smaller `r_shape`) *raises* E[R] and therefore the absolute budget
> `B = b·n·E[R]`. The `r_shape` axis thus mixes **resource inequality** with **budget scale**
> (`k_shape` has no such confound). Two clean options:
> 1. Report it as-is and note the confound (heavier R tail ⇒ more inequality *and* more money).
> 2. Add a companion sweep holding `B_total` fixed: set `b` per cell to `0.5 · 2 / E[R](r_shape)`
>    (= 0.231 at r_shape=1.3; 0.5 at r_shape=2; 0.714 at r_shape=3.5). Isolates inequality from budget.

---

## B. Sharpens the headline findings

### B1 · `info_value` — isolate the information channel ⭐⭐⭐
**Question.** The key new result is the decomposition *forward value = compounding (∝ε) +
information value (survives ε→0)*. It currently rests on a single ε→0 assertion point. At **ε=0**
there is no compounding, so **all** of forward's advantage is information value. Map that channel
across horizon and signal precision.

**Grid.** `epsilon` = 1e-4 (held ≈0) with `T_rounds` {2,3,4,5} × `tau_k` {0.3, 1, 3}. **12 cells.**
Watch: `fwd_vs_myo_PG` (pure info value here) and `signal_fwd`. Expect small-positive PG that
depends on τ_K (info only matters if the signal is readable) and grows modestly with T.

```r
info_value = list(
  name = "info_value", tier = 1,
  description = "Information channel isolated: eps~0 (no compounding) across horizon T and signal noise tau_k.",
  grid_fn = function() cbind(expand.grid(T_rounds = c(2,3,4,5), tau_k = c(0.3, 1, 3)), epsilon = 1e-4),
  varied_params = c("T_rounds", "tau_k"),
  primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                      color_var = "tau_k", title = "Information value (eps->0) across horizon",
                      y_label = "S8 - S5 (pure info value)"),
  secondary_plot = NULL)
```

### B2 · `horizon_long` — does foresight's advantage saturate? ⭐⭐
**Question.** The corrected `fwd_vs_myo_PG` is monotone-increasing in T at high ε. Does it saturate
or keep climbing beyond T=5?

**Grid.** `T_rounds` {5, 6, 7, 8, 10} × `epsilon` {0.3, 0.85}. **10 cells.** Cost grows as O(T²) but
T=10 is still cheap; results are granularity-stable post-fix. Watch: `fwd_vs_myo_PG`, `b_idx_S8`
(does front-loading keep increasing?).

```r
horizon_long = list(
  name = "horizon_long", tier = 1,
  description = "Long horizon T=5..10 at eps in {0.3,0.85}: does the forward advantage saturate?",
  grid_fn = function() expand.grid(T_rounds = c(5,6,7,8,10), epsilon = c(0.3, 0.85)),
  varied_params = c("T_rounds", "epsilon"),
  primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                      color_var = "epsilon", title = "Forward advantage at long horizons",
                      y_label = "S8 - S5"),
  secondary_plot = NULL)
```

---

## C. Interactions & robustness (lower priority)

### C1 · `signal_complementarity` — do the two signals substitute? ⭐⭐
**Question.** `resource_noise` (τ_R alone, at τ_K=1) showed τ_R doesn't matter. But do the grant
and resource signals **substitute** — does a good resource read make the grant signal less valuable?

**Grid.** `tau_k` {0.3, 1, 3} × `tau_r` {0.3, 1, 3}, T=2. **9 cells.** Watch: `signal_fwd` across the
τ_R axis at each τ_K.

```r
signal_complementarity = list(
  name = "signal_complementarity", tier = 2,
  description = "Grant-signal noise tau_k vs resource-signal noise tau_r at T=2: do the signals substitute?",
  grid_fn = function() expand.grid(tau_k = c(0.3, 1, 3), tau_r = c(0.3, 1, 3)),
  varied_params = c("tau_k", "tau_r"),
  primary_plot = list(type = "heatmap", x_var = "tau_r", y_var = "tau_k",
                      fill_var = "signal_fwd_mean", text_var = "signal_fwd_mean",
                      title = "Signal value vs both signal noises", fill_label = "S8-S7"),
  secondary_plot = NULL)
```

### C2 · `horizon_noise_hieps` / `horizon_scale_hieps` — interactions where the effect is large ⭐
**Question.** `horizon_noise` (T×τ_K) and `horizon_scale` (T×b) both run at the base ε=0.1, where
the forward advantage is smallest. The interesting interactions live at *high* ε. Re-run at ε=0.55.

**Grid.** Each is the existing sweep with `epsilon = 0.55` added. `horizon_noise_hieps`: T {1..5} ×
τ_K {0.1,0.3,1,3,10} (25). `horizon_scale_hieps`: T {1..5} × b {0.1,0.3,0.5,1} (20). **45 cells.**
*(Alternative: a single 3-way `T {1..5} × ε {0.1,0.3,0.55} × τ_K {0.3,1,3}` = 45 cells.)*

```r
horizon_noise_hieps = list(
  name = "horizon_noise_hieps", tier = 1,
  description = "Horizon T vs signal noise tau_k at HIGH growth (eps=0.55), where forward's edge is large.",
  grid_fn = function() cbind(expand.grid(T_rounds = 1:5, tau_k = c(0.1,0.3,1,3,10)), epsilon = 0.55),
  varied_params = c("T_rounds", "tau_k"),
  primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                      color_var = "tau_k", title = "Forward advantage across horizon and noise (eps=0.55)",
                      y_label = "S8 - S5"),
  secondary_plot = NULL)
# horizon_scale_hieps: same pattern with b = c(0.1,0.3,0.5,1), varied c("T_rounds","b").
```

### C3 · `pop_horizon` — are the horizon results n-robust? ⭐
**Question.** The bug that was fixed was n-sensitive. Confirm the corrected horizon results are
robust to population size (and see whether the advantage scales with n). Pure due-diligence.

**Grid.** `n` {30, 50, 100} × `T_rounds` {2, 3, 5}, at ε=0.3. **9 cells.** Watch: `fwd_vs_myo_PG`
stable in sign/shape across n.

```r
pop_horizon = list(
  name = "pop_horizon", tier = 3,
  description = "Population size n vs horizon T at eps=0.3: robustness of the horizon result to n.",
  grid_fn = function() cbind(expand.grid(n = c(30, 50, 100), T_rounds = c(2, 3, 5)), epsilon = 0.3),
  varied_params = c("n", "T_rounds"),
  primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                      color_var = "n", title = "Forward advantage vs horizon across population size",
                      y_label = "S8 - S5"),
  secondary_plot = NULL)
```

---

## D. Optional refinement

If the **signal-value law** (`signal_fwd` vs α_K, τ_K) or any regime boundary turns out to be sharp
rather than smooth, refine the grid adaptively — sample densely where the metric changes fastest
rather than on a uniform grid (see the `adaptive-parameter-sweeps` skill). Only worth it if the
coarse `signal_value` / `regime_map` grids reveal a boundary you want to pin down.

---

## Recommended minimum
Do **A1 `tail_map`** (fills the real gap) and **B1 `info_value`** (backs the cleanest new result);
add **B2 `horizon_long`** for the asymptotic story. That's **31 cells ≈ 5–10 min**. The rest is
robustness you can add if reviewers ask.

## Not worth sweeping
`gamma` (linearly rescales output — uninformative), `k_min`/`r_min` (scale params), `M`/`n_steps`
(numerical, already validated). Varying the **production functional form** is a structural
robustness question — a paragraph in the paper, not a manifest entry.
