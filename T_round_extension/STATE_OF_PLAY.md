# State of Play — T-Round Extension

*Housekeeping snapshot of the completed study. Reflects everything built and validated. Read this
before trusting any numbers from earlier drafts — one substantive bug was found and fixed (§5), and
two diligence checks refined the results (§4a). All data here is post-fix.*

> **Note on scope of this file.** This `T_round_extension/` package is a self-contained analysis
> bundle; its `code/` folder mirrors the current repo-root code (the model engine **`model.R`**,
> slim `app.R`, `sweep.R`, `sweep_T.R` — the former `simulate_T.R` was absorbed into `model.R`; see
> §2). For the current code and the path-to-publication, see the repo root `README.md`, `ROADMAP.md`,
> and `RESULTS.md`.

---

## 1. Objective

Extend the 2-round Bayesian grant-allocation model (v5) to **T ∈ {1..5+} rounds**, run the final
parameter sweeps, and get the study analysis-ready and shareable. Prime directive: **generalize the
existing planner, don't reinvent** — the T=2 path must reproduce v5 exactly (it does, bit-identically).

## 2. What was built (all complete)

| Deliverable | Where | Status |
|-------------|-------|--------|
| T-round model engine (`run_simulation_T` + forward/myopic/seed planners) | `model.R` (root) | ✅ |
| Shiny app generalized to T rounds, sourcing the engine | `app.R` (root) | ✅ |
| Manifest runner + output schema + checkpointing | `sweep_T.R` | ✅ |
| **Full parameter sweep — 16 sweeps, 219 cells, 200 trials/cell** | `data/` | ✅ (corrected) |
| Validation suite (anchor, assertions, benchmark) + diagnostics | `validation/`, `diagnostics/` | ✅ |
| Publication figures (8×, PNG+PDF) | root `figures/` | ✅ |
| Statistical digest for the paper | root `RESULTS.md` | ✅ |
| Reproducibility + packaging (reproduce, license, citation, env, OSF plan) | root | ✅ |

The model (production `λ=γKR/(K+R)`, compounding knowledge, non-persistent grants, 9 strategies)
matches v5's equations. The T-round forward planner is a **receding-horizon certainty-equivalent
(CE)** generalization: each round it plans the full remaining horizon (anticipating compounding +
one-step information value), executes the current round, then re-plans with the newly observed pubs.
At T=2 it reduces **bit-identically** to v5.

**Code refactor (post-manifest):** the model was consolidated from `app.R` + `simulate_T.R` into a
single base-R engine **`model.R`**; `app.R` became a slim Shiny layer (`source("model.R")` + UI) with
a `T` rounds slider; all sweeps/tests now `source("model.R")`. T=2 stays bit-identical throughout.

## 3. Validation status (all green)

- **T=2 anchor:** T-round model at T=2 reproduces v5 to floating point (max |Δ| ≈ 4×10⁻¹⁴, exact
  grants) — `validation/test_T2_reduction.R`.
- **§6 assertions:** all pass (`validation/assertions_log.txt`) — T=1 triviality, ε→0 limit, τ_K
  limits, heavy-tail finiteness + ESS, budget conservation, M-convergence, CE-vs-scenario gap, sign-path.
- **CE approximation quality (now checked at high T):** the CE-vs-scenario valuation gap stays
  **within ±0.4% of forward output across T=2–5 and does not grow with T** (`tests/ce_tax_vs_T.R` on
  the fixed model). So the headline results are not CE artifacts — the deferred full SP-lite planner
  was **not needed**.
- **Granularity stability:** forward's advantage is stable as the greedy step `n_steps` is refined —
  for T≤5 at `n_steps=50`, and for T>5 at `n_steps=400` (see §4a).
- **Benchmark:** per-run cost O(T²); the 16-sweep manifest ~30–35 min on ~8 cores at M=400.

## 4. Current results (the corrected story)

Full digest with effect sizes/z in root `RESULTS.md`; figures in `figures/`. In brief:

**Forward beats myopic everywhere, growing with knowledge-growth ε (via front-loading to compound)
and, at high ε, with the horizon.** Grant-signal value is governed by **inequality × precision**
(≈21 at heavy tail + sharp signal → ~0 at light tail + noisy) and **survives K–R correlation**
(large even at ρ_s=0.78). **Uniform seeding never helps.** The forward edge **decomposes** into a
compounding channel (∝ε) plus a small, T=3-peaked information channel. Robustness nulls: resource-
signal noise irrelevant; no finite-n artifact; budget scale doesn't drive the gain.

### 4a. Two diligence findings that refined the results

- **`horizon_long` granularity artifact (caught & fixed).** At the manifest default `n_steps=50`,
  the greedy's per-round granularity thins as T grows and **inflated the T>5 advantage ~4×** (an
  artifact — it oscillated wildly with T). Re-run at converged `n_steps=400`, the true story is
  clean and **monotone: forward's advantage keeps growing slowly with T, no saturation, no collapse**
  (ε=0.3: +0.32→+0.85 over T=5→10; ε=0.85: +2.0→+3.4). T=1–5 was unaffected (granularity-stable at
  50). Convergence verified in `tests/horizon_long_convergence.R`; `horizon_long` config now pins
  `n_steps=400`.
- **The T=3 information peak, isolated.** The `info_value` sweep (ε≈0, no compounding) shows the
  information channel is **small and sharply peaked at T=3** (~+0.22, ~0 elsewhere) — a feature of
  the planner's one-step information anticipation, now characterized rather than a mystery.
- **`tail_map` (resource-side inequality).** Confirms knowledge inequality dominates signal value,
  with resource inequality a secondary modulator (⚠ budget co-varies with the resource tail).

## 5. ⚠️ The bug that was found and fixed (important context)

**Symptom.** The first manifest run (`n_steps=50`) showed `fwd_vs_myo_PG` declining/going negative at
high T — a clean-looking phase diagram. Diagnostics showed it wasn't real: it flipped sign
erratically with population size n and **diverged** (to −7) as `n_steps` was refined. A correct
effect converges; divergence signals a bug.

**Root cause.** `ce_reweight_posterior` did a **stochastic SIR resampling** step when the reweighted
ESS dropped below M/4. The forward planner's marginals are finite differences that call it twice (at
`g₁` and `g₁+dg`); when SIR fired, the two calls drew **different random atoms**, so the marginal was
dominated by O(1) resampling noise instead of the O(dg) signal — negligible at 2 rounds/coarse steps
(why v5 "worked"), catastrophic at long horizons/fine steps. **Grant-signal-specific** (a sharp K
posterior lowers ESS → SIR fires; pubs-only forward S7 was always stable).

**Proof it was a bug, not a real CE effect.** At fine granularity the forward greedy scored *below
the myopic schedule on its own CE objective* — impossible for a correct optimizer.

**Fix.** Made `ce_reweight_posterior` a **deterministic reweight** (no resampling); bumped **M 200→400**
for ESS headroom. After: PG converges, is granularity-stable, T=2 still bit-identical to fixed v5. A
side benefit — the fix revealed the compounding + information decomposition (§4).

**What changed:** only `ce_reweight_posterior` (now in `model.R`). It's a shared primitive, so it
also fixes the live Shiny app's forward strategies — **redeploy the app** (loose end (a)).

## 6. Loose ends & known caveats

- **(a) Redeploy the Shiny app.** `app.R` now sources `model.R` and uses the fixed planner; the live
  shinyapps.io app still has the pre-fix planner. Bundle must include both files. *(Open — yours.)*
- **(b) Fixed params aren't stored per data row** (only varied axes). Context in `DATA_DICTIONARY.md`
  §3–4.
- **(c) T=3 information peak** — now isolated/characterized via `info_value` (§4a); a real feature of
  one-step anticipation, small; report honestly, don't over-generalize to "information grows with T".
- **(d) Heavy-tail cells (α_K→1):** IS ESS/M < 0.1 for ~¼ of researchers at α_K≈1.1; M=400 mitigates.
  Read the heaviest-tail `signal_value`/`regime_map`/`tail_map` cells with mild caution.
- **(e) `horizon_long` needs `n_steps=400`** (§4a); the manifest default 50 is too coarse past T=5.
  If you extend any horizon result beyond T=5, raise `n_steps`. *(The CE-approximation concern this
  used to flag is now resolved — see §3.)*
- **(f) Budget confound in `tail_map`:** the α_R axis co-varies with the budget (B ∝ E[R]); interpret
  it as inequality *and* budget scale, not pure resource inequality.
- **(g) b-axes capped at [0.1, 1]** in `funder_scale`/`seed_value`/`regime_map`.
- **(h) Superseded buggy data** at `sweep_results/legacy/T_round_buggy_M200/` (pre-fix). **Do not
  analyze it** — use `data/` here (the corrected run).

## 7. What's left (path to publication)

The analysis and packaging are done. Remaining items are **yours** (see root `ROADMAP.md` for the
full plan): fill `CITATION.cff` (ORCID, co-authors); **materialize the iCloud files** before any
upload (`docs/OSF_UPLOAD.md`); commit/push to GitHub; redeploy the Shiny app; create + register the
OSF project; write the paper (use `RESULTS.md` + `figures/`, and the handoff prompt in
`docs/PAPER_HANDOFF.md`).

## 8. How to reproduce (from the project root)

```r
source("model.R"); source("sweep.R"); source("sweep_T.R")
main_sweep_T(seeds = 1:200, out_dir = "sweep_results/T_run_fixed")   # ~30–35 min, ~8 cores
```
Or `Rscript reproduce.R` for the full end-to-end regeneration + validation. Environment/versions in
`docs/ENVIRONMENT.md`. `horizon_long` self-pins `n_steps=400`; all other sweeps use the base
`n_steps=50`.
