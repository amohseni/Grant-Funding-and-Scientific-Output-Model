# State of Play — T-Round Extension

*Housekeeping snapshot. Reflects everything built and validated to date. Read this before
trusting any numbers you may have seen in earlier drafts — one substantive bug was found and
fixed late, and it changed the horizon results.*

---

## 1. Objective

Extend the 2-round Bayesian grant-allocation model (`app.R`, "v5") to **T ∈ {1..5} rounds** and
run the final parameter sweep for the paper. Prime directive: **generalize the existing planner,
don't reinvent** — the T=2 path must reproduce v5 exactly.

## 2. What was built (all complete)

| Deliverable | File(s) | Status |
|-------------|---------|--------|
| T-round simulator + T-round forward/myopic/seed planners | `code/simulate_T.R` | ✅ done |
| Manifest runner + output schema + checkpointing | `code/sweep_T.R` | ✅ done |
| Full parameter sweep (13 sweeps, 172 cells, 200 trials) | `data/` | ✅ done (corrected) |
| Validation suite (anchor, assertions, benchmark) | `validation/` | ✅ all pass |
| Visual results report | `report/results_report.html` | ✅ done |

The **§3 model** (production, compounding, non-persistent grants, 9 strategies) matches v5's
actual equations. The T-round forward planner is a **receding-horizon certainty-equivalent (CE)**
generalization of v5's 2-round planner: at each round it plans the full remaining horizon
(anticipating compounding + one-step information value via `ce_reweight`), executes the current
round, and re-plans next round with the newly observed pubs. At T=2 it reduces **bit-identically**
to v5.

## 3. Validation status (all green)

- **T=2 anchor:** `validation/test_T2_reduction.R` — T-round model at T=2 reproduces v5 to
  floating-point (max |Δ| ≈ 4×10⁻¹⁴, exact grants) across 7 param points × 5 seeds × 9 strategies.
- **§6 assertions:** `validation/assertions_log.txt` — all pass. Covers: T=1 triviality, ε→0
  limit, τ_K limits, heavy-tail finiteness + ESS monitoring, budget conservation, M-convergence,
  CE-vs-scenario valuation gap (−0.26% at T=2, near-exact), sign-path.
- **Granularity stability:** after the fix (§5), the forward advantage is **stable as the greedy
  step `n_steps` is refined** (was the tell that exposed the bug). See `diagnostics/`.
- **Benchmark:** `validation/benchmark_report.txt` — per-run cost O(T²); full manifest ~30 min
  on 7 cores at M=400.

## 4. Current results (the corrected story)

See README §2 for the findings and `report/` for figures. In brief:
**Forward beats myopic everywhere; the advantage grows with knowledge-growth ε (via front-loading
to compound) and, at high ε, with the horizon. Grant-signal value is set by inequality × precision
and survives K–R correlation. Seeding never helps.** All T=2 cross-sections are solid.

## 5. ⚠️ The bug that was found and fixed (important context)

**Symptom.** The first manifest run (M=200, `n_steps=50`) showed the forward-vs-myopic gap
`fwd_vs_myo_PG` *declining* and going *negative* at high T for low ε — a clean-looking
"phase diagram." Follow-up diagnostics revealed this was **not real**: the high-T results flipped
sign erratically with population size n, and *diverged* (to −7) as the greedy granularity `n_steps`
was refined. A correct effect converges; divergence signals a bug.

**Root cause.** `ce_reweight_posterior` (in `app.R`) did a **stochastic SIR resampling** step when
the reweighted effective sample size dropped below M/4. The forward planner's marginals are finite
differences that call this function twice (at `g₁` and `g₁+dg`); when SIR fired, the two calls drew
**different random atoms**, so the marginal was dominated by O(1) resampling noise instead of the
O(dg) signal. Negligible at 2 rounds / coarse steps (why v5 "worked"), catastrophic at long
horizons / fine steps. It was **grant-signal-specific** because the grant signal sharpens the
posterior → lower ESS → SIR fires (pubs-only forward, S7, was always stable).

**Proof it was a bug, not a real CE effect.** At fine granularity the forward greedy scored
*below the myopic schedule on its own CE objective* — impossible for a correct optimizer, since the
myopic schedule is a feasible plan. (`diagnostics/` has this test.)

**Fix.** Made `ce_reweight_posterior` a **deterministic reweight** (no resampling) — atoms fixed,
only importance weights change — so the finite-difference marginals are exact. Bumped **M 200→400**
for effective-sample headroom at heavy tails. After the fix: PG converges, is granularity-stable,
and T=2 is still bit-identical to the (fixed) v5.

**What changed in `app.R`.** Only `ce_reweight_posterior` (≈ lines 278–314) — the SIR branch was
removed and the header comment updated to explain why. This is a **shared primitive**, so it also
affects the **live v5 Shiny app's** forward strategies (S7–S9) — see loose end (a).

**A side benefit:** the fix revealed a clean decomposition. At ε→0 (no compounding) the forward
schedule is uniform yet forward still beats myopic by ~0.2 — pure **information value**, previously
masked by the SIR noise. So: *forward value = compounding value (∝ε) + information value*.

## 6. Loose ends & known caveats

- **(a) The `app.R` fix changes your deployed Shiny app.** The forward planner (S7–S9) now uses the
  deterministic reweight; its outputs differ slightly (for the better) from the previously deployed
  version. **Redeploy the Shiny app when convenient** so the interactive tool matches these results.
- **(b) Fixed params aren't stored per data row** (only the varied axes are). Context is documented
  in `DATA_DICTIONARY.md` §3–4; keep it handy when comparing across sweeps. (Minor; a re-run could
  embed them but isn't worth it.)
- **(c) The ε=0.05 row has an isolated `fwd_vs_myo_PG` spike at T=3** (+0.23, z≈23) that returns to
  ~0 at T=4–5. Likely a mild feature of the "one-step info anticipation, CE-beyond" approximation at
  T=3 (exactly one anticipated round + one CE round). Small; worth a look if T=3 is load-bearing.
- **(d) Heavy-tail cells (α_K→1):** importance-sampling ESS/M < 0.1 for ~¼ of researchers at
  α_K≈1.1 (flagged by assertion #6). M=400 mitigates; read the heaviest-tail `signal_value` /
  `regime_map` cells with mild caution.
- **(e) CE approximation not bounded at high T by a full stochastic-programming reference.** The
  CE-vs-scenario gap is −0.26% at T=2 (near-exact) and post-fix results are granularity-stable, so
  CE is well-behaved — but a full SP-lite planner comparison at T≥3 (spec §6.9 extended) was
  **not built** (deliberately deferred). Optional rigor for the paper.
- **(f) b-axes were capped at [0.1, 1]** per instruction: `funder_scale`, `seed_value`, `regime_map`
  dropped the spec's b-values above 1. Note if a reviewer expects the wider range.
- **(g) Superseded buggy data** lives in the project root at `sweep_results/T_run/` (M=200, pre-fix).
  **Not** in this package. Do not analyze it — use `data/` here (from `sweep_results/T_run_fixed/`).

## 7. Suggested next steps

*Analysis (the reason for this package):*
1. Nail the **information-value vs compounding-value decomposition** — the ε→0 intercept vs the
   ε-slope of `fwd_vs_myo_PG`; separates the two channels cleanly for the paper.
2. **Schedule analysis from `_rawlong`** — quantify front-loading (`b_idx`, `alpha_t` profiles)
   vs ε and T; this is the mechanism behind the headline.
3. Firm up the **signal-value law** (`signal_value`, `regime_map`, `signal_precision`): fit
   `signal_fwd` as a function of tail heaviness and signal precision.

*Additional parameter sweeps* — see **`NEXT_SIMULATIONS.md`** for the full proposal (grids +
rationale, ready to paste into `SWEEP_CONFIGS_T`). Highlights: `tail_map` (resource-side
inequality — the one real gap), `info_value` (isolates the information channel, backs the
compounding-vs-information decomposition), `horizon_long` (does the advantage saturate?).
Recommended minimum ≈ 31 cells / 5–10 min.

*Optional model work:*
4. Investigate caveat (c) (the T=3 spike) if it matters.
5. Build the SP-lite planner (caveat e) to bound the CE tax at high T rigorously.
6. If wider budgets are wanted, re-run `funder_scale`/`regime_map` with b > 1 (caveat f).

## 8. How to reproduce (from the project root, not this folder)

```r
# from the project root:
source_app <- function() { src<-readLines("app.R"); eval(parse(text=paste(
  src[1:(grep("^shinyApp\\(",src)[1]-1L)],collapse="\n")), envir=globalenv()) }
source_app(); source("simulate_T.R"); source("sweep.R"); source("sweep_T.R")
main_sweep_T(seeds = 1:200, out_dir = "sweep_results/T_run_fixed")   # ~30 min, 7 cores
```
`validation/launch_manifest_fixed.R` is the exact driver used for the production run.
