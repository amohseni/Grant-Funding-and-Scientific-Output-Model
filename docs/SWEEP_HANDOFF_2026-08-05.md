# Sweep handoff: seed diagnostics, CES robustness, back-loading attribution

*Written 2026-08-05 by the paper session. Self-contained instructions for a fresh session working in
`Grant-Funding-and-Scientific-Output-Model/`. Three packages, independent, runnable in any order.
Package C is the most invasive code change; do it last. Read `docs/ROBUSTNESS_HANDOFF.md` first for
model orientation; read this document top to bottom before writing any code.*

*Design note: these are hypothesis-testing sweeps on small grids with preregistered predictions, not
open exploration, so grid designs are correct rather than adaptive sampling. Two targeted refinement
rules appear where a boundary is the deliverable (B-Tier-B and C).*

---

## What is being tested, in one screen

**A. Does spreading money thinly actually cost anything?** The model currently says a guaranteed
minimum grant for everyone costs almost nothing (under 0.4% of output). But that floor was only ever
imposed on the weakest optimizer (no peer-review signal), only in round 1, with the rest of the
budget re-optimized around it. Test: impose the floor on the strongest optimizer (S8, with the
signal), in the regime where its targeting is worth most (concentrated talent, sharp signal); and
separately make the floor permanent and total, which must reproduce pure uniform funding exactly, a
built-in correctness check. Decides whether the paper says "floors are nearly free" or "floors cost
exactly what targeting is worth."

**B. Does anything depend on the harmonic production function, and is optimal concentration
monotone in substitutability?** Tier A re-runs the four headline sweeps under other members of the
CES family: directions should survive, sizes may move. Tier B, the gate: one static round, best
static funder, fine grid over the substitution parameter, measure how concentrated the optimal
allocation is (Gini, coverage). Prediction: harder substitution, more concentration, monotonically.
Monotone promotes the concentration-dispute framing to a candidate headline; any interior reversal
kills it as stated.

**C. Why does the optimal funder back-load?** Two candidate drivers: paid growth (grants raise
talent, which argues for spending early to compound) and free growth (talent compounds from baseline
resources anyway, which argues for waiting until researchers are more productive per dollar). The
current model controls both with one knob (ε), so no existing sweep can tell them apart. Test: split
ε into ε_free and ε_paid and turn each off separately. If back-loading survives on pure free growth
and flips to front-loading on pure paid growth, the report's stated mechanism is confirmed; if pure
paid growth also back-loads, the mechanism story is wrong and the section must be rewritten.

---

## 0. Global conventions (all packages)

1. **Git hygiene first.** `git status` will show uncommitted edits to `model.R` (a `budget_ref`
   switch) and `sweep_T.R` (a `resource_regime` sweep config), plus `sweep_results/_probe/`. These
   are intentional. Commit them as they are in one commit before starting, so every run in this
   handoff is attributable to a hash. Record the hash in every output directory as `RUN_INFO.txt`
   (hash, date, seeds, R `sessionInfo()`).
2. **Trials and seeds.** 200 trials/cell (seeds 1:200) unless a package says otherwise. Common
   random numbers across strategies within a run are already guaranteed by the engine (one
   population draw per seed); do not change that.
3. **Contrasts are paired.** All reported contrasts are within-seed differences, mean ± SE over
   seeds, |z| > 2 convention, as in `RESULTS.md`.
4. **Normalization.** Report every effect both raw (Σλ units) and as percent of the same cell's
   no-funding output S1, matching the report's convention.
5. **Allocator.** Harmonic runs: `allocator = "smooth"`. Any non-harmonic run (Package B): the
   smooth myopic water-fill (`waterfill_round_t`) hard-codes the harmonic marginal
   γK²/(K+R+g)² and is silently wrong for any other production form. Use `allocator = "greedy"`
   with `n_steps = 400` (800 for Leontief), and add this guard to `run_simulation_T`:
   `if (allocator == "smooth" && ces_gamma != -1) stop("smooth waterfill is harmonic-only")`.
6. **Validation gates after ANY edit to model.R.** (a) `Rscript tests/test_T2_reduction.R` passes;
   (b) one pinned run (`seed=1, T_rounds=5, epsilon=0.3, strategies=1:9, allocator="smooth"`)
   reproduces its pre-edit `total_expected` values bit-identically at default arguments (capture the
   pre-edit values before touching anything); (c) `tests/assertions_T.R` passes. Do not run any
   sweep until all three are green.
7. **Outputs.** One directory per sweep under `sweep_results/` (names below), with the standard
   `_raw.rds / _rawlong.rds / _summary.rds` triple via `run_sweep_T`/`summarize_sweep` machinery,
   plus a summary CSV export. Finish the whole job with a results memo
   `docs/DIAGNOSTICS_RESULTS_2026-08.md` that answers each preregistered prediction below with
   CONFIRMED / REFUTED / MIXED and the numbers. Do not editorialize beyond that; the paper session
   interprets.
8. **Where to run.** This job wants Claude Code on Aydin's Mac, in the repo directory: R and the
   data live there, and the cloud container has no R installed. Estimate cores with
   `parallel::detectCores()` before trusting the runtime numbers below.
9. **Runtime calibration.** The 16-sweep manifest (219 cells × 200 trials) runs ~30-35 min on ~8
   cores at M=400. Use that to sanity-check the estimates below; if any package is tracking >3x its
   estimate, stop and reassess rather than letting it run all night silently.

---

## Package A: seed-floor diagnostics (D3, D4, D1)

**Question.** Is the near-zero cost of a uniform seed floor a general property, or an artifact of
only ever having diluted weak (no-signal) optimizers with a round-1-only floor? Full context:
`for Claude/drafts/grant-funding/notes-seed-floor.md`.

**Preregistered predictions (paper session's, falsifiable):**
- P-A1: seed cost scales with the diluted family's targeting value times the floored budget share,
  with proportionality κ ∈ (0.3, 1]. Concretely at (α_K=1.3, τ_K=0.3, b=0.5, x_seed=0.75, T=2),
  cost of the floor on S8 is 2-8% of S1, not 0.4%.
- P-A2: with a persistent every-round floor at x_seed=1 in the no-signal myopic family, the cost
  equals S4−S2 for the same seeds (exactly: that configuration IS S2, so the contrast is an identity
  check on the implementation; treat any nonzero difference as a bug in the new code).
- P-A3: cost is superlinear (convex) in x_seed at fixed b.
- Alternative outcome that would refute P-A1: signal-family seed cost stays under 0.5% of S1
  everywhere. Then "floors are nearly free" is robust and the paper claims it.

### A1. Code change: seed variants of the signal strategies (D3)

Add strategies S10 and S11, exact mirrors of S6 and S9 but with the grant signal:
- S10 = myopic Bayes (pubs + grant signal + seed): round 1 allocates fraction `x_seed` of the
  round-1 tranche uniformly, remainder by myopic optimization with pubs + σ_R + σ_K; rounds 2..T
  full myopic optimization, no seed. (S10 is to S5 what S6 is to S4.)
- S11 = forward Bayes (pubs + grant signal + seed): same floor, remainder by forward optimization
  with the grant signal. (S11 is to S8 what S9 is to S7.)

Implementation: follow exactly how S6/S9 differ from S4/S7 in the strategy dispatch in `model.R`
(strategy loop inside `run_simulation_T`; grep `x_seed` for the existing seed logic) and how
S5/S8 differ from S4/S7 (`use_grant_signal`). Extend `extract_metrics_T` in `sweep_T.R` with
`out_S10`, `out_S11` and contrasts `seed_sig_myo = S10 − S5`, `seed_sig_fwd = S11 − S8`.
Also add, if the allocation vectors are accessible at metric-extraction time cheaply, the displaced
budget share `disp_S11 = sum(abs(g_S11 − g_S8))/(2B)` per round-1 (this is diagnostic D2, the
direct measure of how much of the floor the re-optimizer absorbs; skip if it requires `detail=TRUE`
plumbing, and note the skip in the memo).

### A2. Sweep D3: `sweep_results/D3_seed_signal/`

Grid (T=2, smooth allocator, strategies `c(1,2,4,5,7,8,10,11)`):

| axis | values |
|---|---|
| b | 0.1, 0.3, 0.5, 1.0 |
| x_seed | 0.1, 0.25, 0.5, 0.75 |
| (k_shape, tau_k) | (2, 1) base; (1.3, 0.3) heavy+sharp |

32 cells × 200 seeds. Estimated runtime: ~5-10 min. Primary readouts per cell: `seed_sig_myo`,
`seed_sig_fwd`, both raw and as % of S1, next to the family targeting value (S5−S2 for the same
cell) so the κ ratio in P-A1 is computable directly.

### A3. Code change + sweep D4: persistent floor, `sweep_results/D4_seed_persistent/`

Add a run option `seed_persistent = FALSE`; when TRUE, seed strategies apply the x_seed floor in
EVERY round's tranche, not round 1 only. Touches only the seed branch of the strategy logic.

Grid (T=2, smooth, strategies `c(1,2,4,6)` with `seed_persistent=TRUE`):

| axis | values |
|---|---|
| x_seed | 0.25, 0.5, 0.75, 0.9, 1.0 |
| (k_shape, tau_k) | (2, 1); (1.3, 0.3) |

10 cells × 200 seeds, minutes. Readout: `seed_myo = S6 − S4` as a function of x_seed, and the
x_seed=1.0 identity check against S2 (P-A2). Plot cost vs x_seed with the linear-interpolation
reference line from 0 to −(S4−S2); the curve's position relative to that line is the envelope
result (P-A3).

### A4. Analysis-only (D1, no runs)

From the existing `sweep_results/T_run_smooth/seed_value_summary.rds`: fit cost vs x_seed at each b;
report whether curvature is convex (P-A3 on the old round-1 design). Include in the memo.

### A5. Analysis-only: the T-scaling exponent (no runs)

The report claims planning value grows "quadratically in T"; that needs a fitted exponent, not an
impression. From `T_run_smooth/horizon_growth_summary.rds` and `horizon_long_summary.rds`: fit
log(PG) on log(T) at each ε with PG > 0 (ε ≥ 0.3), report the exponent and its SE per ε, pooled and
separately for T ≤ 5 and T = 5..10. Include in the memo; the paper will cite this number.

---

## Package B: CES family (σ robustness + the concentration gate)

**Question, two tiers.** Tier A: do the paper's headline results survive across the admissible CES
family (γ_ces ≤ 0)? Tier B: is the concentration of optimal funding monotone in the elasticity of
substitution σ? Tier B gates a candidate paper framing (the concentration dispute as a disagreement
about σ); it is a result in its own right, not a robustness appendix.

**Preregistered predictions:**
- P-B1 (Tier A): the direction of every headline comparative static survives at γ_ces ∈ {0, −3}:
  signal value increasing in talent inequality and precision; seed cost small in the no-signal
  family; forward advantage ≈ 0 at ε ≈ 0.1, T=2. Effect SIZES may move; directions should not.
- P-B2 (Tier B): Gini of the optimal allocation is monotone DECREASING in σ (equivalently,
  concentration increases as γ_ces falls toward Leontief) at every (α_K, b) cell. Any interior
  non-monotonicity refutes P-B2 and kills the Story-5 framing as stated; report its location, do
  not smooth over it.
- P-B3: coverage (fraction funded) is monotone increasing in σ.

### B1. Code change: the CES switch

`model.R` currently hard-codes the harmonic form in exactly two primitives (`lambda_rate` at ~:65,
`update_knowledge` at ~:69); all 28+17 call sites route through them. Replace with:

```r
ces_mean <- function(K, R, gc) {
  if (is.infinite(gc) && gc < 0) return(pmin(K, R))          # Leontief
  if (abs(gc) < 1e-8) return(sqrt(K * R))                    # Cobb-Douglas limit
  ((K^gc + R^gc) / 2)^(1 / gc)
}
lambda_rate      <- function(K, R, gamma)   gamma   * ces_mean(K, R, CES_GAMMA) / 2
update_knowledge <- function(K, R, epsilon) K + epsilon * ces_mean(K, R, CES_GAMMA) / 2
```

with `CES_GAMMA` set at entry to `run_simulation_T(..., ces_gamma = -1)` (assignment at function
entry is fork-safe under `mclapply`; note in a comment that threading it as an argument through all
call sites is the cleaner refactor if anyone has the appetite). **Conventions, do not deviate:**

- The **/2 keeps γ_ces = −1 bit-identical to the current model**: ces_mean(K,R,−1)/2 = KR/(K+R).
  Validation gate 6(b) must pass with defaults after this edit. (The model write-up's factor
  productivity A relates to the code's `gamma` by A = gamma/2; that is a prose-side matter, not a
  code-side one.)
- **Naming: the code's `gamma` is factor productivity; the CES exponent is `ces_gamma`/`CES_GAMMA`
  everywhere.** Never reuse the bare name for the exponent; the report and model write-up already
  collide on this symbol and the code must not join them.
- Same γ_ces governs production and knowledge growth (coupled, the write-up's default). No
  decoupled (γ_λ, γ_K) runs in this package.
- The posterior/likelihood (`loglik_pubs`) calls `lambda_rate`, so the funder's inference model
  tracks the true form automatically. That is intended (the funder knows the technology).

### B2. Tier A: `sweep_results/sigma_tierA_<sweep>/`

Re-run four sweeps from the manifest at γ_ces ∈ {0, −3}, greedy n_steps=400, **100 seeds** (1:100;
directions need less precision than the harmonic headline numbers, and greedy-400 is slow):

1. `signal_value` (α_K × τ_K grid, T=2) — backs the inequality × precision law.
2. `horizon_growth` (T × ε) — backs forward-tracks-compounding and b_idx back-loading.
3. `seed_value` (b × x_seed) — backs the floor result.
4. `correlation` (ρ × α_K) — backs correlation robustness.

Leontief (γ_ces = −Inf) is a boundary check, not an interior member: run `signal_value` only, at
n_steps=800, 50 seeds, and report with the caveat that allocation at the kink is degenerate
(marginal-return equalization has no bite; ties are broken by fill order). Skip the other three at
Leontief.

Estimated runtime: this is the slow package. Greedy-400 costs roughly 8x smooth per run.
Ballpark: several hours on 8 cores. Run overnight; per convention 8, checkpoint per sweep (the
`sweep_one_T` resume logic already skips completed summaries).

### B3. Tier B: the concentration gate, `sweep_results/sigma_tierB_concentration/`

The cleanest object is static: **T=1, strategy S5 only** (myopic + grant signal; at T=1 myopic is
optimal and no dynamics confound Force A). Add metrics to `extract_metrics_T` if absent:
`gini_g_S5` (Gini of the round-1 grant vector), `coverage_S5` (share of researchers with
g > 1e-8), `top10_share_S5` (share of budget to the top decile of grants).

Grid (T=1, greedy n_steps=400, 200 seeds; T=1 is cheap):

| axis | values |
|---|---|
| ces_gamma | 0, −0.25, −0.5, −1, −2, −3, −6, −12 (σ = 1, .8, .67, .5, .33, .25, .143, .077) |
| k_shape | 1.3, 2, 3.5 |
| b | 0.1, 0.5, 1.0 |

72 cells × 200 = 14,400 T=1 runs; ~15-30 min. Leontief endpoint separately at n_steps=800, reported
with the degeneracy caveat.

**Refinement rule (the one adaptive element):** if Gini vs γ_ces is non-monotone in any (k_shape, b)
cell beyond noise (|z| > 2 on the offending pairwise difference), insert two additional γ_ces points
bracketing the reversal and re-run those cells at 400 seeds before reporting. The boundary location
is the deliverable if P-B2 fails.

---

## Package C: back-loading attribution (decoupled free vs paid growth)

**Question.** The report attributes back-loading to free force C (knowledge compounds from baseline
regardless of grants, so later rounds face higher-K researchers and a dollar does more there)
dominating paid force B (grants raise K, which argues for spending early to compound). Both forces
scale with ε in the current model, so the existing sweeps cannot discriminate them. Decouple them.

**Preregistered predictions:**
- P-C1 (pure free: ε_paid = 0, ε_free > 0): the forward planner still back-loads (b_idx_S8 > 0.5),
  more so as ε_free rises. This is the report's stated mechanism in isolation.
- P-C2 (pure paid: ε_free = 0, ε_paid > 0): the planner FRONT-loads (b_idx_S8 < 0.5): with growth
  only through grants, early spending compounds over more remaining rounds and there is no free
  level-shift to wait for.
- P-C3: the diagonal ε_free = ε_paid reproduces the current model's b_idx and PG bit-identically
  (validation identity, not a prediction).
- If P-C2 fails (planner back-loads even with pure paid growth), the report's mechanism story is
  wrong as stated and the saturation-timing interaction needs rethinking; report the schedule
  shapes (rawlong per-round shares), do not interpret.

### C1. Code change (invasive; do last; budget a careful hour)

Current growth entangles the two channels: K_{t+1} = K_t + ε·Λ(K, R0+g)/2, where the increment
divides into the part that happens without the grant and the grant-caused part. Replace with:

```
K_{t+1} = K_t + eps_free * lam_base + eps_paid * (lam_post - lam_base)
  where lam_base = ces_mean(K, R0,     CES_GAMMA)/2     # growth from baseline resources
        lam_post = ces_mean(K, R0 + g, CES_GAMMA)/2     # post-grant
```

`run_simulation_T` gains `eps_free = epsilon, eps_paid = epsilon` (defaults reproduce the current
model exactly: eps_free·lam_base + eps_paid·(lam_post − lam_base) = ε·lam_post).

**The hard part:** `update_knowledge` has ~17 call sites and the decomposition needs the baseline
R0 and the grant g separately at each. Classify every call site before editing: (a) realized
dynamics in the main simulation loop; (b) the forward planner's internal anticipation
(`fwd_researcher_value` / `plan_forward_ce` and the smooth forward path); (c) any CE-reweighting
path. The funder knows the model, so (b) must use the same decomposition as (a); a planner
anticipating the old dynamics while the world runs the new ones is a bug, not an experiment.
After the edit, validation gate 6 plus P-C3 on three seeds at T=5 before any sweep.

### C2. Sweep: `sweep_results/bload_decouple/`

Grid (T=5, smooth allocator, harmonic, strategies `c(1,5,8)`, 200 seeds):

| axis | values |
|---|---|
| eps_free | 0, 0.1, 0.3, 0.85 |
| eps_paid | 0, 0.1, 0.3, 0.85 |

16 cells × 200 seeds at T=5: ~20-30 min. Readouts: `b_idx_S8`, `fwd_vs_myo_PG` (S8−S5), and the
full per-round budget shares from rawlong (the schedule shape is the evidence, not just its
center-of-mass; keep rawlong).

**Refinement rule:** if the b_idx = 0.5 contour falls inside the grid (P-C2 says it does, along the
eps_free axis near zero), add a 5-point transect perpendicular to the contour at the two most
informative crossings (200 seeds each) so the boundary is located to within a factor of ~2 in
eps_free/eps_paid ratio. The contour is the paper's figure if the attribution section gets one.

### C3. Optional corroboration (no code change, 5 min)

The already-configured `resource_regime` sweep (uncommitted, in `sweep_T.R`) maps b_idx over
r_min × ε with a decoupled purse (`budget_ref="K"`). Run it at 200 seeds if time permits; it probes
the same mechanism from the resource side (the paper session's F9) and its config is already
written. Output to `sweep_results/resource_regime/`.

---

## Flags to carry back to the paper session (not this session's job to fix)

1. **Model write-up Eq (10) contradicts the code and the report.** Eq (10) states a per-round
   budget constraint with "budget is not transferable across rounds"; the implemented forward
   strategies allocate the total remaining budget across the horizon (report p.8: round-1 share
   0.08 to 0.37 at T=5). The write-up's constraint describes the non-forward strategies only.
   The spec must be corrected before the model section is drafted.
2. Default-parameter tables disagree between the model write-up (T=10, n=100, B=50) and the code
   base params (T=2, n=50, b=0.5). Reconcile in the spec.
3. Whatever this session finds, `RESULTS.md` at the repo root is stale relative to the exact-
   allocator story and should be regenerated or clearly deprecated after these runs land.

## Completion checklist

- [ ] Pre-edit pinned-run values captured; all validation gates green after each code change.
- [ ] A: D3 (32 cells), D4 (10 cells), D1 analysis; κ ratios computed; identity check P-A2 exact.
- [ ] B: Tier A four sweeps × {0, −3} + Leontief signal_value; Tier B 72 cells + any refinement.
- [ ] C: decoupling edit validated (P-C3), 16 cells + contour transects; resource_regime optional.
- [ ] Every output dir has RUN_INFO.txt; memo `docs/DIAGNOSTICS_RESULTS_2026-08.md` answers
      P-A1..P-C3 with CONFIRMED / REFUTED / MIXED and numbers; everything committed and pushed.
