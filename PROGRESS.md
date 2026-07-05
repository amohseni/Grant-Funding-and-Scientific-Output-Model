# PROGRESS — T-Round Extension + Final Parameter Sweep

> **📦 Organized deliverable: [`T_round_extension/`](T_round_extension/).** That folder is the
> clean, self-contained analysis package — start at `T_round_extension/README.md`
> (model + findings + how to use), `STATE_OF_PLAY.md` (status, the bug story, loose ends,
> next steps), and `DATA_DICTIONARY.md` (data schema). This PROGRESS.md is the raw running
> log kept for provenance.

Spec: T-round generalization of the v5 Bayesian grant-allocation model + final
manifest sweep. Prime directive: generalize `app.R` (v5), don't reinvent; the
T=2 path reproduces v5 exactly (algorithmically).

## Approved decisions (from spec review)
- **Budget normalization:** KEEP v5 **mean-E[R]**, default **b=0.5**
  (`B_total = 2B = 2·b·n·E[R]`, `E[R]=r_min·α_R/(α_R−1)`). [UPDATED — the
  earlier median/b=0.1 choice was reverted by the user.] Consequence: the §6.1
  anchor reproduces v5 EXACTLY (algorithm + numbers); the golden master is the
  direct numeric anchor, no regeneration needed.
- **b-sweeps span [0.1, 1]** per user. Trim spec §5 b-axes that exceed 1
  (funder_scale had →5, seed_value/regime_map had →2/3): clamp to ≤1, note it.
- Other unswept defaults follow spec §2 (ε=0.1, x_seed=0.25) EXCEPT budget.
  Anchor regression is checked at ε=0.3, x_seed=0.5 (golden-master conditions).
- **Forward planner:** **generalize the existing** A/B/`ce_reweight` machinery
  to T rounds (NOT pure CE). Info-anticipation chaining rule: at each round t,
  reweight the posterior via `ce_reweight` for the immediately-next round only
  (one-step info anticipation), CE point-estimates for rounds beyond. At T=2
  this collapses to the current `forward_A`/`forward_B` pair → exact reproduction.

## Corrections to inherited context
- Transfer-message model formulas were WRONG. Actual v5 (app.R:107-114):
  `λ = γ·KR/(K+R)`, `K_next = K + ε·K·R/(K+R)`. Spec §4 matches the real code.
- `simulation.R` is a stale **v3** snapshot, not a headless core — ignore it;
  app.R (v5) is the source of truth.
- Convention docs cited by the spec (`RSHINY_STYLE_GUIDE.md`,
  `PSEUDOCODE_CONVENTIONS.md`) are absent from the repo — inferring house style
  from app.R/sweep.R.

## Resolved: §6.1 PG anchor was wrong (n=50, 200 seeds)
- Spec §6.1 asserts `fwd_vs_myo_PG` (S8−S5) ≈ **−0.33** at (ε=0.3, τ_k=1, b=0.5).
  Measured on current v5: **−0.043 ± 0.038 — NOT distinguishable from 0.**
  The −0.33 is a single-seed artifact (golden-master seed 2 = −0.358), not the
  mean. **Correct the §6.1 assertion** to `|PG| < 0.1, weakly negative`; do NOT
  assert PG ≈ −0.33 (would fail regression on correct code).
- The other three §6.1 anchors validate: signal_fwd = +6.17 ✓ (spec ≈6.5),
  seed_fwd = −0.28 < 0 ✓, alpha_S8 = 0.464 < 0.5 ✓.
- **New-regime anchor** (ε=0.1, new median/b=0.1 budget, x_seed=0.25), the
  reference to write into regenerated §6.1: PG = **−0.026 ± 0.006** (sig, tiny),
  signal_fwd = **+2.72**, seed_fwd = **−0.049**, alpha_S8 = **0.494**.
- §6.10 headline intact: forward ≈ ties myopic at T=2 (~0), so a positive
  turn as T grows is a clean sign-flip story (not "recovery from −0.33").

## Focused PG(T) run (n=50, 50 seeds, tests/pg_focused.R) — headline is REGIME-DEPENDENT
- Earlier 8-seed smoke ("PG grows to +1 at T=5") was NOISE. Proper 50-seed run:
  PG=S8−S5 is NON-MONOTONIC in T and depends on growth ε:
    * ε=0.1: PG turns robustly NEGATIVE at high T (−0.30 z−2.9 T4; −0.49 z−2.7 T5).
    * ε=0.3: peaks +0.36 (z2.7) at T3, then declines to ~−0.2.
    * ε=0.6: stays POSITIVE, peaks ~T3 (+0.43), +0.40 (z2.6) at T5.
  PG generally peaks near T=3. Whether it stays positive is ε-gated.
- ROBUST + clean: signal_fwd (S8−S7) rises monotonically with T for all ε
  (4.3→12.2). b_idx_S8 rises with ε (front-loading, Force B).
- CAVEAT: T=2 bit-identical to v5, but T≥3 is the CE generalization. Negative
  PG at high-T/low-ε may partly be CE/Jensen tax accumulating over rounds
  (§6.9 CE-vs-SP is meant to bound this), not a property of optimal planning.
- Implication: the Tier-1 horizon_growth sweep (T×ε, 200 seeds) is exactly the
  right instrument to map this regime boundary. §6.10 answer = regime-dependent.

## Phase status
- [x] **Phase 0 — Golden master.** `tests/capture_golden_master.R` →
  `tests/golden_master_v5.rds` (30 rows: 6 param points × 5 seeds, per-strategy
  g1/g2/out/alpha). Normalizer-independent regression anchor.
- [x] **Phase 1 — T-round core** (`simulate_T.R`). DONE + validated.
  Gate PASSED: `tests/test_T2_reduction.R` shows T=2 == v5 bit-identically
  (max |Δ out| = 4e-14, Δg1 = 0, Δα = 0; 7 pts × 5 seeds × 9 strat).
  §3 copula was already in v5 (`draw_initial_population`); reused it.
  Bugs found & fixed: S1/S2 allocated 0 (n from NULL posts → derive from
  p_prev); forward last-round dropped K-compounding history; S9 overspent at
  T=1 (seed double-counted in H=1 greedy). Budget conservation ≤3e-14.
  Preview headline: PG(T)=S8−S5 grows 0 → +0.68(z4.2,T=4) → +0.97(z4.3,T=5).
  NOTE for Phase 3: myopic b_idx dips <0.5 at high T (FP leftover, tranche not
  a multiple of delta; α_t still equal) — assert "α_t equal", not b_idx≡0.5.
  NOTE: heavy-tail (α_K→1) rpois can produce NA p_t (huge λ) — inherited from
  v5; §6.6 winsorization to address in Phase 3.
- [x] **Phase 2 — Sweep runner + §7 schema** (`sweep_T.R`). DONE + validated
  end-to-end (resource_noise tiny run: full §7 schema present, raw+rawlong+
  summary+plots saved). 188 cells post b-trim (b clamped to [0.1,1]). CRN
  across cells (seed=trial) instead of §7 hash — reduces cross-cell contrast
  variance (documented in file header). Parallel via mclapply, checkpointed.
- [x] **Phase 4 — Micro-benchmark (§8)** → benchmark_report.txt. per-run
  T=1..5: 0.022/0.070/0.145/0.255/0.403 s (O(T²), cost5/cost2=5.7). Full
  manifest 37,600 runs ≈ 1.5h serial, ~0.25h on 7 cores. GO — no dials.
- [x] **Phase 3 — Assertions (§6)** → assertions_log.txt. ALL PASS (1 anchor,
  2 T1, 3 eps0, 4 tauK limits, 6 heavytail-finite, 7 budget, 8 M-convergence,
  9 CE-vs-SP, 10 sign-path). #9 CE valuation gap at T=2 = −0.24% (near-exact).
  #6 flagged: ESS/M<0.1 for ~25% of researchers at α_K=1.1 (IS degrades at
  extreme tails; heaviest-tail cells read with mild caution).
  FOLLOW-UP: extend #9 to T=3–5 to test if CE tax grows with horizon (would
  explain the PG(T) peak-at-T=3-then-decline as partial CE artifact).

## HEADLINE RESULT (horizon_growth, 200 seeds) — clean (T×ε) phase diagram
PG = S8−S5 (forward−myopic), z in []:
  eps=0.05: T3 +0.09  T4 −0.32[z−4.7]  T5 −0.57[z−7.9]   (forward LOSES at high T)
  eps=0.30: T3 +0.28[z4.9]  T4 −0.14   T5 −0.31[z−4.5]
  eps=0.55: T3 +0.52[z8.8]  T4 +0.31   T5 +0.25[z3.7]    (stays positive)
  eps=0.85: T3 +1.02[z17]   T4 +0.97   T5 +1.05[z13]     (strongly positive)
- Regime boundary between ε=0.30 and 0.55 at high T. PG peaks at T=3 for all ε.
- MECHANISM: b_idx_S8 (front-loading) rises with ε (0.48@ε.05 → 0.59@ε.85);
  forward front-loads to exploit compounding only when growth rewards it.
- signal_fwd (S8−S7) rises monotonically with T for all ε (3.9→11.7); larger
  at low ε (signal matters more when compounding doesn't dominate).
- §6.10 answer: forward beats myopic iff growth ε is high enough to reward
  front-loading; the advantage is a (horizon × growth) regime, not monotone-in-T.
- [x] **Phase 5 — Full manifest** → `sweep_results/T_run/` (13 sweeps, 200
  seeds, 22.3 min wall). All §0 deliverables complete.

## FULL RESULTS (200 seeds) — all 13 sweeps
1. horizon_growth (T×ε): PG phase diagram (above). Forward>myopic iff ε high;
   boundary ε∈[0.30,0.55]; PG peaks T=3.
2. horizon_noise (T×τ_k): signal_fwd ↑ with T, ↓ with τ_k (all clean). PG shows
   same T=3-peak-then-decline at ALL τ_k (ε=0.3 regime).
3. horizon_scale (T×b): same PG shape at all b; decline robust to budget.
4. signal_precision (τ_k): signal_fwd monotone 7.4(τ.05)→0.23(τ20); half≈τ1.5.
5. signal_value (k_shape×τ_k): KEY interaction. signal_fwd 20.9 (heavy tail
   k1.3 + sharp τ.3) → ~0 (light tail k3.5 + noisy). Signal matters most at
   high inequality + precise signal.
6. funder_scale (b): PG≈0 at T=2 all b (horizon, not budget, drives forward gain).
7. seed_value (b×x_seed): seed_fwd ALWAYS <0 (−0.02→−0.70), worse at high b/x_seed.
   Uniform seeding never helps the planner. Clean negative result.
8. alpha_regime (τ_k×ε): b_idx_S8 driven by ε (0.50→0.53), ~independent of τ_k.
   Front-loading = compounding (Force B), not signal precision.
9. pre_rounds (n_pre): signal_fwd RISES with baseline data (6→17.6 at n_pre20).
   CAVEAT: partly a scale artifact (naive pre-round funding inflates K,R,λ so
   absolute contrasts grow); interpret relative, not absolute.
10. regime_map (b×τ_k×k_shape): 3D confirmation of #5 (see raw).
11. correlation (ρ_c×k_shape, τ_k=0.3): PROMOTED RESULT CONFIRMED. signal_fwd ↓
    with ρ_c (22.8→15.9 at k1.3; 8.9→3.1 at k2) but REMAINS LARGE at ρ_s=0.78.
    Signal dominance survives realistic positive K–R correlation. rho_s tracks
    rho_c (−0.46/0.01/0.48/0.78).
12. pop_size (n): signal_fwd scales with n (2.3→42); PG≈0 at T=2 all n (no
    finite-n artifact).
13. resource_noise (τ_r): signal_fwd flat (~6) — resource-signal noise doesn't
    affect grant-signal value (separate channels). Clean null.

## ⚠️ MAJOR FINDING: high-T PG results are a GREEDY-DISCRETIZATION ARTIFACT
- CE-tax diagnostic (tests/ce_tax_vs_T.R, n=40) found CE tax negligible (~0,
  ±0.1% of value, no growth with T) → decline is NOT a CE artifact.
- BUT it also revealed the decline is n-FRAGILE. PG(T=5, ε=0.3) vs n (200 seeds):
  n=30 −0.39[z−8] · n=40 +0.80[z+13] · n=50 −0.31[z−4.5] · n=60 −0.40[z−5] · n=80 +1.94[z+19].
  Sign flips ERRATICALLY with n, all high-|z|. This is the signature of a
  discretization artifact: the forward planner allocates δ=B/n_steps lumps over
  T×n slots; at high T the lumps are coarse vs the slots and land pathologically
  depending on n. Confirmed by granularity_check.R (n_steps sweep) — [pending].
- IMPLICATION: Tier-1 horizon results at T≥3–4 (incl. the phase-diagram decline)
  are NOT reliable as-is. T=2 (bit-identical to v5) is solid. ALL the T=2
  cross-sections — signal_value, signal_precision, correlation, seed_value,
  funder_scale, resource_noise, pop_size@T2 — are UNAFFECTED and stand.
- FIX if confirmed: re-run Tier-1 horizon sweeps at much finer n_steps
  (200–400) + more seeds; re-issue the phase diagram. The visual report's
  Fig 1 must carry this caveat or be regenerated.

## ⚠️⚠️ ESCALATION: high-T forward pathology is DIVERGENT & grant-signal-specific
- convergence_check.R (n=50, T=5, 30 seeds): PG=S8−S5 does NOT converge with
  granularity — it DIVERGES negative: eps0.3 {50:−0.41, 200:−2.20, 800:−5.39};
  eps0.85 {50:+1.07, 200:−2.01, 800:−7.22}. The phase diagram's strongest
  POSITIVE cell (eps0.85) flips to −7.2 at fine granularity. Whole horizon
  headline = coarse-granularity artifact.
- LOCALIZED: PGpubs=S7−S4 (forward vs myopic, NO grant signal) is STABLE across
  granularity (eps0.85: +1.72/+1.55/+0.87; eps0.3: ~−0.5 flat). So the multi-
  round forward logic is fine WITHOUT the grant signal. The divergence is
  specific to the grant-signal-informed posterior in the T≥3 forward planner
  (ce_reweight / fwd_researcher_value path). Divergence w/ granularity ⇒ likely
  a BUG in that path (or a genuine CE instability), NOT a bounded economic effect.
- STATUS: horizon/Tier-1 results RETRACTED pending diagnosis. Visual report Fig 1
  must be pulled/caveated. T=2 + all T=2 cross-sections (signal_value,
  correlation, seed, precision, funder, resource, pop@T2) remain SOLID.
- NEXT: diagnose bug-vs-real in the grant-signal T≥3 forward path (check b_idx_S8
  / schedule at fine granularity; inspect ce_reweight + fwd_researcher_value with
  a sharp K posterior over >2 rounds).

## ✅ ROOT-CAUSED: stochastic SIR resampling in ce_reweight_posterior (app.R:302-311)
- BUG-vs-REAL test (/tmp/diag.R): at n_steps=800 the forward greedy scores LOWER
  on its OWN CE objective (V_fwd=191.1) than the myopic schedule (V_myo=191.8) —
  impossible for a correct optimizer ⇒ BUG, not CE bias.
- ROOT CAUSE: ce_reweight_posterior does stochastic SIR resampling when reweighted
  ESS<M/4. fwd_marginal calls it twice (pb_base at g1, pb_plus at g1+dg); when SIR
  fires they resample DIFFERENT random atoms, so the finite-difference marginal is
  corrupted by O(1) resampling noise that dominates as dg→0. Grant signal sharpens
  the posterior (ESS/M 0.57 vs 0.75) → SIR fires → grant-signal-specific. Amplifies
  with T (more downstream rounds at the noisy reweighted posterior).
- FIX VERIFIED: deterministic reweight (no SIR) → PG converges & is granularity-
  stable. eps0.3 T5: +0.28/+0.28/+0.30 (n_steps 50/200/800); eps0.85: +1.55/+1.91/
  +1.89. REAL RESULT: forward BEATS myopic at high T, more at higher ε (no decline).
- Bug is inherited from v5 (v5 forward uses the same ce_reweight+SIR); mild at
  T=2/coarse, catastrophic at T≥3/fine. Fixing app.R fixes v5 Shiny + T-round both.
- PENDING USER: where to apply fix (app.R shared primitive vs simulate_T override);
  then re-establish anchor vs FIXED v5, re-run manifest, regenerate phase diagram.

## ✅ FIX APPLIED (app.R ce_reweight_posterior = deterministic reweight; M 200→400)
- Decisions: fix in app.R (shared, fixes Shiny v5 too); deterministic reweight
  + bump M to 400 for ESS headroom.
- RE-VERIFIED on fixed code:
  * T=2 anchor vs FIXED v5: still bit-identical (Δ=4.3e-14). Golden master
    re-captured (numbers shift slightly — bug was mild at T=2).
  * Convergence: PG(T=5) now granularity-STABLE & positive. eps0.3 +0.29/+0.27/
    +0.30; eps0.85 +1.62/+1.96/+1.96 (n_steps 50/200/800). Bug gone.
  * All §6 assertions PASS. #3 reframed: at eps→0 schedule is uniform (b_idx≈0.5)
    but forward still beats myopic by ~0.2 = pure INFORMATION value (previously
    masked by the SIR noise). New decomposition: forward value = compounding
    (∝ε) + information value (survives ε→0). CE-vs-SP gap still −0.26% at T=2.
- Manifest re-running (fixed, M=400) → sweep_results/T_run_fixed/. Old buggy run
  kept in sweep_results/T_run/ for comparison.
- TODO after manifest: regenerate phase diagram from T_run_fixed, update artifact
  Fig 1 (now a real, stable, POSITIVE forward-advantage story).

## ✅ COMPLETE: corrected manifest + artifact (sweep_results/T_run_fixed/, M=400, 29 min)
- CORRECTED PG=S8−S5 (200 seeds) — positive EVERYWHERE, grows with ε:
  eps0.05: ~0 (T2) .. +0.23(T3) .. ~0 (T4-5)   [≈ pure info value, peaks T3]
  eps0.30: +0.04 .. +0.37 .. +0.19 .. +0.24
  eps0.55: +0.12 .. +0.60 .. +0.62 .. +0.72
  eps0.85: +0.27 .. +1.08 .. +1.28 .. +1.66[z24]
  Mechanism: b_idx_S8 (front-loading) rises with ε (0.51→0.62) and T. No decline.
- T=2 cross-sections essentially UNCHANGED vs buggy run (bug was mild at T=2):
  signal_value 21.1 (heavy tail+sharp) → ~0; correlation signal survives ρ_s=0.78
  (22.9→15.6 at k1.3); seed_fwd <0 everywhere; precision 7.7→0.09 monotone.
- Artifact updated (same URL, label "corrected-planner"): Fig 1 = sequential
  teal, positive-everywhere; strip/methods/notes reflect the fix. Other figs
  refreshed with M=400 numbers.
- Old buggy run retained at sweep_results/T_run/ for provenance.

## HEADLINE (final): forward planning beats myopic across the board; its edge grows
## with knowledge-growth ε (via front-loading to compound) and, at high ε, with the
## horizon. Decomposes into compounding value (∝ε) + information value (survives ε→0).
## Signal value governed by inequality×precision, survives K–R correlation. Seeding
## never helps. All validated; T=2 bit-identical to (fixed) v5.

## OPEN (superseded): post-T=3 PG decline — real vs CE artifact
- Decline after T=3 is universal in τ_k, b (at ε≤0.3) but ABSENT at high ε
  (ε=0.85 plateaus positive). So there IS real ε-structure, not pure artifact.
- #9 CE-tax near-exact at T=2 (−0.24%). Whether it grows with T is UNTESTED.
- To settle: build SP-lite (scenario-lookahead) planner, compare PG_SP vs PG_CE
  at T=3–5 (spec §6.9 extended). Not yet built — flagged for user decision.
