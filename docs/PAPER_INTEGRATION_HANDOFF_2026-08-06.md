# Paper integration handoff — everything from the 2026-08-04/06 simulation campaign

*For the writing session. Self-contained: what was run, what the statistics say, where every
file lives, and how to double-check any number. All results are on branch `smooth-allocator`,
pushed, across four commits: `e172941` (poverty/exploration suite), `d77c854` (diagnostics code),
`21b0d9a` (Packages A–C data + memo), `662a597` (Packages D–E + bootstrap verification).*

## 0. Due diligence — how to check anything in this document

1. **Every headline number below was re-derived from the data files** (not session memory) by
   `sweep_results/_probe/verify_all_claims.R` on 2026-08-06. Rerun it anytime:
   `Rscript sweep_results/_probe/verify_all_claims.R` — it prints each claim's number straight
   from the canonical RDS/CSV files.
2. Every output directory contains `RUN_INFO.txt` (git hash, date, seeds, sessionInfo).
3. Model-code integrity: three validation gates were green after every one of the six model
   edits — (a) `tests/test_T2_reduction.R` (T-round path reproduces the v5 two-round model
   bit-identically at T=2); (b) the pinned run `tests/pinned_T5_baseline.rds` (seed=1, T=5,
   ε=0.3, all 9 strategies, smooth — bit-identical at default arguments); (c)
   `tests/assertions_T.R` (8 checks). Every new mechanism is behind a default that reproduces
   the published model **bit-identically** (not approximately).
4. The two primary results memos, which this document summarizes and supersedes in organization
   (not in content): `docs/DIAGNOSTICS_RESULTS_2026-08.md` (all preregistered verdicts) and
   `T_round_extension/RESOURCE_REGIME_RESULTS.md` (the poverty/bootstrap story, two parts).
5. Conventions everywhere: 9 strategies S1–S9 as in the paper (+ new S10/S11, see §2); paired
   within-seed contrasts; `PG = S8−S5` (value of planning), `signal_fwd = S8−S7` (value of the
   review signal), `b_idx` = budget-schedule center of mass (0.5 even, <0.5 front-loaded,
   >0.5 back-loaded); "% of S1" = normalized by the same cell's no-funding output.

## 1. The six model extensions (all in `model.R`, all default-dormant)

| Switch | What it does | Used by |
|---|---|---|
| `budget_ref = "R"\|"K"` | `"K"` sets B = b·n·E[K]: funder's purse decoupled from community resources (default ties B to E[R], so r_min→0 would zero the budget — a trivial null) | resource_regime, exploration suite |
| S10/S11 + `seed_persistent` | Seed-floor variants of the signal strategies (S10:S5, S11:S8) and an every-round floor; at x_seed=1 persistent S6 ≡ S2 **exactly** (verified, Δ=0) | Package A |
| `ces_gamma` (+ `ces_mean`) | CES production family; −1 harmonic (bit-identical fast path), 0 Cobb-Douglas, −∞ Leontief; smooth allocator guarded harmonic-only; non-harmonic runs use greedy n_steps=400 (Leontief 800) | Package B |
| `eps_free`, `eps_paid` | Growth decomposed: K += ε_free·λ(K,R₀) + ε_paid·(λ(K,R₀+g)−λ(K,R₀)); threaded through realized dynamics AND all T-round planner anticipation (8 call sites); equal rates take the exact old code path | Package C |
| `tau_k_true`, `tau_k_belief` | World draws the review signal at the true noise; the funder's likelihood weights it at the believed noise | D-2 |
| `oracle`; `A_obs_sdlog`; `latentA_cfg` (+ per-atom `post$gam`) | Point-mass posteriors at truth (D-4); observable heterogeneous productivity A_i (E-1a); latent-A atoms with variance split of log-talent (E-1b) | D-4, E-1 |

Scope note: the ε-decomposition and heterogeneity cover the **T-round engine** (the sweep
engine); the v5 two-round anchor path supports defaults only, by design.

## 2. Results, organized by where they land in the paper

### 2.1 Part V extension — "The bootstrap objection" (poverty, paid information, front- vs back-loading)

**The objection to answer:** shouldn't a funder facing a resource-poor community front-load —
invest early to bootstrap growth and learning? (At R₀=0, no output → no learning, no growth;
money in round T buys nothing for T−1 rounds.)

**The logical correction to present first:** at R₀=0 the two lump schedules are payoff-identical
— all-in-round-1 also allocates blind (spending precedes observation), and its learning/growth
are worthless with no later money to exploit them. The argument establishes a complementarity
(early enables, late exploits): positive round-1 spend, interior schedules — not early mass.
Money follows information, and information arrives late.

**Result 1 — the front↔back reversal is governed by ε, not resources.**
`sweep_results/resource_regime/` (200 seeds; the 64-seed original in
`sweep_results/T_run_smooth_supplement/` matches to |Δb_idx| ≤ 0.0011). The b_idx = 0.5
boundary is essentially vertical at **ε\* ≈ 0.02** (range 0.005–0.05 over r_min 0.001–3); all
10 front-load cells have ε ≤ 0.03. Poverty mutes back-loading and raises the round-1 share
(at ε=0.85: b_idx 0.521→0.652, round-1 share 0.174→0.060 as r_min goes 0.001→3) but never
flips the sign: the funder's own grants manufacture the compounding that rewards spending late.

**Result 2 — thin grants are talent-uninformative (new mechanism).**
With λ = γKg/(K+g), a grant g ≪ K pins output to the grant (λ ≈ g): a K=1 and a K=100
researcher are nearly indistinguishable at g=1/3. Consequence, from `exploration_corner`
(no free signal, poverty, heavy tail α_K=1.3, T=6): at standard depth b=0.5 the funder's
discrimination value S8−S2 is **+0.5** vs **+17.5** with a sharp free signal. Paying for
information requires **depth**: g ≳ K, i.e. b ≳ T/2.

**Result 3 — depth makes paid information work, and it back-loads the money.**
`sweep_results/exploration_200/exploration_depth_summary.rds` (200 seeds): S8−S2 rises
0.5 → **41.2** (b=3) → **56.8** (b=6); the schedule becomes **seed-and-harvest** — round-1
share **0.110** at b=3 (even = 1/6 ≈ 0.167), mass deployed once informed. Across all 32
exploration cells min b_idx = **0.497**: the ε≈0 cells sit statistically below 0.5 at 200-seed
precision but by ≤ 0.003 — economically nil and confined to where planning value ≈ 0.
**Slogan: learning front-loads the observation and back-loads the money.**

**Result 4 — attribution: which growth channel drives back-loading?**
`sweep_results/bload_decouple/` (Package C, 200 seeds, T=5). Pure free growth (ε_paid=0):
b_idx = 0.529/0.587/**0.680** at ε_free = 0.1/0.3/0.85. Pure paid growth (ε_free=0):
b_idx = **0.489/0.481/0.481** — genuine front-loading, strictly below 0.5 (SE ≤ 0.003), and
mildly valuable at high ε_paid (PG = +0.19). The diagonal reproduces the coupled model to 4
decimals (0.5537, 0.6192 vs `horizon_growth`). Transects locate the boundary: back-loading
takes over once **ε_free/ε_paid ≈ 1/8–1/3** (crossing at ε_free ∈ (0.05, 0.10) at ε_paid=0.3;
(0.10, 0.15) at ε_paid=0.85). The free channel dominates at remarkably small relative strength.

**Result 5 — the free forces dominate even when switched off (the B→C, D→E symmetry).**
At depth, the forward planner's deliberate schedule loses to myopic re-deciding (S8−S5 = −16 at
b=3, −45 at b=6): the myopic funder re-decides on posteriors updated by its own funding's
output, capturing paid information automatically. Paid growth (B) regenerates free growth (C);
paid observation (D) regenerates free observation (E).

**Two independent verifications (present both as robustness):**
- *Planner-free certification* (`sweep_results/bootstrap_verify/honest_schedules_*.csv`):
  a grid of two-block schedules executed with the myopic within-round allocator and NO CE
  planner, at the depth corner, 200 seeds. The even schedule is best at ε≈0 (498.2; best
  early-mass schedule loses by 7.8); at ε=0.3 the top schedule ties even within noise and the
  best early-mass loses by 20.8. Seed-and-harvest is not a planner artifact.
- *CE self-consistency* (`bootstrap_verify/ce_self_consistency.csv`): CE(S8's plan) ≥ CE(even-
  tranche plan) under S8's own objective in **20/20 seeds** (+3.33 ± 0.11). The negative PG at
  depth is honest mispricing of information by the CE objective, **not an optimizer bug** —
  state this caveat wherever S8−S5 < 0 is discussed.

Also useful: `exploration_poverty` (ε≈0): with the signal off, the *rich* community defers
(round-1 share 0.133, money-follows-information), the poor community cannot (flat, nothing to
learn from λ≈g output); with the signal on, timing is irrelevant. And the information value of
funding is non-monotone in baseline resources (peaks at r_min≈1: S8−S7 ≈ 10.8, vs 2.3 poor /
6.6 rich): too poor → nothing to observe; too rich → free observation already reveals talent.

### 2.2 Part VI rewrite — seed floors are NOT nearly free (Package A)

Data: `sweep_results/D3_seed_signal/`, `D4_seed_persistent/`. The old "floors cost <0.4%"
number came from flooring only the weakest optimizer (no signal), round 1 only.
- **P-A1 CONFIRMED:** floor the strongest optimizer (S11 vs S8) where targeting matters
  (α_K=1.3, τ_K=0.3, b=0.5, x_seed=0.75, T=2): cost = **6.11% of S1** (raw 5.71). Cost ≈
  κ × (family targeting value) × (floored budget share), κ ∈ [0.10, 0.86] — κ in (0.3, 1] at
  tight budgets/heavy tails, below 0.3 only at ample-budget small-floor cells. The floor
  displaces 2.8–43% of the S8 allocation (`disp_S11`).
- **P-A2 CONFIRMED (exact):** persistent floor at x_seed=1 reproduces S2 with Δ=0 in outputs
  and grant vectors (implementation identity).
- **P-A3 CONFIRMED:** cost convex in x_seed on both designs (round-1-only: quad coef −0.06 to
  −0.24 by b; persistent: −4.86 heavy / −1.76 base tail).
- Paper line: **"floors cost what targeting is worth, times the floored share"** — small only
  where targeting itself is worth little (ample budget, weak signal, light tail).

### 2.3 Part V numbers — the T-scaling exponent (A5)

From `T_run_smooth/horizon_growth` + `horizon_long`: log PG on log T gives exponent
**2.33 (SE 0.12)** at ε=0.3 and **2.05 (0.15)** at ε=0.85 for T ≤ 5 — "approximately
quadratic" is correct as a small-T statement — but **1.44 (0.03)** and **0.78 (0.05)**
respectively over T = 5–10: planning value saturates at long horizons. Cite the exponent, and
scope the quadratic claim to T ≤ 5.

### 2.4 New robustness section — the CES family (Package B)

Data: `sweep_results/sigma_tierA_gc0/`, `sigma_tierA_gcm3/`, `sigma_tierA_leontief/` (Tier A,
greedy n_steps=400, 100 seeds; Leontief 50 seeds, n_steps=800, degenerate-kink caveat),
`sigma_tierB_concentration/` (Tier B, T=1, S5, 200 seeds + 400-seed refinement).
- **P-B1 MIXED — the headline exception is itself a result:** the information story survives
  the whole family (signal value monotone in inequality × precision at Cobb-Douglas, γ=−3, and
  Leontief; seed cost small in the no-signal family everywhere, worst −0.44% of S1; correlation
  robustness holds), **but Cobb-Douglas kills the planning phenomenon entirely**: PG(T=5) ≤
  0.014 at ALL ε ≤ 0.85 and b_idx ≤ 0.512, vs γ=−3 where planning strengthens (PG = 3.77,
  b_idx = 0.660 at ε=0.85). **Back-loading and the forward advantage are complementarity
  phenomena (γ_ces < 0)** — scope Part V accordingly; arguably a headline robustness finding.
- **P-B2 REFUTED (the gate did its job):** concentration of the optimal static allocation is
  NOT monotone in substitutability. The direction is **budget-conditional**: tight budgets
  concentrate toward Leontief (α_K=1.3, b=0.1: Gini 0.916→0.938, Leontief 0.947); ample
  budgets equalize toward Leontief (α_K=2, b=1: 0.354→0.285, Leontief 0.207). And there is a
  genuine **interior maximum** at light tails/ample budgets (α_K=3.5, b=1: peak at γ≈−6,
  falling to −12 with z=3.5, refined at 400 seeds). **The "concentration dispute as a
  disagreement about σ" framing dies as stated.** The replacement is clean: *harder
  complementarity concentrates funding only when money is tight; with an ample budget it
  equalizes (top-up-to-the-bottleneck covers everyone).* P-B3 (coverage) mirrors this exactly.

### 2.5 Part IV additions — information integrity (Package D)

- **D-1, the empirical anchor quantified** (`sweep_results/D_auc_calibration/`, pure Monte
  Carlo): the model's implied τ_K at Fang et al.'s AUC = 0.54 is **> 20** at α_K ∈ {1.3, 2}
  (AUC at τ=20 is still 0.60 at the heavy tail) and 7.4–10.3 at α_K=3.5 — all far above the
  signal-value half-decay elbow (τ ≈ 2.5). Footnote 13's reading survives with a wide margin;
  note the corollary that an AUC as low as 0.54 under heavy tails implies *extreme* noise.
- **D-2, the overtrust caveat curve** (`D_misspecified_trust/`): asymmetry confirmed —
  overtrust costs more than undertrust forfeits — but the sign flip is tail-dependent. At
  α_K=2, overtrusting a noisy signal (truth τ=3, belief ≤1) turns review's value **negative**
  (−0.53/−0.43/−0.14; individually |z|≲1.6, consistent across cells). At α_K=1.3 it merely
  attenuates (13.99 vs 17.20 calibrated): in heavy-tailed fields even misweighted review
  separates the whales — Part IV.2's coarse-signal logic extends to misspecified weighting.
- **D-3, the report's own open item closed** (`D_resource_ablation/`): the resource signal's
  own contribution to S8 is **0.3–1.1% of S1 everywhere** (including ρ=−0.5, its predicted
  best case). The redundancy claim stands as stated. (Side observation left uninterpreted:
  for pubs-only S4 the contribution is larger and flips sign at ρ=0.8.)
- **D-4, the thesis figure Part I→IV** (`D_gap_convergence/`): corr(realized grants,
  max(cK−R,0)) rises monotonically as the signal sharpens — 0.34→**0.97** (α_K=1.3),
  0.22→**0.95** (α_K=2) — recovering footnote 3's r = 0.95–1.00; the output gap to the
  full-information oracle shrinks 17.2%→1.7% of S1; the pubs-only baseline plateaus at
  corr 0.13–0.18. This is "the cost of not knowing the gap," measured.

### 2.6 Limitations/appendix — heterogeneous productivity (Package E-1)

Data: `sweep_results/E1_heterogeneity/`.
- **E-1a (A observable) CONFIRMED:** grants track the productivity-scaled rule
  g*_i = max(K_i(√(γ_i/ν)−1) − R_i, 0) at **r = 0.87–0.98** (≥0.95 wherever τ_K=0.3); all
  headline directions survive (signal value 19.6→18.1 as sdlog 0→0.5 at heavy+sharp; seed
  cost stays negative; PG ≈ 0 at T=2). The closed form generalizes as claimed.
- **E-1b (A latent) — the product-reduction hypothesis FAILS materially:** splitting log-talent
  variance between K (share s) and latent A at fixed E[log T], Var(log T): signal value
  collapses **11.8 → 3.7 (−69%)** as s goes 1→0.4 (heavy tail, sharp signal; −83% at
  moderate/τ=1); who-gets-funded reshuffles (rank corr vs s=1 down to **0.72**); Gini and
  output shift beyond noise. Registered design choice to state: the review signal observes K
  only, so as the K-share of talent shrinks, review sees a shrinking share of what matters —
  that is part of what the reduction would have had to survive. **Conclusion for the paper:
  one-dimensional talent is a substantive commitment, not a convenience; A and K are
  separately identified only by dose-response (output vs funding depth), never by
  cross-sectional output.** This strengthens the empirical-calibration agenda.

### 2.7 Scale insurance (E-2) — cite freely at n=50

`sweep_results/E2_headline_scale{,_n200}/`: the (T×ε) forward-advantage surface at n=200
correlates 0.996 (PG as % of S1) and 0.999 (b_idx) with n=50 across all 25 cells (relative
magnitudes ~10–15% smaller); at ε=0.85 the advantage keeps sign and T-monotonicity at every
b ∈ {0.1, 0.5, 1}, scales roughly with the purse (0.06%→1.52% of S1 at T=5), and the schedule
is b-invariant (b_idx 0.60–0.62 at T=5).

### 2.8 Package F — removed by design (do not re-propose)

Grant persistence is deliberately not modeled: grants are consumed; what persists (skills,
publications, reputation) is exactly K, and the paid-growth channel ε·λ(K, R₀+g) already
carries the durable residue. State the consumable-grant assumption in the setup as a
substantive claim; one limitations line for endowments and capital-intensive fields (the
latter folds into the Leontief end of the CES discussion).

## 3. One-line verdict table (every preregistered prediction)

| ID | Claim tested | Verdict |
|---|---|---|
| P-A1 | Floor on strong optimizer costs 2–8% of S1 | **CONFIRMED** (6.11%) |
| P-A2 | Persistent x_seed=1 ≡ uniform | **CONFIRMED** (exact, Δ=0) |
| P-A3 | Floor cost convex in x_seed | **CONFIRMED** (both designs) |
| A5 | PG ~ T² | **Scoped**: exponent 2.1–2.3 for T≤5; 0.8–1.4 for T=5–10 |
| P-B1 | CES directions survive | **MIXED**: info/seed/corr yes incl. Leontief; planning vanishes at Cobb-Douglas |
| P-B2 | Gini monotone in σ | **REFUTED**: budget-conditional + interior max (γ≈−6, z=3.5) |
| P-B3 | Coverage monotone in σ | **REFUTED**: mirrors P-B2 |
| P-C1 | Pure free growth back-loads | **CONFIRMED** (b_idx→0.680) |
| P-C2 | Pure paid growth front-loads | **CONFIRMED** (b_idx 0.481–0.489; PG +0.19 at 0.85) |
| P-C3 | Diagonal ≡ coupled model | **CONFIRMED** (bit-identical; sweep-level to 4 dp) |
| P-D1 | Implied τ(AUC=0.54) above elbow | **CONFIRMED** (>20 heavy/moderate; 7–10 light; elbow 2.5) |
| P-D2 | Overtrust negative, undertrust attenuates | **MIXED**: asymmetry yes; negative only at moderate tails |
| P-D3 | Resource signal redundant | **CONFIRMED** (≤1.1% of S1 everywhere) |
| P-D4 | Gap-rule convergence as signal sharpens | **CONFIRMED** (r→0.95–0.97; oracle gap→1.6%) |
| P-E1a | Scaled rule tracks; directions survive | **CONFIRMED** (r 0.87–0.98) |
| P-E1b | Product reduction ≈ valid (≤10–15% shifts) | **REFUTED** (−69–83% signal value; rank corr 0.72) |
| P-E2 | n- and b-robustness of headline surface | **CONFIRMED** (corr 0.996/0.999) |
| Boot-1 | Early-mass schedule beats even at depth? | **NO** (planner-free; loses by 7.8–20.8) |
| Boot-2 | Negative PG at depth a bug? | **NO** (honest CE mispricing, 20/20) |
| Boot-3 | 64-seed results stable at 200 | **YES** (|Δb_idx| ≤ 0.004) |

## 4. File map

```
docs/DIAGNOSTICS_RESULTS_2026-08.md      all verdicts, both handoffs (the primary memo)
docs/SWEEP_HANDOFF_2026-08-05.md         Packages A-C spec (preregistrations)
docs/SWEEP_HANDOFF_ADDENDUM_2026-08-05.md Packages D-E spec + F removal rationale
docs/PAPER_PROMPT_resource_regime.md     earlier prompt for the poverty story (superseded by this)
T_round_extension/RESOURCE_REGIME_RESULTS.md  poverty/bootstrap write-up (Parts 1+2 + provenance)
sweep_results/
  resource_regime/                200-seed r_min x eps map (use for figures)
  T_run_smooth_supplement/        64-seed originals (match to <=0.0011)
  exploration_200/                200-seed corner/poverty/depth (use for figures)
  bload_decouple/                 eps_free x eps_paid grid + contour transects
  D3_seed_signal/ D4_seed_persistent/    Package A
  sigma_tierA_gc0/ _gcm3/ _leontief/     Package B Tier A (per-gamma dirs)
  sigma_tierB_concentration/      Tier B + Leontief endpoint + 400-seed refinement
  D_auc_calibration/ D_misspecified_trust/ D_resource_ablation/ D_gap_convergence/  Package D
  E1_heterogeneity/ E2_headline_scale/ E2_headline_scale_n200/                      Package E
  bootstrap_verify/               honest schedules + CE self-consistency
  _probe/verify_all_claims.R      re-derives every number above (read-only)
figures/resource_regime/          5 PNGs incl. exploration_depth_schedules.png (the best figure)
```

## 5. Flags — spec contradictions and scope notes the paper must handle

1. Model write-up Eq (10) says budgets are per-round/non-transferable; the implemented forward
   strategies allocate the total remaining budget across the horizon. Fix before drafting.
2. Default-parameter tables disagree (write-up: T=10, n=100, B=50; code: T=2, n=50, b=0.5).
3. Root `RESULTS.md` predates the exact-allocator story in places; treat
   `DIAGNOSTICS_RESULTS_2026-08.md` + this document as canonical for the new results.
4. Cobb-Douglas kills planning value — scope every Part V claim to complementarity (γ_ces<0).
5. Concentration–substitutability is budget-conditional — replaces the old Story-5 framing.
6. The ε-decomposition and heterogeneity live in the T-round engine only (v5 anchor: defaults).
7. Leontief cells carry the kink-degeneracy caveat (ties broken by fill order; n_steps=800).
8. Wherever S8−S5 < 0 appears (depth cells), attach the honest-CE-mispricing caveat (Boot-2).
9. Precision honesty: "no front-loading" means min b_idx 0.497 across the exploration suite —
   statistically below 0.5 at 200 seeds in ε≈0 cells but ≤0.003 in magnitude and worthless
   (PG ≈ 0 there). Pure-paid front-loading (P-C2, 0.481) is the only strict-and-real case,
   and it requires ε_free < ~⅛–⅓ of ε_paid.

## 6. Suggested integration order

1. Part V: add "The bootstrap objection" subsection (§2.1 here, Results 1–5 + verifications).
2. Part V: cite the T-exponent with the T≤5 scope (§2.3); scope to γ_ces<0 (§2.4 flag).
3. Part VI: rewrite the floor result around P-A1's law (§2.2).
4. Part IV: add D-1 (AUC calibration), D-2 (overtrust caveat), D-4 (gap-convergence figure);
   close Appendix A's open item with D-3 (§2.5).
5. New robustness section/appendix: CES family (§2.4) — the budget-conditional concentration
   law and the Cobb-Douglas boundary are candidate headline results, not fine print.
6. Limitations/appendix: heterogeneity (§2.6) and the consumable-grant assumption (§2.8).
