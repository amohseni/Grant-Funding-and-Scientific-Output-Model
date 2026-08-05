# Diagnostics results — seed floors, CES robustness, back-loading attribution

*Execution of `docs/SWEEP_HANDOFF_2026-08-05.md`, run 2026-08-04/05 on 8 cores.
Code hash `d77c854` (diagnostics code) on top of `e172941`; every output directory
carries `RUN_INFO.txt`. All validation gates (T2 anchor, pinned T=5 run, assertions)
were green after every model edit; the pinned run is bit-identical at defaults.
Verdicts below answer the preregistered predictions; interpretation is left to the
paper session per convention 7.*

---

## Package A — seed-floor diagnostics

**P-A1 — CONFIRMED.** The floor is not "nearly free" when it dilutes a strong optimizer.
At the focal cell (α_K=1.3, τ_K=0.3, b=0.5, x_seed=0.75, T=2): cost of the floor on S8
(forward + signal) = **5.71 raw = 6.11% of S1** (prediction: 2–8%, vs the 0.4% no-signal
story). Costs scale with the family's targeting value × floored budget share:
κ ∈ **[0.10, 0.86]** across the 32 cells, inside the predicted (0.3, 1] at tight budgets
and heavy tails, falling below 0.3 only at ample budget with small floors (b=1, x_seed≤0.25,
base tail). Myopic (S10−S5) and forward (S11−S8) costs are near-equal everywhere.
The displaced-budget diagnostic D2 was cheap to add (`disp_S11`): the floor displaces
2.8%–43% of the S8 allocation, rising with x_seed. Data: `sweep_results/D3_seed_signal/`.

**P-A2 — CONFIRMED (exact).** With the persistent every-round floor at x_seed=1, S6
reproduces S2 with **zero** difference in output and in every grant vector, both parameter
regimes. The new seed code path is correct by construction. Data:
`sweep_results/D4_seed_persistent/`.

**P-A3 — CONFIRMED (both designs).** Cost is convex (accelerating) in x_seed:
- Round-1-only design (existing `T_run_smooth/seed_value`): quadratic coefficient
  −0.064 (b=0.1) → −0.237 (b=1.0), all |z|>2.
- Persistent design (D4): quad coef −4.86 (α_K=1.3) and −1.76 (α_K=2.0); cost at
  x_seed=1 equals the full targeting value S4−S2 (−2.89 and −1.33 raw), the envelope
  endpoint, reached convexly from below the linear reference.

**A5 — the T-scaling exponent (paper cites this number).** log(PG) on log(T), PG=S8−S5>0:
- T ≤ 5: exponent **2.33 (SE 0.12)** at ε=0.3; **2.26 (0.16)** at ε=0.55; **2.05 (0.15)**
  at ε=0.85 — "approximately quadratic" is right in the tested range.
- T = 5..10: **1.44 (0.03)** at ε=0.3; **0.78 (0.05)** at ε=0.85 — clear saturation;
  "quadratic" must be stated as a small-T statement.

---

## Package B — CES family

Code: `ces_mean` + `CES_GAMMA` (exponent named `ces_gamma` everywhere; the code's `gamma`
remains factor productivity), harmonic fast path keeps γ_ces=−1 bit-identical; smooth
allocator guarded harmonic-only; all non-harmonic runs greedy n_steps=400 (Leontief 800).

**P-B1 — MIXED: three of four headline directions survive across the family; the planning
results do not survive at the Cobb-Douglas boundary.**
- *Signal value increasing in inequality and precision*: **CONFIRMED** at γ_ces=0, −3, and
  Leontief. Monotone in both axes in every variant; sizes move (γ_ces=0 shrinks ~2.5×;
  −3 and Leontief grow ~1.2×; e.g. top cell 20.9 harmonic → 7.5 CD → 24.6 at −3 → 24.3 Leontief).
- *Seed floor cheap in the no-signal family*: **CONFIRMED**; worst cell −0.44% of S1
  (harmonic −0.40%, CD −0.31%, γ_ces=−3 −0.44%).
- *Correlation robustness*: **CONFIRMED**; signal value declines with ρ and remains
  substantial at ρ=0.8 in all variants (CD shows a within-noise uptick at ρ=0.8, 100 seeds).
- *Forward ≈ myopic at low ε, T=2*: **CONFIRMED** everywhere (|PG| ≤ 0.03).
- **Exception (new result): at Cobb-Douglas (γ_ces=0) the planning phenomenon vanishes
  entirely** — PG(T=5) ≤ 0.01 at ALL ε up to 0.85, and b_idx_S8 stays 0.498–0.512
  (vs harmonic 0.508–0.619, γ_ces=−3 0.511–0.660). Back-loading and the forward
  advantage are complementarity phenomena: present for γ_ces<0, strengthening as
  substitution hardens, absent at σ=1. Directions survive on the complementarity side;
  the boundary member kills them.

**P-B2 — REFUTED.** Concentration of the optimal static allocation is NOT monotone in σ,
in two ways (T=1, S5, 200 seeds; refinement at 400 seeds):
1. **The direction depends on the budget.** Tight budgets (b=0.1): Gini rises toward
   Leontief, P-B2's direction (α_K=1.3: 0.916 → 0.938, Leontief 0.947). Ample budgets
   (b=1.0): Gini FALLS as substitution hardens (α_K=2: 0.354 → 0.285 at γ_ces=−12,
   Leontief 0.207; α_K=1.3: 0.631 → 0.528, Leontief 0.448) — the opposite.
2. **An interior maximum exists** in light-tail/ample-budget cells: at α_K=3.5, b=1.0
   the refined curve peaks at γ_ces≈−6 (0.1964) and falls by z≈3.5 to 0.1852 at −12,
   continuing down to 0.166 at Leontief; at α_K=3.5, b=0.5 the peak is at γ_ces≈−9
   (0.2482; the −9→−12 drop is z≈1.4, with Leontief 0.223 continuing the decline).
   Refinement data: `sigma_tierB_refine_*` in `sweep_results/sigma_tierB_concentration/`.

The Story-5 framing (the concentration dispute as a disagreement about σ) is killed *as
stated*. The replacement pattern is coherent: **harder complementarity concentrates
funding only when the budget is tight; with an ample budget it equalizes** (top-up-to-the-
bottleneck logic covers everyone). Boundary: the sign flips between b=0.1 and b=0.5–1.0,
tail-dependent.

**P-B3 — REFUTED as stated, same structure.** Coverage is monotone in σ in P-B3's
direction only at tight budgets (α_K=1.3, b=0.1: 0.19 → 0.12, Leontief 0.10); at ample
budgets coverage rises toward Leontief (α_K=2, b=1.0: 0.964 → 0.977, Leontief 0.98+;
α_K=3.5, b=0.1: 0.81 → 0.86). Mirror of P-B2.

Leontief throughout carries the degeneracy caveat (allocation at the kink; ties broken by
fill order; n_steps=800).

---

## Package C — back-loading attribution

Code: growth decomposed as K_{t+1} = K_t + ε_free·λ(K,R0) + ε_paid·(λ(K,R0+g) − λ(K,R0)),
threaded through the realized dynamics AND all T-round planner anticipation paths (8 call
sites; the funder knows the model). The v5 two-round path retains coupled growth and is
untouched (defaults only). Equal rates take the exact historical code path.

**P-C3 — CONFIRMED (validation identity).** Explicit ε_free=ε_paid=0.3 reproduces the
plain run bit-identically (3 seeds, T=5, output and b_idx). At the sweep level the
diagonal reproduces the coupled model's b_idx to 4 decimals (0.5537 at 0.3; 0.6192 at 0.85
vs `horizon_growth` 0.5537/0.6192).

**P-C1 — CONFIRMED.** Pure free growth (ε_paid=0) back-loads, increasingly:
b_idx_S8 = 0.529 / 0.587 / **0.680** at ε_free = 0.1 / 0.3 / 0.85; PG up to 3.97.

**P-C2 — CONFIRMED.** Pure paid growth (ε_free=0) FRONT-loads: b_idx_S8 = 0.489 / 0.481 /
0.481 at ε_paid = 0.1 / 0.3 / 0.85 — strictly below 0.5 (SE ≤ 0.003), with round-1 share
0.205–0.233 (> 1/5). It is mildly valuable at high ε_paid (PG = +0.19 at 0.85). The
report's stated mechanism is confirmed: **back-loading is driven by the free channel
(wait for researchers to become more productive per dollar); the paid channel alone
argues for early spending (compound over more remaining rounds), but is dominated
whenever both are active** (full 4×4 matrix in `sweep_results/bload_decouple/`).

**Contour transects (refinement rule): the free force dominates at remarkably small
relative strength.** Two 5-point transects at 200 seeds bracket the b_idx=0.5 contour to
within the required factor 2:
- ε_paid=0.30: crossing at ε_free ∈ (0.05, 0.10), interpolated ε_free* ≈ 0.08 —
  ratio ε_free*/ε_paid ≈ **0.27**.
- ε_paid=0.85: crossing at ε_free ∈ (0.10, 0.15), interpolated ε_free* ≈ 0.11 —
  ratio ≈ **0.13**.
Back-loading takes over once free growth is roughly **1/8 to 1/3 the strength of paid
growth** (the ratio falling as ε_paid rises). Data: `bload_transect_*` in
`sweep_results/bload_decouple/`.

**C3 (optional corroboration, decoupled purse) — DONE.** `resource_regime` at 200 seeds
(`sweep_results/resource_regime/`) reproduces the 64-seed supplement run to within
|Δ b_idx| ≤ 0.0011 in every cell: 10 front-load cells, all at ε ≤ 0.03; back-loading at
ε=0.85 runs 0.521 (r_min=0.001) → 0.652 (r_min=3). Same mechanism seen from the resource
side: with the coupled ε, even a zero-resource community back-loads because the funder's
own paid growth resupplies the free channel — exactly the C2 attribution (paid-only
front-loads; a whiff of free growth, ratio ≥ ~1/8..1/3, restores back-loading).

---

## Flags for the paper session (carried from the handoff + new)

1. Model write-up Eq (10) (per-round, non-transferable budgets) contradicts the code and
   report for the forward strategies — spec must be fixed before drafting (handoff flag).
2. Default-parameter tables disagree between write-up and code base params (handoff flag).
3. `RESULTS.md` needs regeneration/deprecation notes after these runs (handoff flag).
4. **New:** Cobb-Douglas eliminates planning value and back-loading (P-B1 exception).
   The paper's Part V claims should be scoped to complementarity production (γ_ces<0);
   this is arguably a headline robustness result in itself.
5. **New:** the concentration–substitutability relationship is budget-conditional
   (P-B2/P-B3 refutations) — a candidate replacement for the Story-5 framing.
6. Implementation note: the ε_free/ε_paid decomposition covers the T-round engine only;
   the v5 two-round anchor path supports equal rates only (by design, defaults).
7. Tier A output layout: one directory per γ_ces value (`sigma_tierA_gc0/`, `_gcm3/`,
   `_leontief/`), four (one) sweeps inside each — trivial deviation from the handoff's
   per-sweep naming, chosen so `sweep_one_T` checkpointing works unchanged.

## Completion checklist (handoff)

- [x] Pre-edit pinned values captured (`tests/pinned_T5_baseline.rds`); all gates green
      after each of the three code changes.
- [x] A: D3 (32 cells), D4 (10 cells), D1 + A5 analyses; κ computed; P-A2 exact.
- [x] B: Tier A four sweeps × {0, −3} + Leontief signal_value; Tier B 72 cells +
      Leontief endpoint + 4-cell refinement at 400 seeds.
- [x] C: decoupling edit validated (P-C3 exact); 16 cells + both contour transects done;
      resource_regime C3 done (matches supplement run to |Δ| ≤ 0.0011).
- [x] RUN_INFO.txt in every output directory; this memo.
- [x] Final commit + push.
