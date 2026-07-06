# Handoff — replace the greedy allocator with a smooth optimizer

*Context transfer for continuing the "smooth allocator" work in a fresh session. Written 2026-07-05.
Everything you need to finish is here or linked from here. Read this top-to-bottom once before touching code.*

---

## 0. TL;DR — where we are in one paragraph

The T-round study was **finished and packaged** (see `T_round_extension/STATE_OF_PLAY.md`, `ROADMAP.md`,
`RESULTS.md` — all 16 sweeps run, figures made, docs written). Then a diligence pass found a small,
isolated **spike in the forward-vs-myopic contrast at b=0.3, T=2, ε=0.1** (`fwd_vs_myo_PG` ≈ z24,
localized to `fwd_vs_myo_P` = S7−S4). We **diagnosed it as a greedy-discretization *artifact*, not a
bug** (proof below), and decided (**"Approach A"**) to **replace the discrete greedy allocator with a
smooth continuous optimizer** — water-filling for the single-round problem, projected-gradient for the
forward planner — which removes the discretization step δ entirely and is more faithful. **Phases 1 and
2 are both built, validated, and integrated into `model.R` behind an `allocator="smooth"` switch. The
decisive go/no-go test PASSED (GO).** The `"greedy"` default remains bit-identical, so nothing published
has changed yet. **Next task: flip the default to `"smooth"` and re-run the 16-sweep manifest + figures +
prose (§8)** — the "re-run everything" the user signed up for.

---

## 1. The problem: what the "resonance" is

The allocator, `greedy_round_t()` (model.R:943), fills a round's budget in **equal discrete steps**
`δ = B / n_steps` (default `n_steps = 50`), each step handed to the researcher with the current highest
marginal `∂λ/∂g`. Because the Pareto knowledge draw creates a *ladder* of researchers with similar
marginals, at **special budget values the δ-grid lands differently for the forward (S7) vs myopic (S4)
schedules**, opening a gap between two allocations that should be nearly identical. That gap is the
"resonance." It is an **artifact of the step size δ**, not a real property of the model: refine δ
(raise `n_steps`) and it vanishes.

This is the *second* granularity artifact found in this project. The first (`horizon_long`, T>5) is
documented in `STATE_OF_PLAY.md` §4a and was fixed by pinning `n_steps=400` there. This b=0.3 one is
the same disease at T=2. (Distinct from — and downstream of — the earlier **SIR bug** in
`ce_reweight_posterior`, `STATE_OF_PLAY.md` §5, which was a genuine bug, already fixed.)

## 2. Proof it's an artifact, not a bug (the diagnostics)

Script: `tests/diag_b03_resonance.R` · full output: `tests/diag_b03_resonance.log`. A **correct** effect
converges as δ→0; an **artifact** decays to 0. It decays to 0:

- **Test 1** (ε=0, T=2, b=0.3, pubs-only; at ε=0 the continuous optimum forces S7−S4≡0):
  `fwd_vs_myo_P` = **+0.111 (z10.5)** at `n_steps=50` → **+0.058** (100) → **−0.004** (200) → **flat ~0**
  through 1600. Monotone decay to zero = discretization artifact. ✓
- **Test 2** (ε=0.1, fine b-grid): a *comb* of resonances at `n_steps=50` (spikes at b=0.24, 0.30, 0.32:
  +0.10 to +0.14), **all flat (≤0.01)** at `n_steps=800`. It's a grid phenomenon, not physics. ✓
- **Test 3** (ε=0, L1 allocation distance ‖g1_S7 − g1_S4‖₁/B): small everywhere (~0.02–0.03), i.e. the
  schedules differ by only a few δ-steps — consistent with a discretization gap.

**Verdict: confirmed greedy-δ resonance. Proceed to replace the allocator (Approach A).**

## 3. The decision and the plan (Approach A)

User's decision, verbatim: *"Let's try approach A. If it works, wonderful. We'll set ourselves up to
re-run all the parameter sweeps. If it doesn't work, maybe we can do the more granular end steps."*
Agreed de-risking: **build the single-round water-filling first and prove it's exact + removes the
resonance, before touching the forward planner.**

- **Approach A** = replace greedy with a smooth optimizer (no δ). This is the chosen path.
- **Fallback** (if A's Phase 2 proves too hard) = just raise `n_steps` on the affected sweeps (the
  "more granular end steps"), as was already done for `horizon_long`. Cheaper, less faithful, but safe.

---

## 4. STATUS TABLE

| Step | Status | Artifact |
|---|---|---|
| Part 1 — diagnose the resonance | ✅ **Done**, artifact confirmed | `tests/diag_b03_resonance.R` + `.log` |
| **Phase 1 — single-round water-filling** | ✅ **Built & validated** (exact, fast) | `tests/waterfill_round_t.R` |
| **Phase 2 — forward continuous optimizer** | ✅ **Built & validated** (dominates greedy, KKT ✓) | `tests/plan_forward_smooth.R` |
| Decisive go/no-go test (ε=0 ⇒ S7−S4→0 at low `n_steps`) | ✅ **PASS → GO** (100 seeds) | `tests/gonogo_smooth.R` + `.log` |
| Integrate into `model.R` (myopic + forward paths) | ✅ **Done** — behind `allocator="smooth"` switch | §7 (model.R:1061 block) |
| Re-run 16 sweeps · regen figures · update docs | ⏳ **NEXT** — flip default + re-run | §8 |

**Go/no-go result (2026-07-05, `tests/gonogo_smooth.log`):** A) at ε=0, `max|smooth S7−S4| = 0.0056`
(all |z|≤1.4 ≈ 0) vs the greedy-ns50 comb (+0.096/+0.121/+0.142 at b=0.24/0.30/0.32). B) at ε=0.1,
`max|smooth − greedy800| = 0.0127` vs greedy50-vs-greedy800 = 0.1507 (12× closer). **Approach A works.**
At ε=0 greedy-ns800 still shows a faint positive lean at the resonance b's while smooth is dead-flat — the
exact optimizer is if anything *more* correct than fine greedy.

**Reversibility: still total.** The `allocator="greedy"` default is **bit-identical** to before
(T2 anchor PASS, max |Δ| = 4e-14); `sweep_results/T_run_fixed/`, `T_round_extension/data/`, figures, and
`RESULTS.md` are all **untouched**. To ship Approach A, flip the default to `"smooth"` and run §8; to
abandon it, delete the new `tests/` files and the model.R smooth block.

---

## 5. Phase 1 — single-round water-filling (DONE, validated)

**File:** `tests/waterfill_round_t.R` (function `waterfill_round_t()` + a self-checking harness). Run
`Rscript tests/waterfill_round_t.R` to re-verify at any time. Latest run:

```
water-fill: obj=46.9173483 spend=30.0000000 time=0.0205s/call
vs fine greedy ns=3200: obj=46.9173300  wf-greedy=+1.83e-05 (>=0)  ||dg||1/B=0.0038
water-fill = 6.1x greedy-ns50 (0.00335s), 0.3x greedy-ns800 (0.063s)
```

- **Exact:** it is the analytic KKT optimum — matches the fine-greedy (`n_steps=3200`) limit to 1.8e-5,
  **dominates** greedy at every `n_steps`, spends the budget exactly, and every funded researcher sits
  at a common marginal (water level). The +1.8e-5 gap is greedy still being slightly *sub*-optimal.
- **Fast:** 0.02 s/call — **2.5× faster than the *correct* alternative** (`n_steps=800` greedy, 0.063 s).
  (It's ~6× the *buggy* `n_steps=50` greedy, but that baseline is the thing we're removing.)
- **Faithful to the objective (critical subtlety):** the model's round-t objective is the
  **posterior-EXPECTED** rate `Σ_m w_m·λ(K_m, R_m+g)` averaged over the M importance-sampling atoms —
  **NOT** λ evaluated at the posterior mean. So water-filling needs the **mixture-aware** condition
  `Σ_m w_m·γK_m²/(K_m+R_m+g)² = ν` (a nested solve per researcher), which the prototype implements.
  A single-atom closed form `g_i(ν)=K√(γ/ν)−K−R` would be wrong. Verify this if you ever simplify it.
- **Speed history** (for context, don't repeat): naive nested bisection 3.4 s → vectorized 0.80 s →
  inner-Newton 0.65 s → **outer safeguarded-Newton + warm-started inner Newton 0.02 s** (current). The
  cost was the *count* of n×M matrix passes (outer × inner), not flops; Newton on both loops fixed it.

## 6. Phase 2 — the forward continuous optimizer (NEXT TASK, the crux)

The forward CE objective is **smooth** in the grant matrix G (H rounds × n), and its per-slot gradient
**already exists**: `fwd_marginal()` (model.R:994) is exactly ∂(CE objective)/∂g_{s,i}. So keep the
entire planner — receding horizon, `ce_reweight_posterior` information anticipation, K-compounding — and
**replace only the inner loop of `plan_forward_ce()`** (model.R:1017): swap the "add δ to the global
argmax slot" `while` loop for **projected-gradient ascent on the budget simplex** {ΣG = B_rem, G ≥ 0},
using `fwd_marginal` as the gradient and recomputing `bar1 = ce_reweight_posterior(π0, g₁)` as g₁ moves
(it already does this incrementally). Euclidean simplex projection is standard O(Hn log Hn).

- **Why it should work at ε=0** (the target): no compounding ⇒ rounds couple only through the budget;
  the water-fill value per round is concave in that round's budget, so the optimal split is uniform
  (Jensen) — *identical to the myopic equal-tranche + water-fill*. Hence S7−S4 → 0 **without** large
  `n_steps`. With ε>0, compounding makes early rounds worth more ⇒ the optimizer front-loads (the
  `b_idx > 0.5` story in `RESULTS.md`), reproducing the real forward advantage exactly.
- **Risk:** the `ce_reweight` information term may make the objective **non-concave**. If projected
  gradient stalls or is path-dependent, fall back to **Frank-Wolfe** (linear oracle = all remaining
  budget to the single top-marginal slot; no projection needed) or multi-start. Flagged in the original
  prescription too.
- **Reuse, don't rewrite:** `fwd_researcher_value()` (objective, model.R:972) and `fwd_marginal()`
  (gradient, model.R:994) are the pieces you need; the water-fill from Phase 1 is the natural
  within-round subroutine.

**The decisive go/no-go test** (run this the moment Phase 2 compiles): with S4 = water-fill and
S7 = continuous-forward, at **ε=0** the contrast `fwd_vs_myo_P` (S7−S4) must be ~0 at **`n_steps`-free /
low granularity** across all b — i.e. reproduce Test 1's *converged* row without needing `n_steps=800`.
AND at ε>0 the results must match the fine-greedy limit (so the published science is preserved). If both
hold, Approach A works; integrate. If not, use the fallback (§3).

## 7. Integration map (once Phase 2 passes)

Exact edit sites in `model.R`:

1. **Myopic path** — `allocate_round()` (model.R:1067): the three `greedy_round_t(...)` calls for
   **S4/S5 (line 1079)** and **S6 (lines 1083, 1086)** become `waterfill_round_t(...)` (drop the
   `delta` arg; keep `g_hist`, `init_g` for the S6 seed floor). Add `waterfill_round_t` into `model.R`.
2. **Forward path** — `plan_forward_ce()` (model.R:1017): replace the `while (remaining >= delta)`
   greedy loop (lines 1043–1057) with the projected-gradient / Frank-Wolfe solve. Keep the `H == 1`
   branch (lines 1021–1028) as a single-round water-fill. Keep `bar1` recomputation.
3. **δ becomes vestigial:** `delta`/`n_steps` (model.R:626, 641, 1209, 1219) can stay as no-ops for
   signature stability, or be removed. `horizon_long`'s `n_steps=400` pin becomes irrelevant.

**Consequence for validation:** the **T=2 ⇄ v5 bit-identical anchor breaks by design** (v5 uses greedy).
It becomes "agree to solver tolerance." Update `tests/test_T2_reduction.R` expectations accordingly and
note it in `STATE_OF_PLAY.md` §3. This is expected, not a regression.

## 8. Downstream consequences (the "re-run everything" the user signed up for)

Once integrated, the allocator change touches every forward/myopic number:

- **Re-run the full 16-sweep manifest** → a NEW output dir (do **not** overwrite `T_run_fixed/`; keep it
  for provenance). Command pattern: `source("model.R"); source("sweep.R"); source("sweep_T.R");
  main_sweep_T(seeds = 1:200, out_dir = "sweep_results/T_run_smooth")`. Budget ~1–1.5 h at Phase-1 speed
  for the myopic side; Phase-2 forward cost TBD (the forward planner does more solves — benchmark first).
- **Regenerate all figures** (`figures/make_figures.R`) from the new `.rds`.
- **Diff every summary cell** that moved > 2 SE vs `T_run_fixed/`; the headline story should *survive*
  (the resonance was tiny and local) but confirm it, don't assume it.
- **Update the prose:** `RESULTS.md`, `T_round_extension/STATE_OF_PLAY.md` (§4a — the b=0.3 note; §5 —
  add the allocator change), `DATA_DICTIONARY.md`, and add a **one-paragraph methods note** on the
  water-filling / projected-gradient allocator + why (faithfulness, removes δ artifacts).
- **Re-point** `T_round_extension/data/` to the new canonical run; move `T_run_fixed/` to `legacy/`
  labeled "greedy allocator, superseded."

## 9. Key model facts you must not re-derive wrong

- **Objective = posterior-EXPECTED** over M atoms, not λ at posterior mean (see §5). This is *the* easy
  mistake; `post_lambda_round_t` (model.R:923) is the ground truth: `sum(post$w * lambda_rate(K, R0+g))`.
- **Production:** λ(K,R) = γKR/(K+R). **Knowledge compounding:** K ← update_knowledge(K, R0+g, ε).
  **Grants non-persistent:** round-t resources = R0 + g_t only. `update_knowledge`, `lambda_rate`,
  `draw_initial_population`, `draw_signals`, `build_posteriors_hist` are all in `model.R`.
- **Budget:** B_total = 2·b·n·E[R], E[R] = Pareto mean r_min·α_R/(α_R−1); per-round tranche =
  B_total/T for all non-forward strategies; the forward planner allocates B_rem across the horizon.
- **9 strategies:** S1 none, S2 uniform, S3 naive∝pubs, **S4/5/6 myopic** (pubs / +grant-signal / +seed),
  **S7/8/9 forward** (pubs / +grant / +seed). Contrasts: `fwd_vs_myo_P` = S7−S4, `fwd_vs_myo_PG` = S8−S5.
- **Forward = receding-horizon certainty-equivalent** planner: each round plan the full remaining
  horizon (anticipate compounding + one-step info via `ce_reweight_posterior`), execute round 1, re-plan.
- **M = 400** importance-sampling atoms (bumped from 200 during the SIR-bug fix for ESS headroom).

## 10. How to resume (fresh session, from project root)

```r
# 1. See the validated Phase-1 optimizer + re-run its checks:
Rscript tests/waterfill_round_t.R

# 2. Re-confirm the artifact diagnosis (optional, ~a few min on 8 cores):
Rscript tests/diag_b03_resonance.R          # reproduces tests/diag_b03_resonance.log

# 3. Build Phase 2: prototype a continuous plan_forward_ce (projected gradient
#    reusing fwd_marginal), then run the decisive go/no-go test (§6):
#    ε=0  ⇒  S7−S4 ≈ 0 at LOW n_steps across b   (currently needs n_steps=800 with greedy)
#    ε>0  ⇒  matches the fine-greedy limit
# 4. If it passes: integrate per §7, re-run manifest per §8.
# 5. If it fails: fall back to raising n_steps on affected sweeps (§3).
```

## 11. Pointers to the rest of the project (unchanged, still valid)

- `T_round_extension/STATE_OF_PLAY.md` — the finished-study snapshot (bug history, results, caveats).
- `ROADMAP.md` — path to publication; the user's remaining items (CITATION.cff ORCID/co-authors,
  materialize iCloud files, GitHub push, Shiny redeploy, OSF, write paper).
- `RESULTS.md` — the statistical digest the paper is written from (numbers will shift after §8).
- `docs/PAPER_HANDOFF.md` — the paper-writing handoff prompt.
- `model.R` — the engine; all function line numbers above are current as of this writing.

**Bottom line for the next session:** Phase 1 is a clean, banked win. The whole task reduces to **§6 —
build the continuous forward planner and run the ε=0 go/no-go test.** That single result decides whether
Approach A ships (then §7–§8 mechanical) or we take the fallback.
