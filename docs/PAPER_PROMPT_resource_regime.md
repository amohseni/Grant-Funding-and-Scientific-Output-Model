# Prompt for the paper-writing project — the resource-poverty / bootstrap results

Paste everything below the line into the writing project.

---

We have a new set of results to integrate into the paper (the T-round grant-funding model,
smooth/exact allocator). They answer the most natural objection to the Part V back-loading result:
**"Surely a funder facing a resource-poor community should front-load — invest early to bootstrap
growth and learning?"** The answer is no, and the reasons produce two new mechanisms worth a
subsection. Everything below is self-contained; all numbers are from 64-seed sweeps (T-round model,
exact allocator). Where I give a metric: **b_idx** = budget-schedule center of mass (0.5 = even
split, <0.5 = front-loaded, >0.5 = back-loaded); **S8−S5** = forward vs myopic planning value;
**S8−S2** = optimal vs uniform (the value of discrimination); shares compare to the even split 1/T.

## One modeling addition (must be stated in the setup)

In the base model the budget is normalized to community resources, B = b·n·E[R] — so a poor
community mechanically has a poor funder, and r_min→0 is a trivial null. We added a decoupled
normalization B = b·n·E[K] (`budget_ref="K"`): the funder keeps a fixed purse while baseline
resources vary. At default parameters E[R]=E[K], so this changes nothing in any existing result;
it only opens the poverty regime to study.

## The motivating conjecture (present it, then correct it)

At R₀=0: no funding → no output → no learning (publications are the funder's signal) and no
knowledge growth (K ← K+ελ). So money in round T buys nothing for T−1 rounds, while money in
round 1 enables output, learning, and growth for the whole rest of the horizon — round-1-heavy
funding "must dominate."

**The gap:** at R₀=0 the two lump schedules are *payoff-identical*. All-in-round-1 also allocates
blind (spending precedes observation), and the learning and growth it generates are worthless with
no later money to exploit them. What the argument really establishes is a **complementarity** —
early money *enables*, later money *exploits* — which guarantees strictly positive round-1 spend
and rules out all-late schedules, but does not put the *mass* early. The information logic in fact
cuts the other way: information's value accrues to dollars spent *after* it arrives. **Money
follows information, and information arrives late.**

## Result 1 — the front↔back reversal is governed by ε, not by resources

Sweep `resource_regime`: r_min ∈ {0.001…3} × ε ∈ {1e-4…0.85}, T=5, fixed purse. The b_idx = 0.5
boundary is essentially **vertical at ε\* ≈ 0.02** (range 0.005–0.05 across r_min). Poverty *mutes*
back-loading — at ε=0.85, b_idx runs 0.521 (r_min=0.001) → 0.651 (r_min=3), and the round-1 share
runs 0.174 → 0.059 — so poorer communities do get relatively more early money (the conjecture's
comparative static is right). But the sign never flips: even at r_min=0.001 the funder's own early
grants raise recipients' K, manufacturing the compounding that rewards spending late. Front-loading
exists only at ε ≲ 0.02 and is worthless there (S8−S5 ≈ 0). **Back-loading is over-determined: the
paid bootstrap resupplies the free-growth condition that drives it.**

## Result 2 — thin grants are talent-uninformative (new mechanism)

With λ = γKg/(K+g), a grant g ≪ K pins output to the grant, not the talent: λ ≈ g (at g=1/3, a K=1
researcher yields λ=0.25 and a K=100 researcher λ=0.33 — indistinguishable). In the pure-exploration
corner (no peer-review signal, poverty, ε≈0, heavy talent tail α_K=1.3, T=6), at standard depth
(b=0.5) the funder's discrimination value S8−S2 is **+0.9** — versus **+17.5** when a sharp free
signal is available. Its schedule is flat and planning is worthless. **In poor fields, paying for
information is not "fund them" — it is "fund them deeply enough that output becomes talent-limited":
g ≳ K, i.e. budget scale b ≳ T/2.**

## Result 3 — depth makes paid information work, and it back-loads the money

Sweep `exploration_depth` (b ∈ {0.5, 1.5, 3, 6} × signal on/off, same corner): at b=3 the
discrimination value S8−S2 jumps from 0.9 to **39** (b=6: **56**) — funded output now reveals
talent. The optimal schedule becomes **seed-and-harvest**: round-1 share **0.107** at b=3 (even =
0.167), with the mass deployed from round 2 on, once the funder knows who is who. b_idx *rises* to
0.51–0.53; across all 54 cells of the suite it never falls below 0.5 beyond noise. **Learning
front-loads the observation and back-loads the money.**

## Result 4 — the free forces dominate even when switched off

At depth, the forward planner's deliberate schedule *loses* to the myopic even-tranche funder
(S8−S5 = −16 at b=3, −45 at b=6): the myopic funder re-decides each round on posteriors updated by
its own funding's output, capturing the paid information automatically. This completes a symmetry
with the growth story and is the five-force scaffold's cleanest vindication: **paid growth (B)
regenerates free growth (C); paid observation (D) regenerates free observation (E). The free forces
dominate even when their exogenous sources are switched off, because the paid forces resupply
them.** Caveat to state: the CE forward planner misprices information at depth, so its schedule is
not a certified optimum; the claim rests on flat-tranches-plus-re-deciding (S5) beating deliberate
deferral (S8), and on no force favoring early mass (statics favor even; information favors late).

## Suggested placement and framing

A subsection extending Part V — e.g. **"The bootstrap objection: poor communities and paid
information"** — structured as: (1) the conjecture, stated sympathetically; (2) the lump-equivalence
correction; (3) Results 1–4; (4) the closing symmetry (B→C, D→E) tying back to Part II's scaffold.
It strengthens the paper's central thesis rather than qualifying it.

## Figures available (figures/resource_regime/)

1. `resource_regime_heatmap.png` — b_idx over r_min × ε; near-vertical boundary at ε\*≈0.02.
2. `resource_regime_bidx_vs_eps.png` — b_idx vs ε (log) by r_min; curves fan out above ε\*.
3. `exploration_depth_schedules.png` — full round-by-round schedules at depth: the seed-and-harvest
   shape (below-even round 1, mass once informed). Probably the single best figure.
4. `exploration_depth_vsunif.png` — S8−S2 vs depth, signal on/off: paid information switching on.
5. `exploration_poverty_r1share.png` — round-1 share vs r_min: the tilt requires learnable output.

## Provenance

Sweeps `resource_regime` (30 cells), `exploration_corner` (16), `exploration_poverty` (8),
`exploration_depth` (8); 64 seeds/cell; manifest entries in `sweep_T.R`; data
`sweep_results/T_run_smooth_supplement/`; full write-up
`T_round_extension/RESOURCE_REGIME_RESULTS.md`; model change `budget_ref` in `model.R`
(`run_simulation_T`), bit-identical to the published model at defaults.
