# Resource-poverty and the funding schedule: does a poor community flip the planner to front-loading?

**Question (motivating hypothesis).** The paper's Part V result is that the forward planner *back-loads*
(spends more in later rounds) because knowledge compounds from baseline resources anyway (free force C),
so later rounds face higher-K researchers whose marginal product is higher. Conjecture: in a
**resource-poor community** the free forces switch off — with baseline resources `R₀≈0` there is no output,
so no free knowledge growth (C) and no free observation of talent (E). The only way to generate *any*
growth or *any* information is to **pay** for it (forces B and D). So the forward-looking Bayesian planner
should **front-load** — invest early to bootstrap growth and learning — and this should **reverse** to
back-loading as the community's own resources rise. Where is that reversal?

**Metric.** `b_idx_S8` = center-of-mass of the forward planner's budget schedule, `Σ t·αₜ / (T+1)`.
`0.5` = even split; **`<0.5` = front-load** (spend early); **`>0.5` = back-load** (spend late).

---

## Setup

Posing the question requires a **decoupled budget**. Under the model's default normalization
`B = b·n·E[R]` the funder's purse is tied to baseline resources, so the injection-to-baseline ratio is
`2b/T`, *independent of `r_min`*; driving `r_min→0` also zeroes the budget — a trivial null. We added a
back-compatible switch `budget_ref` to `run_simulation_T`: `"K"` sets `B = b·n·E[K]`, holding the funder's
purse **fixed** while baseline resources `r_min` vary. (At base params `E[R]=E[K]=2`, so `"K"` is
bit-identical to the default `"R"` at `r_min=1`; they diverge only when `r_min≠1`.)

Sweep **`resource_regime`** (added to `sweep_T.R`, tier 1): `r_min ∈ {0.001, 0.03, 0.3, 1, 3}` ×
`epsilon ∈ {1e-4, 0.01, 0.03, 0.1, 0.3, 0.85}`, at `T=5`, `budget_ref="K"`, smooth (exact) allocator,
base `tau_k=1, k_shape=2`. 30 cells × 64 seeds. Data: `sweep_results/T_run_smooth_supplement/`;
figures: `figures/resource_regime/`.

---

## What the data say

**1. The reversal is real, but it is governed by the compounding rate ε, *not* by baseline resources.**
The `b_idx=0.5` boundary is essentially **vertical** in (ε, r_min) space, at **ε\* ≈ 0.02** (range
~0.005–0.05 across r_min). Below ε\*: front-load; above: back-load, rising steeply with ε.

**2. Baseline resources shift the schedule's *strength*, not its *direction*.** Scarcity does push
spending forward — the user's directional intuition is correct — but not enough to cross into strict
front-loading at any realistic ε. At `ε=0.85`, `b_idx_S8` runs 0.521 (`r_min=0.001`) → 0.651 (`r_min=3`),
and the round-1 share falls from 0.174 → 0.059. Richer baseline ⇒ more free compounding (force C) ⇒ later
rounds more productive ⇒ **back-load harder**. This quantifies force C.

**3. Where front-loading occurs, it is economically negligible — foresight buys nothing there.** In
*every* front-load cell (`ε ≲ 0.01–0.03`), the planning value `fwd_vs_myo_PG = S8−S5 ≤ 0` (zero within
noise, occasionally slightly negative), and the round-1 share is within ~0.005 of the even split (0.20).
With no compounding there is nothing for a schedule to exploit, so front-vs-back is indeterminate and
myopic ties the planner.

**4. Why the simple hypothesis fails.** Even at `r_min≈0`, the funder's *own* early (paid) grants raise
those researchers' K, recreating exactly the higher-K-later condition that rewards back-loading. The paid
bootstrap manufactures its own compounding. So **paid knowledge (force B) robustly dominates paid
information (force D) in setting the schedule**, and both, once compounding is present, point at
back-loading. There is **no regime in which forward-looking front-loading is both strict and valuable.**

**5. Secondary — information value is non-monotone in baseline resources.** `signal_fwd = S8−S7` peaks at
intermediate `r_min≈1` (≈10.8) and is lower at both extremes (≈2.3 at `r_min=0.001`, ≈6.6 at `r_min=3`):
too poor ⇒ little output to observe; too rich ⇒ free observation (E) already reveals talent, so the paid
signal adds less.

---

## Interpretation for the paper

This *sharpens* Part V rather than overturning it. The clean statement is:

> Back-loading is **over-determined**. It arises whenever knowledge compounds (ε>0), whether the
> compounding is free (from baseline resources) or paid (from the funder's own early grants). Resource
> poverty removes the *free* channel and mutes back-loading toward an even split, and it does push
> spending marginally forward — but it cannot flip the planner into strict front-loading, because the
> planner's own bootstrapping resupplies the compounding. The only front-loading force, paid information
> (D), is too weak to reverse the schedule and is worthless precisely where it would act (ε≈0).

The reversal the conjecture sought exists only as a **degenerate ε→0 edge**, not as a resource-driven
regime. The paper can state the boundary crisply: **the forward planner front-loads only when ε ≲ 0.02,
and only trivially; for any ε above that it back-loads at every level of community resources.**

---

## Part 2 — The dominance argument and the pure-exploration corner

**The argument.** With `R₀=0` for everyone: no funding → no output → no learning and no growth. So
funds given only in round T buy nothing for T−1 rounds, while funds in round 1 enable output,
learning, and growth for all later rounds — "this has to dominate logically," suggesting a regime
where round-1-heavy allocation is optimal.

**Where the logic slips.** At `R₀=0` the two *lump* schedules are payoff-identical: all-in-round-1
also allocates blind (spending precedes observation), and the learning and growth it buys are
worthless with no later money to exploit them. What the argument really establishes is a
**complementarity** — early money *enables*, later money *exploits* — which guarantees (i) strictly
positive round-1 spend at the optimum and (ii) interior schedules dominating both lumps. It does
*not* determine where the *mass* sits. And the information logic cuts the other way: information's
value accrues to the dollars spent *after* it arrives. Money follows information; information
arrives late.

**Three sweeps** (`exploration_corner`, `exploration_poverty`, `exploration_depth`; T=6, heavy tail
`k_shape=1.3`, decoupled purse, smooth allocator, 64 seeds; S2 uniform included as the
no-discrimination benchmark) test the corner the argument points to: poverty (`r_min=0.001`, initial
pubs ≈0), free signal off (`tau_k=100`, uninformative), ε≈0 vs ε>0, and — decisively — funding
**depth** `b ∈ {0.5, 1.5, 3, 6}`.

**Finding 1 — thin grants make funded output talent-uninformative.** With `λ = K·g/(K+g)` and
per-capita grants `g ≪ K`, output is resource-limited: `λ ≈ g` regardless of talent (K=1 → λ=0.25;
K=100 → λ=0.33 at g=1/3). At standard depth (b=0.5) in the no-signal poverty corner, the funder
learns nothing from its own funding: S8 beats uniform by only ~0.9 (vs +17.5 with a sharp free
signal), the schedule is flat (b_idx≈0.50), and planning value is ≈0. **Paying for information is
not "fund them" — it is "fund them deeply enough that output becomes talent-limited (g ≳ K, i.e.
b ≳ T/2)."**

**Finding 2 — depth makes paid information work, and it back-loads the money.** Raising b from 0.5
to 3 in the no-signal corner lifts the discrimination value S8−S2 from 0.9 to **39** (b=6: 56) —
funded output now reveals talent. The optimal schedule becomes **seed-and-harvest**: round-1 share
falls to **0.107** (b=3; even = 0.167) with the mass deployed from round 2 on, once informed —
`b_idx` *rises* to 0.51–0.53. Nowhere in any sweep does `b_idx_S8` fall below 0.5 by more than
noise. **The hypothesized front-loading regime does not exist in this model: learning front-loads
the *observation*, and back-loads the *money*.**

**Finding 3 — even paid information is captured "for free" by re-deciding.** In the depth cells the
forward planner's deliberate schedule *loses* to the myopic even-tranche funder (S8−S5 = −16 at b=3,
−45 at b=6): the myopic funder re-decides each round on posteriors updated by its own funding's
output, capturing the information without planning for it. This extends the paper's central pattern
to the information channel under poverty: just as the paid knowledge bootstrap (B) regenerates the
free-growth condition (C) that favors back-loading, paid observation (D) regenerates the free-
observation condition (E) that the myopic funder exploits automatically. **The free forces dominate
even when their exogenous sources are switched off, because the paid forces resupply them.**

(Caveat: S8's schedule is the CE planner's choice, not a certified optimum — the negative S8−S5 at
depth shows the CE approximation misprices information there. But that same fact means flat-even-
with-re-deciding is already ≥ the planner's deferred schedule, and no channel favors early mass:
statics favor even, information favors late. The conclusion stands on S5 vs S2: +55.7 at b=6.)

Figures: `figures/resource_regime/exploration_depth_{bidx,vsunif,schedules}.png`,
`exploration_poverty_r1share.png`.

---

## Reproduce

```r
source("model.R"); source("sweep.R"); source("sweep_T.R")
# full manifest-quality run (200 seeds, all strategies):
sweep_one_T("resource_regime", seeds = 1:200,
            base_params = list(allocator = "smooth"),
            out_dir = "sweep_results/T_run_smooth_supplement")
# figures + boundary table:
source("sweep_results/_probe/plot_resource_regime.R")
```

The results above use 64 seeds and focal strategies {5,7,8} (the contrasts `b_idx_S8`, `S8−S5`, `S8−S7`);
the sign/shape of every claim is stable across the probes at 12–16 seeds and 120–150 posterior samples.

**Provenance upgrades (2026-08-05/06):** `resource_regime` re-run at 200 seeds
(`sweep_results/resource_regime/`, |Δ b_idx| ≤ 0.0011 vs the 64-seed run) and all three
exploration sweeps re-run at 200 seeds (`sweep_results/exploration_200/`, |Δ| ≤ 0.005) — use these
for main-text figures. Two independent verifications (docs/DIAGNOSTICS_RESULTS_2026-08.md, bootstrap
items): (1) an honest fixed-schedule sweep with no CE planner certifies that no early-mass schedule
beats even-or-late at the depth corner (best early-mass loses by z=4.9 at ε=0.3) — seed-and-harvest
is planner-free; (2) the negative S8−S5 at depth is honest CE mispricing, not an optimizer bug
(S8 maximizes its own objective in 20/20 seeds). Part 2's caveat stands as written.
