# The seed-floor puzzle: resolution and diagnostics

2026-08-05. Aydin's objection: "if a seed floor is not impactful, uniformly apply all money and no
strategy should do anything." Verdict: the inference is invalid, the model's numbers cohere, but the
instinct is vindicated at the design level: the decisive experiment has never been run.

## The three claims

- C1 (model result): round-1 floor, re-optimized remainder, no-signal family: costs 0.03-0.40% of S1.
- C2 (model result): full uniformity (S2) vs optimal+signal (S8): costs 9-25 points of S1.
- C3 (Aydin's inference): C1 implies allocation is irrelevant, contradicting C2. INVALID; see below.

## Why C1 and C2 cohere

1. Design facts: the seed floor is round-1 only (25% of total budget at T=2, x_seed=0.5, not 50%);
   the remainder re-optimizes around it; S6/S9 carry NO grant signal, so the floor has only ever
   diluted the weakest optimizer in the set.
2. The right yardstick is targeting value within the no-signal family: S4-S2 = 1.6 pts (base),
   3.1 (heavy tail). Linear dilution bound at x_seed=0.5: 0.25 x 1.6 = 0.40. Observed: 0.24.
   At x_seed=0.75: bound 0.60, observed 0.36-0.40. Below bound but of its order. Nothing anomalous.
3. Mechanisms for the sub-linear part:
   (a) Inframarginal absorption: the optimum is a top-up rule; for any researcher with g* >= floor,
   re-optimization restores the identical final allocation. Only dollars to g* < floor researchers
   are misallocated, and they still earn that researcher's positive marginal product; the loss is the
   wedge nu - MP_i, not the dollar.
   (b) Envelope logic: near a concave optimum, small feasible-set perturbations cost second-order;
   abandoning optimization costs first-order. Locally flat, globally steep. C1 and C2 are two ends of
   one curve; S2 is the x_seed=1, all-rounds, no-reoptimization endpoint and the model prices it.
4. Aydin's "small fraction zeroes all gaps" hypothesis, corrected form: the budget never zeroes gaps
   (c adjusts to exhaust B); rather, most floor dollars are inframarginal. Directly computable (D2).

## Where the instinct is right

"Bayesian optimization plus dilution" is unpriced: no seed variant of the SIGNAL strategies exists.
And seed_value ran only at base tails. If cost scales with the diluted family's targeting value, the
seed result is the targeting result read backwards, and one quantity (dispersion of marginal returns
= the value of knowing the gap) prices peer review, targeting, and floors alike. Thesis-reinforcing.

## Diagnostics, ranked

- D1 (free): convexity of cost in x_seed from existing seed_value summary. Predict superlinear.
- D2 (one detail=TRUE run): inframarginal fraction of floor dollars. Predict >2/3 at base.
- D3 (DECISIVE, small rerun): S5+seed and S8+seed on the b x x_seed grid at alpha_K in {2, 1.3}.
  Predict: cost scales with family targeting value; several points at heavy tail + signal.
  If instead it stays <0.5 pts, "floors nearly free" is robust and both parties update.
- D4 (small rerun): persistent every-round floor, x_seed -> 1, no-signal family. Predict cost
  interpolates to S4-S2. Closes the C1<->C2 circuit.
- D5: seed_value at alpha_K=1.3, tight b.

## Paper consequence

Hold the seed-floor section until D3/D4. Target claim, if D3 comes back proportional:
"A uniform floor is cheap exactly where targeting is cheap; where peer review earns its keep, floors
are costly in proportion." State the section as a price schedule, not a verdict; engage the lottery
literature where the price is lowest.

Rerun notes: smooth myopic waterfill is harmonic-only (fine here); D3/D4 need a seed flag on signal
strategies + a persistence flag, small edits to the strategy table in model.R, no new planner.
