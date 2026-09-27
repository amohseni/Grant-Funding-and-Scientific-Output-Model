# Section 8 (Seed grants and lotteries): function analysis and preparation (2026-09-08)

Not a draft. The floors half is fully determined by verified results; the lotteries
half has genuine design decisions that are yours (and plausibly Kevin's and Simon's)
before drafting. This file: the section's function, the verified inputs, a proposed
structure, preliminary pricing results, and the decisions I need.

## 1. Function

What the paper has promised \S8, gathered:

- Abstract: "uniform seed grants and lotteries, which forgo targeting, are cheap
  exactly where review is worth little."
- Intro: "\S8 prices seed grants and lotteries"; the paper answers "what egalitarian
  alternatives cost."
- \S2 (the sharpest promise): the existing case for lotteries prices what a lottery
  SAVES (review costs, wasted proposal effort: Gross and Bergstrom's contest model);
  "what remains unpriced is what a lottery forgoes. Our model supplies that half of
  the ledger: the value of targeting, which runs from negligible, where budgets are
  ample, signals uninformative, or capability evenly spread, to substantial, where
  capability is heavy-tailed and budgets tight." And: "when lotteries are defensible."
- \S7 bridge: "Next we price two schemes that recur in policy debates: uniform seed
  grants and lotteries."

So \S8's function is composition, not novelty: convert the machinery already built
(the gap rule and the value of targeting from \S4; review's field-dependent value and
the fine-distinctions corollary from \S6) into prices for the two egalitarian schemes,
field by field, and close the ledger against the savings side the literature has
priced. It is the paper's policy payoff. No new mechanisms; every claim should reduce
to an established result applied to a scheme.

## 2. Verified inputs (all re-derived this session; verify_s8.R, price_lotteries.py)

Seed floors (Package A, D3/D4 summaries, 200 seeds, T=2):
- Focal cost (floor on the strongest funder where targeting matters: heavy tail,
  sharp review, b=0.5, three-quarters floored): 6.1 percent of no-funding output.
- The law: cost ~ (value of targeting) x (floored share) x kappa, kappa in
  [0.10, 0.86]; kappa in (0.3, 1] at tight budgets and heavy tails.
- Cost convex in the floored share everywhere (e.g. heavy, b=0.5: 0.6 / 1.5 / 3.6 /
  6.1 percent of no-funding output at floored shares 0.1 / 0.25 / 0.5 / 0.75).
- Exact endpoint: a full persistent floor IS uniform funding (P-A2, delta = 0).
- Caveat inherited from the campaign: the old "floors cost less than 0.4 percent"
  number came from flooring the WEAKEST funder; the price depends on whose targeting
  is displaced. State as a price, never a prohibition (standing do-not-claim).
- D3 regime pairing caveat: the heavy-tailed cells run with sharp review (tau = 0.3),
  the default cells with default review (tau = 1); cross-field comparisons from D3
  carry that difference. The single-round computation below varies one thing at a
  time.

Value of targeting, the relative object (D3, T=2): informed targeting's gain over
uniform funding, as a fraction of uniform funding's own gain over no funding,
decreases in the budget: 1.17 to 0.72 (b = 0.1 to 1) in the heavy+sharp regime; 0.74
to 0.25 in the base regime. This is the object whose decrease matches \S2's promise
and Corollary 2's asymptote; the absolute gain (percent of no-funding output)
increases over this budget range (8 to 37 percent heavy+sharp; 5 to 11 base), because
b <= 1 is still on the rising side of the curve. The draft must pick objects
explicitly; I propose the relative object as the headline and the absolute as
footnoted context.

Lotteries, analytic: for strictly concave per-researcher output (Lemma 1), any
lottery is weakly outperformed in expected output by the equal split of the same
budget over the same pool, researcher by researcher (Jensen: the lottery gives each
pool member the split's grant in expectation, as a two-point gamble). In particular a
full lottery at any winner count is outperformed by uniform funding, so its price is
at least uniform funding's price: the full value of targeting. Three-line proof;
appendix-ready.

Lotteries, computed (NEW single-round pricing computation, price_lotteries.py: n=50,
2000 populations per cell, expected output exact given the allocation; screen
s = K + noise(tau); funded share q = 0.2; schemes: screened equal split = top-q by
screen at B/qn each; screened lottery = pool of top-2q, half funded at random at the
same grant, the design New Zealand and the German line run; uniform; full lottery.
Metric: share of complete-information targeting's value captured):

| field, budget | screen tau | scr. split | scr. lottery | uniform | full lottery |
|---|---|---|---|---|---|
| heavy (1.3), b=0.1 | 0.3 / 1 / 3 / 10 | 0.78 / 0.75 / 0.65 / 0.51 | 0.60 / 0.56 / 0.51 / 0.45 | 0.44 | 0.38 |
| default (2), b=0.1 | 0.3 / 1 / 3 / 10 | 0.79 / 0.71 / 0.58 / 0.49 | 0.62 / 0.57 / 0.51 / 0.46 | 0.52 | 0.43 |
| even (3.5), b=0.5 | 0.3 / 1 / 3 / 10 | 0.67 / 0.56 / 0.47 / 0.43 | 0.54 / 0.49 / 0.45 / 0.42 | 0.83 | 0.41 |

What the computation shows, in claim form:
- L1: a full lottery is the most expensive scheme everywhere (Jensen made visible).
- L2: where capability is heavy-tailed and the screen at least moderately
  informative, screened schemes capture most of targeting's value, and the screen
  degrades gracefully with noise (the \S6 whale-separation mechanism, now doing
  policy work).
- L3: randomizing within the screened pool has a real price: 0.05 to 0.19 of
  targeting's value below the same screen's equal split, across the cells above.
  "Randomize among the fundable" is not free in output terms; its case must rest on
  the savings side.
- L4: where capability is spread evenly and the budget not tight, the SCREEN is the
  mistake: uniform funding captures 0.83 of targeting's value and beats every
  concentrated scheme. There, a lottery among all applicants is nearly free because
  targeting itself is nearly worthless. This is the precise content of "cheap exactly
  where review is worth little," and it is stronger than I expected: spreading is
  close to optimal, not merely cheap.

## 3. Proposed structure (five paragraphs, two figures)

P1. Re-introduce the two schemes as policies inside the model (floors; lotteries full
and screened, as funders run them), and the scope sentence: the model prices what a
scheme forgoes in expected output; what a scheme saves (review costs, proposal
effort, bias) is priced by others and enters the ledger at the end.
P2. The unifying price: both schemes forgo targeting; restate the value of targeting
in place and give its field dependence (relative object headline). Thesis: the
conditions that make review worth little are the conditions that make forgoing
targeting cheap.
P3. Seed floors: the law, the convexity, the exact uniform endpoint, price-not-
prohibition. FIGURE A: floor cost (percent of no-funding output) vs floored share,
curves for the two regimes x two budgets (D3 data).
P4. Lotteries: the Jensen proposition (appendix pointer); then the computed pricing:
L1-L4. FIGURE B: share of targeting's value captured vs screen noise; curves =
screened split, screened lottery, uniform, full lottery; two panels or two line
families for heavy vs even fields (design below).
P5. The ledger: forgone value (ours) against savings (the contest model's); when
lotteries are defensible, answered field by field; bridge "Next we vary the
production technology..." (\S9).

## 4. Decisions I need before drafting

1. THE PRICING COMPUTATION AS GROUND. It evaluates policy schemes on the model's
   primitives in a single round, exactly (the mechanism-test precedent, C6/C7), not
   through the full T=2 Bayesian simulation. Alternative: implement the schemes as
   strategies in model.R and run them through the sim (heavier; a model-code change
   that plausibly belongs to the three of you). I recommend the computation, with its
   single-round scope stated in the generating context.
2. SCREENED-LOTTERY DESIGN. I used: pool = top 40 percent by screen, fund half at
   random, grant size equal to the split's (the NZ shape). Knobs: pool size, funded
   share q (0.2 fixed here), grant size. Also whether to add the same-pool equal
   split as a displayed line to isolate the randomization price (I lean yes; it is
   the cleanest reading of L3).
3. JENSEN PROPOSITION: state in body with proof in Appendix A (as a second
   proposition), or footnote it? I lean body + appendix: it is the section's one
   analytic claim and the figure-or-proof rule wants it in the appendix.
4. THE ABSTRACT SENTENCE. "Cheap exactly where review is worth little" is licensed in
   the field dimension and (via the relative object) asymptotically in the budget;
   L4 sharpens it (spreading is nearly optimal there, and the screen itself is the
   mistake). Flag: revisit the abstract's wording when \S8 locks.
5. FUNDED-SHARE HONESTY. Uniform beats the q=0.2 screened split at default tails with
   b=0.5 and tau >= 1 (0.69 vs 0.70 at tau=1 is a tie; uniform wins at tau=3). The
   funded share is a design knob the section holds fixed; a fuller treatment would
   optimize q per scheme. I propose stating the q=0.2 scope plainly and not
   optimizing; say the word if you want the q sweep.

## 5. Do-not-claim (standing + new)

- Floors as prohibition; state as price. The 0.03-0.4 percent legacy number only with
  its flooring-the-weakest provenance.
- No cross-field D3 comparisons without the tau pairing caveat.
- Nothing about review costs, bias, or proposal effort as results of ours; the
  savings side is cited, not modeled.
- L4's "spreading nearly optimal" scoped to even spread + the tested budgets; the
  computation is single-round.
