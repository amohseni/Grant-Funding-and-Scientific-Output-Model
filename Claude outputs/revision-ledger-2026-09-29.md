# Revision ledger (2026-09-29): what changed in response to the review

Every number below was computed in-session from the staged model (model.R, smooth
allocator) or from the staged sweep summaries; scripts sens.R, sens2.R, byround.R in
the session scratchpad, results in sens_results.csv, sens2_results.csv, byround_results.csv
(copied to for-claude/analysis/ with this ledger).

## Package A: corrections (both versions)
- Lemma 1: second derivative -2A K^2/S^3 (was -4A after the factor-2 removal).
- Corollary 2: vanishing only; monotone decrease stated as a numerical finding (ratio to
  uniform's gain fell monotonically on b in [0.02, 3] in 300 of 300 random populations).
- Lower bound on capability now stated for expected output; realized output bounds less.
- Methods and SI state the likelihoods (Gaussian signals, Poisson output), the 200-particle
  importance sampler, and the planner (certainty-equivalent receding horizon; does not
  value the information its own grants generate).
- Lottery comparison: divisibility and diminishing returns stated at use; minimum-viable-
  grant remark added to Appendix B / SI Text 2; indivisibility added to limitations.
- Scope: investigator-level funding (model and limitations).
- Sorzano and Pueche-Granados engaged in the Discussion (contrast: observed impact vs
  latent capability and resources).
- Modal: "Within the model, no funding mechanism is best in every field"; "the model's
  best use of the budget..."; "under complementarity, observation comes first and money
  second"; "the timing of spending is worth less than the signal...".
- Timing moved after spreading in both versions; roadmap and four-questions order match.
- Parameter justification: productivity tails are an imperfect guide to capability.
- Allocator sentence in the settings table: greedy allocator is the only one implemented
  for the general family (a rerun with the water-filling allocator is not possible; it is
  harmonic-only).

## Package B: the general gap rule
Proposition 1 now holds for F(K, R) = A K h(R/K), h increasing and strictly concave
(constant returns, diminishing returns to resources); c = (h')^{-1}(nu/A); harmonic case
c = sqrt(A/nu) - 1. Corollary 2 requires h bounded (true for gamma < 0 in the power-mean
family, false for Cobb-Douglas, where the value of targeting grows without bound in B for
unequal capabilities: proof sketch by Cauchy-Schwarz in the remark). Leontief excluded
(not strictly concave); treated numerically in Appendix C.

## Package C: sensitivity runs (new Appendix E / SI Text 5, Figs 15 / S7, Tables)
- Mean-normalized capability sweep (E[K] = 2 for every tail) REPLACES the Fig 4 / 2C data.
  Retention of review's maximum value: at default noise 82 / 73 / 46 percent (heavy /
  default / even); at three times the default noise 51 / 32 / 8 (was 82 / 44 / 10 at
  fixed scale, which confounded inequality with a threefold difference in mean K).
  Records-only shortfall as share of complete-information gain: 45 / 27 / 11 (was 40 / 27 / 12).
  Review's maximum value as share of the gain from funding: 43 / 24 / 8 (was 39 / - / 9).
  Grant correlations with the gap rule: default 0.20 to 0.94 (was 0.22 to 0.94); heavy
  0.27 to 0.96 (was 0.34 to 0.96). Shortfall range: default 26 to 3.7 (was 25 to 3.6);
  heavy 40 to 3.2 (was 31 to 2.8). Error bars (100 populations) now drawn.
- Matched informativeness (rank correlation 0.8 / 0.55 / 0.3): heavy 0.91 / 0.78 / 0.48;
  default 0.83 / 0.63 / 0.30; even 0.64 / 0.45 / 0.17. Ordering unchanged, gaps smaller.
- Capability-resource correlation rho in {-0.5, 0, 0.5, 0.8}: targeting value heavy 0.96
  to 0.69, default 0.56 to 0.25, even 0.25 to 0.06; review share heavy 0.89 to 0.84,
  default 0.64 to 0.54, even 0.22 to 0.00. Ordering unchanged at every rho; magnitudes
  fall as resources align with capability. The referee's two-agent counterexample
  reproduces (0.054 vs 0.633).
- Resource tail alpha_R in {1.3, 2, 3.5}: changes targeting by at most 0.1, review by less.
- Pool size n in {25, 50, 100, 200}: targeting rises with n in the heavy field (0.81 to
  1.15); review's share rises with n in the even field (0.14 to 0.30); ordering unchanged.
- Review contaminated by resources, S = K + beta R + eta, beta in {0, .25, .5, 1}, funder
  knows beta: heavy field retains 0.60 to 0.65 of the shortfall at beta = 1 (0.77 to 0.91
  clean); default 0.38 to 0.45; even 0.04 to 0.11.
- Gap inequality as the organizing quantity: over all 159 settings, log(value of
  targeting) on the Gini of the optimal grants has R^2 = 0.96 (0.99 with log(E[K]/E[R]));
  logit(review share) on Gini + logit(rank correlation) has R^2 = 0.86, and adding the
  capability tail as a factor adds < 0.01.
- By-round shortfall (T = 20, default field): records-only 47 to 20 percent, with review
  9 to 4 percent (heavy, mean 2: 61 to 16 vs 7 to 3). Fig 5 / 2B now plot this; the grant
  correlations (0.08 to 0.41 vs 0.81) stay in the text.
- NOT run: persistent reviewer error; project quality; calibration; the robustness rerun
  on the water-filling allocator (harmonic-only).

## Package D: the thesis
Second feature is now "how unequally the capability-resource gaps are distributed",
with capability inequality as its main driver and the capability-resource association
as the second, working through the same channel. Applied in the Significance (117
words), abstract, introduction, Discussion, conclusion of the PNAS version and in the
abstract, introduction, Discussion, conclusion of the long version.

## Still yours
Grant number; duplicate period in the funding line; bibliography CHECK notes; the
Nicholls range; Li and Agha funded-only confirmation; the repository URL. Decide whether
timing stays in the main text (kept, after spreading) or moves to the SI.
