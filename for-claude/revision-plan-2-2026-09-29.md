# Revision plan 2 (after the second review), 2026-09-29

Applies to both documents (tex/main.tex, tex/pnas/pnas-main.tex). Status: awaiting Aydin's approval.

## Verdict on the review

Right and worth acting on: points 1 (three notions of optimal), 2 (demote timing), 3 (signal model, tau_R, the "independent capability signal" reading), 4 (particle audit, heavy tails, rho's definition), 5 (qualify the two-features thesis; latent vs observable), 6 (lottery reframed around lost ranking), 7 (output-maximizing; resources do not persist), 8 (literature; regime map with real programs), and every item under organization and textual issues. Declined: renaming K (one concept, one label; the title names the capability-resource gap), truncated Pareto (medians and a copula definition of rho answer the concern), an exact dynamic program for the certainty-equivalent planner and a terminal value for capability (moot once timing moves to the SI; both stated there as limitations).

## Checks already run (scratchpad/renorm)

1. dynopt.R: two-round complete-information dynamic optimum (free choice of round-1 recipients and the round-1 spending share; gap rule on the remaining budget in round 2) against the paper's benchmark (gap rule each round, budget split evenly). Over alpha_K in {1.3, 3.5}, epsilon in {0.1, 0.85}, b in {0.2, 1}: the benchmark loses at most 1.0% of funding's gain (0.0% at epsilon = 0.1); round-1 grants correlate with the dynamic optimum's at >= 0.99; the dynamic optimum spends 35-49% of the budget in round 1 (later than half, because grown capability makes round-2 dollars more productive). So the one-round gap rule is the right benchmark and the SI can say so with numbers.
2. ess.R: importance sampling at M = 200 vs 2000, one round of output plus both signals. Median effective sample size at tau_K = 0.05 is 19 (10th percentile 1.4; some posteriors collapse to one particle); at tau_K = 0.3, 79; at tau_K >= 1, > 120. Posterior means of top-decile capability at M = 200 are biased low by 8.6% (alpha_K = 1.3, tau_K = 0.3) against 1.5% at M = 2000; at tau_K >= 1 the bias is prior shrinkage and is the same at both M. The referee is right that the low-noise end of Fig. 2C rests on thin posteriors; the direction of the bias understates review's value in the heavy-tailed field.
3. Cross-references: after timing moved behind spreading, the SI settings table and SI Methods still hard-code Fig. 3 = timing and Fig. 4 = spreading. Confirmed swapped. rho is a Gaussian-copula parameter in the code (draw_initial_population), not Pearson.

## Package A: optimality language and the dynamic benchmark (text + one SI table)

- Proposition 1 is the one-round, complete-information, output-maximizing allocation. Abstract, Significance, Results heading, Discussion: "the optimal allocation" -> "the output-maximizing allocation of a round's budget funds the capability-resource gap"; the sequential problem uses it as the benchmark.
- New SI paragraph + Table (SI Text 1): dynopt results at T = 2 on the full grid (alpha_K x epsilon x b, 20 populations), plus a T = 5 check restricted to spending shares (gap rule recipients each round).
- Everywhere a prescription is stated as advice (Significance, Discussion, Conclusion): "output-maximizing" once per paragraph; the objective is aggregate expected research output and nothing else, said once in Model and once in Discussion.
- Model: capability is a state variable (expertise, accumulated research capital), changed by doing research; resources recur and do not persist, and the asymmetry is named as an assumption with its consequence (static results unaffected; dynamic results depend on it).

## Package B: demote timing (structure)

- PNAS: remove "Strategic timing of funding" and Fig. 4 from Results; one sentence in Discussion ("Planning the budget over time adds little beside choosing whom to fund; SI Text 6"). New SI Text 6 "Strategic timing of funding": the current text, Fig. S (both panels), with three stated limits: the planner is certainty-equivalent and does not value the information its grants generate; its gain over the myopic funder combines timing and recipient choice, which we do not separate; there is no terminal value for capability, which favors early spending, so the late-spending finding is conservative.
- Long version: Section "Strategic timing" moves to an appendix (Appendix F) with the same limits; Discussion sentence as above.
- Terminology: "optimal timing" -> "the forward-looking funder's timing"; "timing gain" / "value of timing" -> "gain from planning ahead" (Aydin to approve: this changes a locked term because the referee is right that the comparison is not timing alone).
- SI cross-references: all figure numbers by \ref; settings table rows re-labelled.

## Package C: review-signal robustness and resource information (runs + SI Text 5 additions)

Runs (T = 2, b = 1, fixed-mean fields, 50 populations, M = 200 unless Package D changes M):
- C1 tau_R in {0.1, 0.3, 1, 3, 10} x alpha_K in {1.3, 2, 3.5}: value of targeting under records-only and with review; how much a resource signal is worth beside review.
- C2 multiplicative review noise, log S = log K + eps, sd set so that the rank correlation with K matches the additive default (and two other levels): value of review by field.
- C3 persistent review error: S_i drawn once and reused every round (a reviewer's fixed misjudgment). T = 2 by field, and the T = 20 by-round shortfall of Fig. 2B with persistent against fresh review.
Text: Results paragraph on review closes with the referee's reading, adopted: "what review supplies is a signal of capability separate from resources; peer review is one source of such a signal, and its value is greatest where the gaps are most unequal." One sentence each on C1-C3 in Results; tables in SI Text 5.

## Package D: numerical audit (runs + SI Methods)

- D1 particles: tau_K in {0.3, 1, 3} x alpha_K in {1.3, 2, 3.5} at M in {200, 1000, 5000} (50 populations; 25 at M = 5000), T = 2: value of review and of targeting, with standard errors; median and 10th-percentile ESS reported. Fig. 2B (T = 20) rerun at M = 1000. Decision rule: if any headline number moves by more than its standard error, rerun Fig. 2C at the larger M and report that; otherwise the audit table justifies M = 200.
- D2 heavy tails: alpha_K = 1.3 headline cells at 200 populations, medians beside means; SI Methods notes the infinite variance at alpha < 2 and the realized means.
- D3 Methods: rho is the correlation parameter of a Gaussian copula (rank correlation stated for the values used); "holding mean capability fixed, the tail parameter varies dispersion and tail inequality" replaces "varies inequality alone".

## Package E: thesis qualification, regime map, operationalization (runs + text + one new main figure)

- E1 regime-map grid: b in {0.05, 0.1, 0.2, 0.5, 1, 2, 3} x alpha_K in {1.2, 1.3, 1.5, 2, 2.5, 3.5, 5}, fixed mean, T = 2, 50 populations: value of targeting (complete information over uniform, share of funding's gain) and value of review at tau_K = 1. New PNAS Fig. 4 (replacing timing): x = budget relative to the field's resources (log), y = capability inequality (Gini of the capability distribution, computed from alpha_K), contours of the two values; shaded boxes placing programs by order of magnitude (a national biomedical funder with a large share of its field's support; a national agency in a field it funds thinly; a European excellence council; a private foundation), with b estimated from budgets against field research spending (web-sourced at execution, cited, with the uncertainty stated) and y given as a range from productivity-tail estimates. Long version: same figure in the Discussion.
- E2 Discussion: distinguish the sufficient statistic (inequality of the optimal gaps, latent, budget-dependent) from observables a funder can estimate (its budget against the field's resources; productivity tails corrected for resources; review-score dispersion; seed or pilot grants, which reveal gaps directly). Shorten the recapitulation of Results to make room.
- E3 Significance and Discussion: the two features govern the value of targeting and of information; concentration depends in addition on complementarity (Fig. S6). Remove "concentrate or spread" from the two-features sentence in Significance.

## Package F: lottery reframing and literature (text)

- F1 Results (spreading): lead with the partial-lottery finding, that most of the loss comes from discarding the ranking within the fundable pool, and state the equal-division comparison second, with the minimum-viable-grant condition in the same sentence. SI Text 2 unchanged.
- F2 Intro/Discussion citations: Carnehl, Ottaviani & Preusser (Designing Scientific Grants, NBER w32668 / arXiv 2410.12356); the extended lottery taxonomy (Research Evaluation 2024, doi 10.1093/reseval/rvae025); the New Zealand randomized trial of winning funding (Science and Public Policy 51(6):1042, 2024); Bol, de Vaan & van de Rijt 2018 already cited, now contrasted with capability compounding. Bib entries verified against the DOIs at execution.
- F3 Aydin's items unchanged: grant number, doubled period, bibliography CHECK notes, repository citation.

## Order and cost

Runs first, in the background, in this order: D1 (decides M for everything after), E1, C1-C3, dynopt full grid. Roughly 6-8 hours on two cores. Text packages A, B, F, E2-E3 proceed while the runs go; C, D, E1 text and figures follow their results. Ledger: revision-ledger-2-2026-09-29.md, every changed number and claim.
