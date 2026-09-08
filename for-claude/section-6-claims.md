# Section 6 claims ledger (updated 2026-09-08, full verification pass)

Every claim in \S6 v4, formulated precisely, with its modal force, its ground, and its
verification status. "Mechanism test" = direct computation on the model's primitives
(verify_s6_mechanisms.py). "Re-derived" = recomputed 2026-09-08 in container R from
model.R (hash 21b0d9a), after validating bit-identity against the canonical
D_gap_convergence cells (max diff 1e-14); new runs: fine tau grid
(gap_convergence_fine.csv, 17 tau x 3 alpha x 50 seeds), multi-round trajectories
(verify_s6_rounds.R, T = 8 and 20), paired sharp-end test. Nothing locks until its row
is green.

| # | Claim (as formulated in the text) | Modal force | Ground | Status |
|---|---|---|---|---|
| C1 | Review is a noisy signal of capability (tau_K clause now removed from the opening sentence; tau_K introduced in \S3) | definitional (\S3 spec) | model definition | OK |
| C2 | Value of review = expected value of its reduction of the funder's uncertainty regarding the capability-resource gap; quantified as expected output with review minus without | definitional | definition + scoring convention (\S3) | OK |
| C3 | A more informative signal brings allocation closer to the optimal allocation, output closer to the complete-information benchmark | directional over the sweep, with a footnoted sharp-end reversal | fine grid: corr_S5 and value increase as tau decreases, except tau = 0.05 vs 0.3 under heavy tails: value ~2% lower, paired z = 4.5 (defaults z = 0.1) | VERIFIED (re-derived); reversal footnoted |
| C4 | With a sharp signal the funder recovers MOST of the records-only shortfall (was "nearly all"; recovery 86% defaults, 93% heavy); corr 0.22 -> 0.94 defaults, 0.34 -> 0.96 heavy; shortfall 10.4 -> 1.6% of no-funding output defaults, 17.2 -> 1.7% heavy | measured at stated parameters | fine grid + canonical D-4 (old footnote had 0.95/0.97 and misattributed 17.2/1.7 to defaults; both corrected) | VERIFIED (re-derived, corrected) |
| C5a | Records-only allocation improves slowly across rounds and stays far behind: corr 0.08 -> 0.41 over 20 rounds defaults (0.07 -> 0.43 heavy) vs review-informed 0.81 (0.92) in round one | measured, T = 20, tau_K = 1, 50 seeds; "does not catch up" scoped to simulated horizons | verify_s6_rounds.R (the old claim "no closer across rounds, 0.13-0.18" was FALSE: that range was across capability distributions at round 1) | VERIFIED (re-derived, claim corrected) |
| C5b | Structural cause: capabilities compound, resources do not accumulate, so output approaches the resource-limited level and carries less and less information about capability | analytic (limit + derivative) | lambda -> 2AR; dlambda/dK -> 0; per-round Fisher info 7.9e-2 -> 6.1e-7 over 200 rounds; cumulative info converges (~0.44 over 20k rounds) | VERIFIED (mechanism test) |
| C5c | (Not claimed in text) mechanism of the residual slow improvement in C5a not isolated (candidates: finite-information accumulation; growing K dispersion) | n/a | text attributes nothing | OPEN, flagged to Aydin |
| C6 | Where capability is heavy-tailed, most of the output a funder can add comes from resourcing the few most capable | "most," heavy-tailed case | top-10%-by-K share of optimal-allocation gain at b=0.1: 0.92 / 0.86 / 0.72 (alpha 1.3 / 2 / 3.5); 400 draws each | VERIFIED (mechanism test) |
| C7 | The most capable are easiest to identify; even a noisy signal separates them | comparative | s = K + N(0,tau): top-10% recovery at tau=3: 0.68 heavy vs 0.21 light; middle-pair ordering 0.50-0.54 | VERIFIED (mechanism test) |
| C8 | Even a signal of moderate informativeness captures most of review's value under heavy tails; sharpening beyond a modest level adds little | measured | fine grid: value retained at 3x default noise: 82% / 44% / 10% (alpha 1.3 / 2 / 3.5); half-value noise ~ tau 10 / 2.5 / 0.75 (old "elbow at tau ~ 2.5" was the defaults' half-value point, not a kink; rescoped) | VERIFIED (re-derived, rescoped) |
| C9 | Corollary: fine distinctions among middling applications add comparatively little; heavy investment there MAY be misallocated | "in our model" + "may" | follows from C6 + C7 + C8; uniform-tau caveat stands (rank-local discernment not separately simulated) | MECHANISM VERIFIED; caveat retained |
| C10 | Where capability is spread more evenly, review of any informativeness adds little | comparative + measured at alpha 3.5 | records-only shortfall 3.5% of no-funding output at alpha 3.5 (11.4 defaults, 24.6 heavy); max recovery 72% of that; review's max value 2.5% of no-funding output vs 23% heavy | VERIFIED (re-derived; was PARTIAL) |
| C11 | Where the budget is ample, information about whom to fund has little left to add | analytic limit + measured | Corollary 2 (proven, appendix) + targeting-value curve | OK (analytic) |
| C12 | AUC 0.54 corresponds to tau_K >~ 20 at the default parameters; the accuracy-to-noise mapping is field-dependent (same tau, higher AUC under heavy tails) | calibration mapping, scoped to defaults | auc_grid.csv: alpha 2, tau 20: AUC 0.56 (top decile) / 0.54 (top quintile); alpha 1.3, tau 20: AUC 0.66 | VERIFIED |
| C12b | At tau = 20, review recovers 31% / 9% / 0.1% of the records-only shortfall (alpha 1.3 / 2 / 3.5) | measured | fine grid | VERIFIED (re-derived) |
| C12c | Closing sentence recast to the licensed pair: same informativeness -> several times the value; same accuracy -> different informativeness ("same accuracy worth a great deal in one field" was NOT licensed: AUC 0.54 under heavy tails maps to tau >> 20, off the sweep) | comparative | C12 + C12b | VERIFIED as recast |
| C13 | Overtrust can turn review's value negative at moderate tails; attenuates at heavy tails | "can"; cell-noise caveat in footnote | D-2: alpha 2, true tau 3, belief <= 1: -0.53 / -0.43 / -0.14, each within 2 SE of zero, all negative; alpha 1.3: 14.0 vs calibrated 17.2 (drop z ~ 2.7); 200 seeds | VERIFIED (caveat kept) |
| C14 | Overtrust costs more than undertrust forfeits | asymmetry, measured | tenfold overtrust (true 3, belief 0.3, alpha 2): -0.43 vs calibrated 3.38 (loses more than full value); tenfold undertrust (true 0.3, belief 3): 4.21 of 7.63 retained | VERIFIED |

Notes.
- C5 history: Aydin's sanity check (2026-09-06) prompted the reformulation recorded in
  the previous ledger; this pass then falsified the interim "no closer across rounds"
  claim by direct measurement (see C5a). The lesson both times: measure the trajectory,
  do not infer it from the mechanism.
- Sharp-end reversal (C3): real but small; kept as footnote pending Aydin's verdict
  (keep / drop / investigate).
- All re-derived numbers come from runs whose scripts and outputs are in the repo
  (for-claude/run_D4_fine.R, for-claude/gap_convergence_fine.csv,
  for-claude/verify_s6_rounds.R); canonical coarse-grid cells match
  sweep_results/D_gap_convergence/gap_convergence.csv to 1e-14.
- Figure 4 carries C3, C4, C8, C10, and the main result; see draft note 8 for the
  ranked alternatives and why no AUC marker appears on it.
