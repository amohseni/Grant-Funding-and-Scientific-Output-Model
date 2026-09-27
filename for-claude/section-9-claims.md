# Section 9 claims ledger (v1, 2026-09-08)

Grounds: Package B sweeps (sigma_tierA_gc0 / _gcm3 / _leontief: T=2, greedy allocator
n_steps=400, seeds 100/100/50; sigma_tierB_concentration + leontief endpoint +
400-seed refinement: T=1, S5, 200 seeds), all re-derived in-session from the staged
canonical files (verify_s9.R + verify_s9_OUTPUT.txt); the Cobb-Douglas planning cells
re-derived earlier this session (verify_s7.R). Within-family comparisons only: the
greedy-allocator magnitudes are not quoted against smooth-run magnitudes.

| # | Claim (as formulated in the text) | Modal force | Ground | Status |
|---|---|---|---|---|
| G1 | At every technology tested (CD, gamma=-3, Leontief), the review signal's value increases with capability inequality and with informativeness | monotone across the tested grids | signal_value summaries: all three gamma, monotone in k_shape at every tau and in tau at every k_shape (where value distinguishable from zero) | VERIFIED, figure-backed (fig12) |
| G2 | The signal's value increases with complementarity: 6.8 / 29.3 / 33.7 % of no-funding output at CD / gamma=-3 / Leontief (heavy + sharp) | measured ordering; same ordering across cells | signal_value summaries | VERIFIED |
| G3 | Cobb-Douglas: scheduling gain <= 0.005% of no-funding output (T=5, all compounding rates), schedule within 0.012 of even; gamma=-3: gain ~1.5% and center of mass 0.66 at the highest rate | measured, scoped | sigma_tierA_gc0 / _gcm3 horizon_growth (verify_s7.R re-derivation) | VERIFIED |
| G4 | Tight budget + heavy tail: concentration of the informed funder's grants increases toward Leontief (Gini 0.916 -> 0.938, Leontief 0.947) | measured | tierB k=1.3 b=0.1 row + leontief endpoint | VERIFIED, figure-backed (fig13) |
| G5 | Ample budget: concentration decreases toward Leontief (default field: 0.354 -> 0.285, Leontief 0.207; heavy field at b=1: 0.631 -> 0.528, Leontief 0.448) | measured | tierB b=1 rows | VERIFIED, figure-backed (fig13) |
| G6 | Evenly spread + ample: interior maximum near gamma ~ -6 (0.196 vs 0.185 at -12; peak rise z = 3.5 at the 400-seed refinement); declines toward both ends | measured, z stated | tierB k=3.5 b=1 + refine (-9: 0.195, -4: 0.193) | VERIFIED, figure-backed (fig13) |
| G7 | (Available, not claimed) family seed-floor costs: worst re-derived cost 0.32% (CD) / 0.44% (gamma=-3) of no-funding output; exact strategy pair of the seed_myo column unverified at column level | held out of the body (draft note 5) | seed_value summaries | PARTIAL, not cited |
| G8 | (Available, not claimed) correlation robustness across the family | not re-derived | correlation_summary staged, unread | OPEN, not cited |
| G9 | "The technology does not settle the concentration question by itself; the budget and the field interact with it at every point" | summary of G4-G6 | G4-G6 | OK given G4-G6 |
| G10 | Mechanism glosses (bottleneck top-up; well-vs-adequately placement) | interpretive, flagged | draft note 4; directions verified, mechanisms not isolated | FLAGGED for Aydin |

Notes.
- Middle-budget (b=0.5) rows are mixed by tail (heavy decreases, even increases);
  the body's claims are stated at the swept extremes (tight b=0.1, ample b=1) that
  the figure displays; the b=0.5 rows are in verify_s9_OUTPUT.txt.
- Leontief kink caveat (ties by fill order, finer allocation step) footnoted in the
  generating context and figure comment.
- Story-5 framing ("the concentration dispute is a disagreement about
  substitutability") is refuted by G4-G6 and is presented only as the tempting story
  the section corrects (campaign flag honored).

Addendum (2026-09-09, restructure): G1-G2 now live in Appendix C (fig12 renders as
Fig 13); G3 now carried by Table 2 (Appendix C) with the full T=5 rows re-read from
sigma_tierA_gc0/_gcm3 horizon_growth summaries: CD gain %S1 by eps = -0.001, -0.002,
0.001, 0.001, 0.004 (SE <= 0.0015), b_idx 0.499-0.512; gamma=-3 gain 0.01, 0.13, 0.45,
1.07, 1.54 (SE <= 0.09), b_idx 0.511-0.660. The S7 footnote that quoted these beside
smooth-allocator magnitudes is replaced by the appendix reference (cross-allocator
quoting removed; comparison is CD vs gamma=-3 within tier A). G4-G6, G9 now in S8 body
(fig13 renders as Fig 12); G6 stated with 0.196 peak vs 0.185 at gamma=-12 (my chat
remark that the effect was 0.003 compared the peak with its NEIGHBORS, not the
endpoint; the endpoint gap is 0.011, z=3.5). G10a (bottleneck top-up) kept as
"a pattern consistent with" (hedged, not asserted); G10b cut. Provenance gap: the
sigma_tierA_gcm3 folder has no RUN_INFO.txt; gc0 hash d77c854 recorded.

## Discussion (S9) claims, H1-H14 (2026-09-09; draft: section-9-discussion-draft-v1.md)

No new results; every row cites a body figure, an appendix, or a cited paper.

| # | Claim | Modal force | Ground | Status |
|---|---|---|---|---|
| H1 | Value of targeting is set by budget scale and capability inequality | structural, "in every field we examine" | Cor 2 (App A), Fig 2; S6 Fig 5 | OK |
| H2 | Records obtain little and further rounds add little | measured | S5 Fig 4 (0.08 -> 0.41 by round 20) | VERIFIED |
| H3 | Review's share rises with the two features and informativeness | measured | S6 Fig 5 | VERIFIED |
| H4 | Spending later adds less than the review signal at every horizon and rate swept | measured, range | S7 Fig 9 (PG/signal 0.000-0.66 over 32 cells) | VERIFIED |
| H5 | Deep grants obtain far more than thin in a poor field | measured | S7 Fig 8 (0.29 -> 10.0 %S2 over b 0.5 -> 6) | VERIFIED |
| H6 | Seed grants and lotteries lose most where the stake is largest | measured | S8 Figs 10-11, E4 factor law | VERIFIED |
| H7 | Complementarity raises the stake and downstream quantities | measured within tier A | App C Fig 13 (6.8/29.3/33.7) | VERIFIED |
| H8 | Output distribution understates capability inequality where grants are thin | derived from S5 (output ~ 2AR on thin grants) | S5 Fig 3 | OK |
| H9 | Small-foundation case: coarse review recovers most value; seed grants lose most; partial lottery loses ~1/6 of selection's gain | measured, cell-matched | S6 (82% at 3x noise, alpha 1.3); S8 E5/E11 (0.14-0.18) | VERIFIED |
| H10 | Large-agency case: review adds little; seed grants/lotteries lose little; thin grants to many beat a lottery among few | measured | S6 (alpha 3.5 max 2.5 %S1); S8 E10/E12 (uniform 0.83) | VERIFIED |
| H11 | Efficacy studies measure prediction among funded grants; Fang: top-fifth scores, AUC 0.54 | citation | Fang et al. 2016 text (fetched 2026-09-09); Li & Agha 2015 funded-only by recall, CONFIRM | PARTIAL (Li & Agha) |
| H12 | Therefore S6's accuracy-to-noise mapping is conservative | inference (truncation lowers measured discrimination for a valid signal) | H11 + S6 footnote | OK as "likely"; needs S6 clause |
| H13 | Overtrust loses more than undertrust | measured | S6 Fig 6 | VERIFIED |
| H14 | Three predictions (proxies predict better on deep grants; validity higher in heavy-tailed fields over the whole pool; lottery loss larger in heavy-tailed fields) | comparative predictions derived from S5, S6, S8 | Fig 3; Fig 5 + mapping footnote; Fig 11 | OK (derived) |

Scope paragraph: each held-fixed item checked against S3 (single funder; no strategic
response; scalar Poisson rate; unbiased Gaussian review; fixed n; grants consumed;
harmonic form). Timing scoped to complementarity: App C Table 2.
