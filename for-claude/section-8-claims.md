# Section 8 claims ledger (v1, 2026-09-08)

Grounds: Package A (D3_seed_signal, D4_seed_persistent; 200 seeds, T=2; verify_s8.R,
all re-derived in-session) and the single-round scheme-pricing computation
(price_lotteries.py: n=50, funded share 0.2, screen s = K + N(0, tau), 2000
populations per cell, max MC SE 0.002; expected output exact given the allocation).
Metric for the computation: a scheme's gain over no funding as a share of the optimal
(complete-information gap-rule) allocation's gain.

| # | Claim (as formulated in the text) | Modal force | Ground | Status |
|---|---|---|---|---|
| E1 | Informed targeting adds more than uniform funding's entire gain over no funding (heavy tail + sharp review + tight budget) down to about a quarter of it (default field, b=1) | measured, relative object, ranges over b | D3: (S5-S2)/(S2-S1) = 1.17 -> 0.72 heavy+sharp; 0.74 -> 0.25 base; regime pairing caveat footnoted | VERIFIED |
| E2 | Thesis: the conditions that make review worth little make forgoing targeting cheap | comparative, field dimension; budget dimension via Corollary 2 asymptote | E1 + S6 results + Corollary 2 | OK as stated |
| E3 | Floor cost increases faster than the floored share (convex) | measured, all 8 (regime x b) curves | D3 cost columns; P-A3 | VERIFIED |
| E4 | Floor cost scale set by targeting's value: ~ informed gain over uniform x floored share x factor 0.10-0.86 (above 0.3 at tight budgets + heavy tails) | approximate law, range stated | kappa computation in verify_s8.R | VERIFIED |
| E5 | Focal floor cost: ~6% of no-funding output (heavy+sharp, b=0.5, three-quarters floored); ~1% default field | measured | D3 (6.11 / 1.22) | VERIFIED, figure-backed (fig10) |
| E6 | Full floor = uniform funding exactly | implementation identity | P-A2 (delta = 0) | VERIFIED (exact) |
| E7 | Any lottery over a pool is weakly outperformed by the equal division of the same money over that pool; full lottery <= uniform funding; equality only at winners = pool | analytic (Jensen + Lemma 1(ii)) | Proposition 2, Appendix B (tex/appendix-lottery.tex) | OK (proof written) |
| E8 | Heavy+tight: equal division among the top fifth captures ~3/4 of the optimal gain at sharp screens, ~1/2 at very noisy; degrades gracefully | measured | computation: 0.776 -> 0.512 over tau 0.3 -> 10 | VERIFIED, figure-backed (fig11 left) |
| E9 | Screened lottery's shortfall vs the division among the top fifth is mostly the discarded ranking, not the gamble | measured decomposition | at tau=1 heavy: total 0.18 = 0.14 ranking + 0.04 gamble (top-fifth division 0.748, doubled-pool division 0.605, lottery 0.564) | VERIFIED |
| E10 | Even spread + ample budget: ordering inverts; uniform captures > 4/5 of the optimal gain (0.83) and beats every concentrated scheme; spreading is what pays | measured | computation, alpha=3.5, b=0.5 rows | VERIFIED, figure-backed (fig11 right) |
| E11 | Defensibility (comparative): replacing the ranking with chance costs 0.01-0.07 of the optimal gain where capability even or screen very noisy; 0.14-0.18 where heavy-tailed with sharp-to-moderate screens | measured, range over cells | computation: lottery-vs-ranked gaps by cell | VERIFIED |
| E12 | "Where capability is evenly spread, the model's answer is not a lottery among a fundable few but thinner grants to many" | prescriptive gloss of E10, single-round scope | E10 | OK given E10; scope stated |

Notes.
- The ABSOLUTE defensibility claim ("lotteries are cheap at even spread") is FALSE in
  the computation and is not made: even at alpha=3.5 the screened lottery captures
  0.42-0.54 vs uniform's 0.83. E11 is the licensed comparative claim. Abstract's
  "cheap exactly where review is worth little" flagged for the lock-time abstract
  pass (draft note 3).
- Terminology decision open (draft note 2): "ranked/wide division" vs "ranked/wide
  split" between body and figure.
- All computation cells single-round; the funded share held at one fifth by decision
  (stated, not optimized). D3 rows are T=2 with the heavy=sharp-review pairing.

Addendum (2026-09-08, plainness rewrite): claims unchanged; wording aligned to the
rewritten body (E3 now stated concretely: half the budget as seed grants loses more
than twice the output a quarter loses; verified on all eight D3 curves, ratios
2.1-2.8). Terms: "fraction of the budget given out as seed grants" (not "floored
share"), "selection by review" / "partial lottery" / "lottery over all researchers"
(not "screened lottery"/"full lottery"), "output lost" (not "price"/"cost").
