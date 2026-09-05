# Section 4 (The optimal allocation) notes ledger

Started 2026-08-14. Content banked for \S4 as it is displaced from other sections.

## Method precedents (moved from \S2 per Aydin's edict: method labels live with the derivation)

Place beside the derivation of g* = cK - R:
- The optimality condition (fund to a common marginal return) is a water-filling solution of
  the kind familiar from information theory (Cover and Thomas, 2006). Cite at the derivation,
  one sentence; positioning per notes-gap-rule-novelty.md (precedent legitimates the
  derivation; no mathematical-novelty claim).
- Funding a shortfall from a target level parallels the poverty-gap transfer in development
  economics [cite: canonical poverty-gap source; check Foster-Greer-Thorbecke]. The
  capability-independent special case (fund the resource gap alone) IS the poverty fill; our
  rule scales the fill level by capability.

## Also banked for \S4 (from earlier intro compressions)

- The funding frontier (K = R/c): term debuts here, not in the intro.
- The formula g* = cK - R debuts here (intro is formula-free by decision 2026-08-14).
- The two intuitive options refuted by dissociation (grant-and-redirect): track record
  rewards output, which ample resources can produce without a gap; funding the
  under-resourced rewards scarcity, which low capability can accompany without a gap.
- Figure: optimal grant across the population (report Fig 1); coverage vs budget (report
  Fig 2, the budget conditional).

## Banked FROM \S4 (2026-08-14, numbers-out decision): for the quantitative sections

Aydin's call: \S4 is qualitative only; every claim there must be derivable from the gap
rule (now all are: Examples 1-2 and Corollary 2 in gap-rule-proof-v1.md). The displaced
simulation numbers, for use where parameter context exists (likely \S8 for the targeting
value curve; strategy comparisons wherever track record is measured):
- Track-record-proportional vs uniform, output gain over no funding at default
  parameters: 23.5 vs 26.1 percent; at heavier capability tails: 30.1 vs 29.8 percent
  (report Table 2; re-derive via verify_all_claims.R before use).
- Value of targeting (optimal over uniform, relative to uniform's gain over no funding):
  12 percent at b = 0.1 falling to 5 percent at b = 1.
- Report Fig 2 (output gain vs budget scale for optimal and uniform; the narrowing gap is
  the value of targeting): banked for \S8, 200-seed regeneration still owed.
- The under-resourced option has no implemented simulation strategy; its refutation is
  analytic only (Example 2). Flag if a quantitative section ever claims measurement.
