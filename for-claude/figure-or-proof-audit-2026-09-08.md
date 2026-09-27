# Figure-or-proof audit (Aydin's rule, 2026-09-08)

Rule audited against: a substantive claim is backed by a proof in the appendix or by
simulation results illustrated in a figure; footnotes alone do not suffice. Footnotes
keep their role for generating context, caveats, and secondary magnitudes.

Determinations and actions, by section:

## \S4 (optimal allocation)
- All claims analytic. Support = proof. ACTION: gap-rule-proof-v1.md converted to
  Appendix A (appendix-gap-rule.tex): Lemmas 1-3, Proposition 1 restated and proven,
  Corollaries 1-2, Examples 1-2, Remarks 1-3. All four Appendix [X] pointers in \S4
  resolved to \ref targets (proof, two examples, Corollary 2); fig2's caption now
  refs Corollary~\ref{cor:targeting-vanishes}. The two remaining [X] pointers (\S3
  normalization, \S3 specification/parameter table) belong to the owed
  simulation-specification appendix, marked with TODO comments.
- Figures 1-2 already carry the frontier and targeting-value content. No new figures.

## \S5 (track records)
- Underdetermination argument: analytic, carried in text (identical-output pairs; the
  lambda/2A lower bound). OK as is.
- Thin-grants mechanism: Figure 3. OK.
- Cost of underdetermination (was footnote numbers only): now carried by \S6's
  records-rounds figure; \S5's display promise recast to match ("beside a funder
  holding the review signal").

## \S6 (peer review)
- Value by field, flattening, even-spread case, near-complete recovery: Figure
  (review-value), already placed.
- Records do not catch up (was footnote): NEW records-rounds figure; footnote
  removed, numbers to caption.
- Overtrust (was footnote): NEW overtrust figure with 1-SE error bars (the honest
  choice given the within-noise caveat); slim footnote keeps only the tenfold
  asymmetry, in relative terms.
- AUC calibration (C12): stays footnoted. Determination: it is a calibration of an
  external estimate into the model's noise scale, not a model result; the results it
  feeds (C12b) are carried by the review-value figure's right edge.
- Structural cause (C5b): analytic mechanism stated in text; the cumulative-Fisher
  computation is available as an appendix proposition if wanted (flagged, not built).

## \S7 (timing, draft)
- Main schedule result: NEW depth-schedule figure.
- Compounding boundary / poverty muting: NEW timing-boundary figure.
- Worth-of-scheduling calibration: NEW schedule-vs-signal figure (parity line; all
  cells below it).
- Lump-schedule identity at zero resources: analytic, in text (ledger T2 records the
  derivation).
- Planner-free early-mass check (T3): stays footnoted as robustness behind the
  analytic correction. Determination: verification, not headline.
- Channel attribution (T7): the one substantive \S7 result still footnote-only. A
  fourth figure (center of mass over the free-x-paid rate grid) would carry it;
  flagged for Aydin rather than built, to keep \S7 at three figures.

## Cross-cutting (same pass)
- All results recast in percent terms over their relevant comparisons (no-funding
  output, uniform funding's output, the even schedule's output); no publication
  counts remain in body or footnotes of \S\S5-7 drafts/tex.
- Single-number summaries of sweep-ranging quantities eliminated (the one-thirtieth
  case); replaced by range statements over the tested cells with the extreme named.
