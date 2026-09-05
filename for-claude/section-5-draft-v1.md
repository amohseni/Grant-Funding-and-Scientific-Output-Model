# Section 5 (Track records), draft v1 (2026-09-05)

Function analysis. \S4 closed: "We turn now to how such estimation can work: from the
track record alone (\S5)." So \S5 must deliver: (1) the analytic reason output cannot
separate capability from resources (cross-sectional underdetermination, stated exactly);
(2) the measured cost of relying on records alone (the stake \S6 will price); (3) the
thin-grants mechanism, which says where the record is least informative and what
identification requires; (4) the bridge to review as the missing capability signal.
Drafted plain-first; every formal gloss checked by the substitution test.

---

## 5. What the track record can and cannot supply

The funder observes track records: each researcher's realized output, round by round. A
track record is evidence about the gap, and the question for this section is how far that
evidence goes.

Consider first what output can identify in principle. Expected output is one number,
produced by capability and resources jointly: every pair of capability and resources with
$2AKR/(K+R) = \lambda$ produces the same expected output $\lambda$. A researcher of high
capability and modest resources and a researcher of modest capability and abundant
resources can therefore produce identical records. Output does carry some information
about capability, a lower bound: since expected output never exceeds $2AK$, a researcher
who produces $\lambda$ has capability above $\lambda/2A$. Beyond that bound, the record
cannot separate the two researchers above, and these are exactly the researchers the gap
rule treats most differently: the first may have the largest gap in the population, the
second may have none. Cross-sectional output, however carefully recorded, underdetermines
what the funder most needs to know.

Our model measures what this underdetermination costs. A Bayesian funder relying on
records alone allocates far from the gap rule and remains far: across rounds, the
correlation between its grants and the optimal grants stays between 0.13 and 0.18 at the
default parameters, and its output falls short of the complete-information benchmark by
17.2 percent of no-funding output. This shortfall is the span from records alone to
complete information. \S6 asks how much of that span a review signal recovers.

The record is least informative exactly where the gap is largest. When a researcher's
resources are small relative to their capability, expected output is approximately
$2AR$: proportional to resources, nearly independent of capability. On a thin grant, a
researcher of capability 1 and a researcher of capability 100 produce nearly
indistinguishable output. Small grants therefore buy almost no information about who is
capable; separating capability from resources by observing output requires funding depth,
grants comparable to capability itself. In our simulations of a resource-poor community
with no review signal, a funder allocating at standard depth gains almost nothing over
uniform funding, while the same funder with deep grants recovers most of the value of
discrimination (\S7 measures this and its consequences for timing). For nascent and
resource-poor fields the implication is direct: thin grants cannot reveal who is capable.

What the record cannot supply, then, is a signal of capability separate from resources.
That is what peer review can provide. \S6 measures what such a signal is worth.

---

## Notes for Aydin

1. NUMBERS AND THEIR CONTEXT: 0.13-0.18 (allocation correlation, records-only, across
   rounds at defaults) and 17.2 percent of no-funding output (gap to the
   complete-information oracle) are from D-4 (PAPER_INTEGRATION_HANDOFF, canonical); the
   defaults context is supplied by \S3's parameter footnote. Both owed re-derivation via
   verify_all_claims.R (Mac session; no R here) before lock. Do-not-claim check: the
   17.2 percent is described as the span from records alone to complete information,
   never as what review is worth (the licensed description).
2. STRATEGY-MAPPING FLAG, needs verification before tex: the plateau numbers are the
   handoff's "pubs-only baseline." Our \S3 table's Bayesian without-review strategies
   use track record + resource signal. D-3 shows the resource signal contributes
   0.3-1.1 percent everywhere (redundant), but the correlation numbers may differ
   between pubs-only and records+resource-signal funders. The draft says "records
   alone"; if the data is for pubs-only S4 strictly, either match the strategy named in
   \S3 or verify the numbers for the without-review strategy.
3. The thin-grants paragraph keeps the quantified corner QUALITATIVE ("gains almost
   nothing over uniform funding" / "recovers most of the value") because its numbers
   (S8-S2 = +0.5 at b = 0.5 vs +41.2 at b = 3) come from the exploration corner, whose
   parameters (no free signal, poverty, alpha_K = 1.3, T = 6) differ from the defaults;
   stating them would need the full context mid-paragraph. Options: (a) as drafted,
   numbers to \S7 where the corner is set up; (b) state them here with a
   context-supplying footnote. My lean: (a).
4. The K = 1 vs K = 100 at a thin grant example is the handoff's; "capability 1" and
   "capability 100" use model units, consistent with \S3's parameter conventions.
5. The cross-sectional qualifier in P1 is load-bearing: a record at VARIED funding
   levels does identify capability (dose-response), which is exactly the depth point in
   P3 and the E-1b identification result. "However carefully recorded" is scoped by
   "cross-sectional output" to keep both true.
6. The lower-bound observation (K > lambda/2A) is new in this draft: it makes P1's
   underdetermination claim exact rather than rhetorical (the record is not worthless,
   it bounds; it cannot separate). Analytic, from Lemma 1(iii)'s bound.
7. Substitution test run on every sentence; "stays between 0.13 and 0.18" replaces
   "plateaus" (no figurative plateau).
8. Terminology check: "track record" for realized output (locked); "records alone" as
   the compressed form after first use; "review signal" per \S3; "complete-information
   benchmark" aligns with \S4's "complete information."
