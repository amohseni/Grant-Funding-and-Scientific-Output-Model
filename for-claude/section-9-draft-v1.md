# Section 9 (Production technology), draft v1 (2026-09-08)

SUPERSEDED 2026-09-09 by Aydin's decision: no \S9. Paragraphs 2-3 (signal robustness,
timing scoping) -> tex/appendix-technology.tex (Appendix C, with Table 2 replacing the
timing footnote); paragraph 4 (concentration) -> closing movement of \S8 in main.tex,
retitled "Spreading funds: seed grants, lotteries, and concentration"; non-monotone
clause kept at his word. Text below is historical.

Function analysis. Promises to pay: intro, "\S9 varies the production technology and
revisits the concentration question"; \S2, "how widely funds should be spread (\S9)"
and \S1's concentrate-versus-spread dispute; \S8's bridge, "Next we vary the
production technology and ask which of our results survive." The section's function
in the whole: it is the robustness audit of the paper's one deep modeling commitment
(capability and resources are complements in production) and the payoff of the
spreading dispute. Two jobs that reinforce: (1) sort the paper's results into
technology-robust (the information results: the whole of \S\S5-6 and \S8's pricing
logic) and technology-scoped (the timing results of \S7); (2) show that the
concentration question has no one-dial answer: the technology's effect on how widely
the informed funder spreads is budget-conditional and non-monotone, which replaces
the tempting story that the dispute is at bottom a disagreement about
substitutability. Every number re-derived in-session (verify_s9.R); figures 12-13
built. No new mechanisms; the section re-runs established comparisons across the
family. Claims ledger: section-9-claims.md.

---

## 9. Varying the production technology

Throughout, research output has taken one form: expected output is
$2AKR/(K + R)$, in which capability and resources are complements and output is
limited by the scarcer of the two. This section asks which of our results depend on
that choice. We embed our form in a family of production technologies indexed by an
exponent $\gamma$. At $\gamma = 0$, the Cobb-Douglas case, a shortfall in either
input can be offset by more of the other; our form sits at $\gamma = -1$; as
$\gamma$ decreases the inputs grow harder to substitute for one another; and in the
limit, the Leontief case, output depends only on the scarcer input. We re-run the
model's central comparisons across this family.\footnote{Two rounds, the greedy
allocator, one hundred seeds per cell (fifty for Leontief); the concentration runs
below are single-round with two hundred seeds. Leontief allocations are computed
with a finer allocation step, and ties at the kink are broken by fill order; the
Leontief points carry that caveat.}

The information results survive the whole family
(Figure~\ref{fig:signal-robustness}). At every technology tested, from Cobb-Douglas
to Leontief, the value of the review signal increases as capability grows more
heavy-tailed and as review grows more informative: the delineation of review's value
by field is not an artifact of our production form. What changes with the technology
is the scale. The more complementary the inputs, the more review is worth: under
heavy-tailed capability with sharp review, the signal's value increases from about
seven percent of no-funding output at Cobb-Douglas to about thirty percent at
$\gamma = -3$ and a third at Leontief.\footnote{6.8, 29.3, and 33.7 percent; the
same ordering holds at every tested tail and noise level in which the value is
distinguishable from zero.} Complementarity raises the stakes of targeting, and with
them the value of the information that targeting requires.

The timing results do not survive the family's substitutable end. At Cobb-Douglas
the value of deliberate scheduling essentially vanishes and the optimal schedule
stays within a hundredth of even; at $\gamma = -3$ both the value of scheduling and
the late lean strengthen relative to our form.\footnote{At five rounds: the
scheduling gain is at most five thousandths of one percent of no-funding output at
Cobb-Douglas across compounding rates, with the schedule's center of mass within
0.012 of even; at $\gamma = -3$ and the highest compounding rate the gain is about
one and a half percent and the center of mass 0.66.} In our model, when inputs are
substitutable the difference in output between placing a dollar well and placing it
adequately shrinks, and with it the value of waiting to learn where to place it.

The concentration question remains. Whether funding should be concentrated on a few
researchers or spread across many is an old dispute, and one might hope the
production technology settles it: the harder it is to substitute resources for
capability, the more the money should chase the few who can use it. The relation is
not that simple (Figure~\ref{fig:concentration}). With a tight budget in a
heavy-tailed field, harder complementarity does concentrate the informed funder's
grants. With an ample budget it does the opposite: as complementarity hardens, the
optimal move is to top every researcher up toward their bottleneck, and the
allocation spreads. And where capability is evenly spread, the relation is not even
monotone: concentration peaks at an intermediate complementarity and declines toward
both ends.\footnote{The peak's rise over the $\gamma = -12$ level is 3.5 standard
errors at the four-hundred-seed refinement.} The technology does not settle the
concentration question by itself; the budget and the field's capability distribution
interact with it at every point.

What, then, does the production commitment buy and cost? The value of information
and its dependence on the field are robust to the technology; the timing results
hold where capability and resources are complements and vanish where they are
substitutable; and how widely funds should be spread is set jointly by technology,
budget, and field, with no one of the three deciding it alone. Next we take stock:
what the model implies for funders, and what its assumptions leave out.

---

## Notes for Aydin (v1)

1. TITLE candidates: (i) "Varying the production technology" (used); (ii) "The
   production technology and the concentration question"; (iii) "Robustness: the
   production family". Yours.
2. FIGURE 13 REDESIGNED once during construction: absolute Gini levels hid the
   shapes (the three curves' movements are small against the 0-1 scale), so the
   figure shows the CHANGE from each curve's Cobb-Douglas level, with the absolute
   levels in the caption. Flag if you prefer absolute levels (or three small panels
   with independent scales).
3. OBJECT DISCIPLINE: the concentration runs measure the REVIEW-INFORMED funder's
   allocation (single round), not the complete-information optimal allocation; the
   text and caption say "the informed funder's grants," never "the optimal
   allocation."
4. MECHANISM SENTENCES flagged for your judgment: (a) "topping every researcher up
   toward their bottleneck" as the ample-budget equalization mechanism (the
   campaign memo's reading; direction verified, mechanism not separately isolated);
   (b) the P3 closing sentence ("the difference in output between placing a dollar
   well and placing it adequately shrinks") is an interpretive gloss of the
   measured Cobb-Douglas null, marked by "In our model"; say the word to cut either.
5. LEFT OUT deliberately: the family's seed-floor costs (the re-derivation covers
   Cobb-Douglas and gamma = -3, worst 0.44 percent of no-funding output, but the
   flooring comparison's exact strategy pair needs a column-level check before
   quoting; and the point duplicates \S8); the correlation-robustness runs (not yet
   re-derived; can add on your word). Ledger rows G7-G8 record both as available.
6. The tierA/tierB ground uses the GREEDY allocator (the sweeps predate the smooth
   allocator); stated in the generating-context footnote. Within-family comparisons
   only; no cross-quoting against smooth-run magnitudes.
7. Bridge assumes \S10 = discussion (implications + limitations) per the intro.
