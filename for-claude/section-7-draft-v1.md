# Section 7 (Timing), draft v1 (2026-09-08)

STATUS: SUPERSEDED BY tex/main.tex \S7 (typeset 2026-09-08 on Aydin's approval of v2;
"no more figures needed" = T7 attribution figure declined, decomposition stays
footnoted). One wording change at typeset: the boundary footnote's "half a percent of
the budget's center of mass" corrected to "the schedule's center of mass is within
0.005 of even" (a center of mass is a number, not a budget share). Edit main.tex from
here on.

Function analysis. The intro promises: "\S7 examines the timing of funding, including
the strongest argument for spending early." \S6 hands off with "Next we turn to when to
fund: how funding should be spread across rounds, and how much that choice matters."
\S5 promised: "\S7 measures this [deep grants recovering the value of discrimination]
and its consequences for timing." So \S7 must do five things, in order: (1) re-introduce
the round structure and the two growth facts (capabilities compound through research;
grants are consumed) and pose the schedule question; (2) state the strongest argument
for spending early (resource-poor fields) fairly, then give the correction: the argument
establishes a complementarity, not early mass; (3) state the main result: observation
comes first, money after; spending late is rewarded wherever capabilities compound
appreciably, and poverty mutes but does not reverse this; (4) pay off \S5: thin grants
reveal nothing a funder can exploit later, so paying for information requires depth, and
at depth the informed funder recovers most of the value of discrimination while keeping
a below-even early share; plus the attribution (the free growth channel drives late
spending); (5) calibrate honestly: what the schedule choice is worth, which is little
next to the review signal, growing with the growth rate and horizon, saturating beyond
short horizons, and vanishing under Cobb-Douglas production. Bridge to \S8 (seed grants
and lotteries) in "Next we" form. Do-not-claim constraints observed: no strict
no-front-loading claim (min b_idx 0.497); pure-paid front-loading is real but mild;
no S8-S5 < 0 depth cells in the body (CE-mispricing caveat would be owed); every
planning claim scoped to complementary production. All numbers below re-derived
2026-09-08 from the canonical sweep files (verify_s7.R + OUTPUT in this folder).

---

## 7. The timing of funding

Funding in our model is spread over rounds. The funder holds a fixed budget for the
whole horizon, and between rounds two things move: researchers produce output, which
the funder observes, and capabilities compound through the research that output
represents. Grants themselves are consumed within their round. Up to now the funder has
spent its budget in equal parts each round; this section asks whether it should not:
whether to spend early, evenly, or late, and how much the choice matters.

The strongest argument for spending early concerns resource-poor fields. Researchers
without resources produce nothing; where nothing is produced, nothing is observed and
nothing compounds; so money held for later rounds buys neither evidence nor growth in
the meantime. On this argument a funder facing a poor field should spend at once, to
start the engine. The argument fails on inspection. A funder that spends everything in
the first round allocates before observing any output, exactly as blind as one that
spends everything in the last; in the extreme case of a field with no resources at all,
the two schedules produce identical expected output. What the argument does establish
is a complementarity: early spending enables observation and growth, and late spending
exploits what early spending revealed. That favors putting some money early. It does
not favor putting most of it early.\footnote{In a planner-free check over a grid of
two-block schedules executed by a funder re-deciding each round, the best schedule that
puts an above-even share in the early rounds loses to the even schedule, by 1.6 percent
of the even schedule's output at negligible compounding and by 3.5 percent at the
default compounding rate (six rounds, resource-poor community, two hundred seeds).}

Our simulations bear this out and say which way the asymmetry runs: observation comes
first, money after. Wherever capabilities compound at more than a token rate, the
optimal schedule holds a below-even share in the first round and deploys the bulk once
output has revealed who is capable; spending late is rewarded because late money is
informed money, and because early output compounds on its own
(Figure~\ref{fig:timing-boundary}). Resource poverty raises the early share and
decreases the intensity of the effect, but does not reverse
it.\footnote{Figure~\ref{fig:timing-boundary} shows three of the five baseline
resource scales in the map; the other two lie between the curves shown. Every cell
whose schedule leans early has a compounding rate $\epsilon \leq 0.03$, and there by
at most half a percent of the budget's center of mass.}

[FIGURE: tex/fig8-timing-boundary.tex + fig8-boundary.dat. Schedule center of mass vs
compounding rate (log x), three resource scales, dashed even line. Carries T4.]

How much a funder can learn from its own grants depends on their size. A researcher on
a thin grant produces output nearly proportional to the grant and nearly independent of
capability, so a funder whose grants are thin learns almost nothing it can exploit
later: in a resource-poor field with no review signal, such a funder gains almost
nothing over uniform funding. Deep grants change this. As grants approach the scale of
capability itself, output separates the capable from the rest, and the funder recovers
most of the value of discrimination, on a schedule that puts observation before money:
a small early share to make output observable, the mass deployed once informed
(Figure~\ref{fig:depth-schedule}).\footnote{In a resource-poor community with
heavy-tailed capability, no review signal, and six rounds (two hundred seeds): the
informed funder's gain over uniform funding increases from a third of one percent of
uniform funding's output at standard depth to nine percent at six times that depth and
ten percent at twelve times, roughly the gain a sharp review signal delivers at
standard depth (nine to eleven percent). In these runs the budget is scaled against
capability rather than resources: at the six-times depth, the average grant per
researcher per round is half of mean capability.}

[FIGURE: tex/fig7-depth-schedule.tex + fig7-schedule.dat. Round shares at standard vs
deep funding, T=6 poverty corner; even at standard depth, below-even first round then
informed mass at depth. Carries T6's schedule half; the gain numbers stay footnoted.]

The pattern has a single main cause. Capabilities compound through all research output,
funded or not; grants add to output and so add to compounding. Separating the two
channels: with only the grant-fed channel active, the optimal schedule leans mildly
early, since earlier grants compound longer. With only the free channel active, it
leans strongly late. The free channel dominates once it reaches even a fraction of the
grant-fed rate, so in any field where research compounds regardless of who is funded,
the funder's reason to wait outweighs its reason to hurry.\footnote{With the free
channel off, the schedule's center of mass is 0.48 to 0.49 across grant-fed rates,
strictly below even; with the grant-fed channel off, 0.53 to 0.68. The free channel's
rate takes over at roughly one eighth to one third of the grant-fed rate (two hundred
seeds, five rounds).}

Finally, what is the schedule choice worth? Less than information, in every setting we
test. Across all horizons and compounding rates in our sweeps, the gain from deliberate
scheduling over even installments re-decided each round stays below the review signal's
value: at the default compounding rate and below it is smaller by an order of magnitude
or more, and it approaches the signal's value only at the joint extreme of the highest
compounding rate and the longest horizon we consider
(Figure~\ref{fig:schedule-vs-signal}).\footnote{As a fraction of no-funding output,
the scheduling gain ranges from indistinguishable from zero to 0.7 percent across the
thirty-two tested cells; the review signal's value in the same cells ranges from 0.8
to 7.5 percent.} The gain increases roughly with the square of the horizon over short
horizons and saturates beyond them.\footnote{Fitted exponents 2.3 and 2.1 for horizons
up to five rounds at $\epsilon = 0.3$ and $0.85$; 1.4 and 0.8 over five to ten
rounds.} And it exists only where capability and resources are complements, as in our
production form: under Cobb-Douglas production the value of scheduling essentially
vanishes.\footnote{With Cobb-Douglas production, the gain at five rounds is at most
five thousandths of one percent of no-funding output across compounding rates, and the
schedule stays within 0.012 of even; at a complementarity midway to Leontief the gain
is about one and a half percent of no-funding output and the schedule's center of mass
0.66.} For a funder thinking about how to maximize their impact, choosing whom to fund
outweighs choosing when.

[FIGURE: tex/fig9-schedule-vs-signal.tex. Scheduling gain as a share of the review
signal's value vs horizon, five compounding rates, dashed parity line. Carries T8-T9.]

Next we price two schemes that recur in policy debates: uniform seed grants and
lotteries.

---

## Notes for Aydin (v1)

1. DESIGN as approved in chat (your "let's setup for section 7"): objection-first for
   paragraph 2 because the intro promises "the strongest argument for spending early";
   everything after is results-first. Say the word if you want the main result before
   the objection.
2. TITLE candidates, ranked: (i) "The timing of funding" (used); (ii) "The timing of
   funding: money follows information" (carries the thesis but reads as a slogan);
   (iii) "When to fund". Word-level, yours.
3. NAMING decision: the schedule pattern (small early share, mass once informed) is
   described plainly in the body. The campaign memo's name for it is "seed-and-harvest";
   I did not introduce the coinage. If you want a name for repeated reference in \S8+
   (the seed-grant section will gesture back at it), "seed-and-harvest" is the
   candidate; otherwise the plain description repeats.
4. VERIFICATION: every number re-derived this session from the canonical sweep files
   (verify_s7.R, verify_s7_OUTPUT.txt); all match the integration handoff. One
   refinement: the handoff's "timing is small next to the signal" is stated
   $\epsilon$-dependently in the draft (one-thirtieth at $\epsilon = 0.3$, one quarter
   at $0.85$), since the flat "1/30" holds only at the default rate.
5. MODAL discipline applied: (a) "wherever capabilities compound at more than a token
   rate" carries the $\epsilon^* \approx 0.02$ boundary; no strict "never front-loads"
   (min center of mass 0.497, economically nil); (b) the pure-paid mild front-loading
   is in the attribution paragraph at its measured size; (c) no S8$-$S5 depth cells
   (the deliberate planner's loss to re-deciding at depth is real but needs the
   CE-mispricing caveat; it is left out of the body entirely, say the word if you want
   it footnoted); (d) Cobb-Douglas scope stated in body, not fine print.
6. TERMINOLOGY introduced: "schedule" (how the budget is spread across rounds),
   "center of mass" only in footnotes (b_idx gloss: 0.5 = even), "grant-fed" vs "free"
   compounding channels (candidate locks; the memo says "paid/free" but "paid"
   collides with paying for information one paragraph earlier). "Depth"/"thin grants"
   reused from \S5 as locked.
7. The "identical expected output" claim for the two lump schedules at zero resources
   is analytic (with no resources, unfunded researchers produce nothing, so nothing is
   observed or compounds between rounds either way; both lumps are allocated on the
   prior); ledger row T2 records the derivation.
8. RESOLVED (v2): the depth runs scale the budget against capability, not resources
   (budget_ref = "K": B = b n E[K]); the footnote now says so plainly and gives the
   per-researcher gloss (at b = 3, average grant per round = half of mean capability).
   My earlier "fraction of baseline resources" gloss was wrong for these cells.
9. RESOLVED (v2): three \S7 figures BUILT under your figure-or-proof rule:
   fig7-depth-schedule (round shares, standard vs deep; the b = 3 schedule re-derived
   in-session, 24 seeds, SEs < 0.006, round-1 share 0.104 vs canonical 200-seed
   0.110), fig8-timing-boundary (center of mass vs compounding, three resource
   scales), fig9-schedule-vs-signal (scheduling gain / signal value vs horizon, five
   compounding rates, parity line). Placeholders mark their positions.
10. v2 (your 2026-09-08 directives): both your sentences adopted verbatim ("Resource
    poverty raises the early share and decreases the intensity of the effect...";
    "For a funder thinking about how to maximize their impact..."). All results
    recast in percent terms over their comparisons (even schedule's output, uniform
    funding's output, no-funding output); no publication counts remain. The
    one-thirtieth sentence replaced by the sweep-ranging statement (below the
    signal's value in all 32 cells; order of magnitude or more at default rate and
    below; approaches parity only at the joint extreme), carried by fig9. Remaining
    footnote-only result: the attribution paragraph's channel decomposition (T7);
    a fourth figure (center of mass over the eps_free x eps_paid grid) would carry
    it. Flagged rather than built to keep \S7 at three figures; say the word.
