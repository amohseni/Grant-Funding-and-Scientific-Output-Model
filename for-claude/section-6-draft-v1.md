# Section 6 (Peer review), draft v1 (2026-09-06)

Function analysis. \S6 is the paper's center of gravity: the section where the thesis is
stated as a definition and then measured. Its inputs are all in place: \S4 defined the
gap and the optimal allocation; \S5 showed that records alone leave most of the span to
complete information unrecovered and promised that review can supply the missing
capability signal "under the right conditions"; \S2 promised that the efficacy dispute
measures accuracy while our model measures value, and that no single number settles the
question. \S6 must now do five things, in order: (1) define the value of review precisely
(the licensed formulation) and say how the model measures it; (2) show the central
result, convergence: as the review signal's informativeness increases, allocation
approaches the gap rule and output approaches the complete-information benchmark, while
the records-only baseline stays flat (the D-4 figure, the paper's thesis figure);
(3) state what the value depends on, the conditional thesis: capability dispersion,
budget, and informativeness itself, with the flattening in informativeness that makes
modest review enough; (4) pay off \S2 by reading the empirical dispute inside the model
(the AUC calibration): the same measured accuracy is worth much in one field and little
in another; (5) state the conditions promised by \S5's bridge, chiefly calibrated trust
(the overtrust caveat), and hand \S7 the move from whom to when. Numbers policy per
\S5's pattern: qualitative claims in the body, numbers in footnotes with their context.
Drafted plain-first; substitution test run on every sentence about a formal quantity.

---

## 6. Peer review: what a capability signal is worth

Peer review enters our model as \S3 specified: a noisy signal of capability, with
informativeness governed by the noise parameter $\tau_K$. The value of review, for a
Bayesian funder, is the expected value of its reduction of the funder's uncertainty
regarding the gap. Our model measures this value in output: the additional expected
output of a funder allocating with the review signal over the same funder allocating
without it.

The central result is convergence. As the review signal's informativeness increases, the
funder's allocation approaches the gap rule, and its output approaches the
complete-information benchmark: with a sharp signal, nearly the whole span of \S5 is
recovered.\footnote{At the default parameters, the correlation between the funder's
grants and the optimal grants increases from 0.22 at the noisiest signal in our range to
0.95 at the sharpest; under heavy-tailed capability ($\alpha_K = 1.3$), from 0.34 to
0.97. The output gap to the complete-information benchmark falls from 17.2 to 1.7
percent of no-funding output.
% TODO: re-derive via verify_all_claims.R before lock.
} The records-only baseline of \S5, displayed alongside, stays flat: further rounds of
records do not substitute for a signal of capability, however many accumulate.

[FIGURE 4 placeholder: the convergence figure (D-4). Allocation correlation with the
optimal grants as informativeness increases, records-only baseline alongside; data in
the model repo's sweep_results/D_gap_convergence/; build in the fig1-3 pgfplots style
once the .rds data is exported to .dat from a Mac session.]

How much this convergence is worth depends on the field, in three ways. First, on the
dispersion of capability. Where capability is heavy-tailed, most of the output at stake
sits with a few researchers of very high capability, and identifying them does not
require a sharp signal: even a coarse signal captures most of review's value. Where
capability is spread more evenly, less depends on whom the funder picks, and review's
value is smaller. Second, on the budget, inherited from \S4: targeting matters less as
the budget increases, and information about whom to target inherits exactly this
comparative static. Review is worth most where money is tight. Third, on informativeness
itself: review's value is nearly constant across sharp and moderately noisy signals and
declines only beyond that, so modest informativeness suffices, and sharpening review
past that point buys little.\footnote{At the default parameters the decline begins near
$\tau_K \approx 2.5$; our sweep runs from near-perfect review ($\tau_K = 0.05$) to
nearly uninformative review ($\tau_K = 20$).}

The efficacy dispute of \S2 can now be read inside our model. Among funded NIH grants,
percentile scores predict subsequent productivity barely better than chance; in our
model, a signal of that accuracy corresponds to very noisy review, well beyond the point
where the value of further informativeness has flattened.\footnote{An AUC of 0.54
\citep{Fang2016} corresponds to $\tau_K > 20$ at the default parameters.
% TODO: verify the share-of-value language for tau at this level (heavy vs moderate
% tails) against D_gap_convergence data before lock.
} Whether review of that quality is worth having then depends on the field: under
heavy-tailed capability it still adds real value, because even a coarse signal separates
the researchers who matter most from the rest; under more evenly spread capability it
adds little. This is why no single number settles the dispute. The accuracy that both
sides measure is one input to review's value; the field's capability distribution and
the funder's budget are the others, and the same accuracy is worth a great deal in one
field and nearly nothing in another.

The value of review also presupposes that the funder holds the signal at its true
informativeness. A funder that overtrusts review, treating a noisy signal as sharp, can
do worse than one that ignores review altogether: with moderately spread capability,
overtrust turns review's value negative; under heavy tails it reduces the value without
reversing it.\footnote{With true noise $\tau_K = 3$ and the funder's belief at
$\tau_K \leq 1$: at $\alpha_K = 2$, review's measured value turns negative, consistently
across cells though within simulation noise cell by cell; at $\alpha_K = 1.3$ it
remains positive at a reduced level.
% TODO: verify before lock; state the z-caveat honestly if the numbers are kept.
} Overtrust costs more than undertrust forfeits. These are the conditions of \S5's
promise: review helps a funder that weights it correctly, and helps most where
capability is dispersed and money is tight.

Review tells the funder whom to fund. \S7 turns to when: how funding should be spread
across rounds, and how much that choice matters.

---

## Notes for Aydin

1. The definition paragraph uses your licensed formulation verbatim ("the expected value
   of its reduction of the funder's uncertainty regarding the gap") and then gives the
   operational measure in output terms. One concept, one term: "informativeness" for the
   signal property (aligned in \S5's telegraph sentence too, which said "precision";
   fixed); "sharp" and "coarse" appear only as level descriptions after informativeness
   is established. Veto if you want them out.
2. The three determinants paragraph is the conditional thesis of the abstract
   ("signal's value increases with the dispersion of capability... under heavy-tailed
   capability distributions even a coarse signal suffices") cashed out; the budget
   determinant explicitly inherits \S4's comparative static rather than being introduced
   fresh.
3. Numbers all footnoted with context per the \S5 pattern. Three TODO flags in the tex
   comments: re-derivation of D-4 numbers; verification of the share-of-value language
   at tau > 20 (the safe body claim is "still adds real value... adds little," not a
   share); the overtrust cell-level z-caveat (the handoff reports the negative values as
   consistent across cells but individually within noise; the body says "can do worse,"
   which the data supports as a pattern, and the footnote carries the honesty).
4. FIGURE 4 (convergence, the thesis figure) is a placeholder: its data lives in
   sweep_results/D_gap_convergence/ as .rds, which needs a Mac session (or R here) to
   export to .dat for the pgfplots family style. Design sketch, for your verdict before
   building: x = informativeness (noise decreasing rightward, labeled in words to avoid
   axis-direction confusion), y = correlation of the funder's grants with the optimal
   grants; one curve per tail parameter (defaults and heavy), records-only baseline as a
   flat gray reference labeled directly; no legend. Alternative y: output shortfall to
   the complete-information benchmark (percent of no-funding output), which matches the
   span language of \S5; my lean is correlation for the main panel because it shows the
   mechanism (allocation approaching the rule), with the output version in \S6's text.
   Your call.
5. "The efficacy dispute of \S2 can now be read inside our model" uses "read" in the
   literal reading sense; your ban covers "reads" standing for "accounts for." Flag if
   you want it out anyway ("interpreted inside our model").
6. The overtrust paragraph states the funder-side condition only. The other conditions
   (dispersion, budget) are in the determinants paragraph; "states the conditions" from
   \S5's bridge is thereby discharged across the two paragraphs. If you want one
   explicit conditions sentence gathering them, say the word.
7. Forward-only check run: "the span of \S5" has its antecedent; "the records-only
   baseline of \S5" likewise; no unresolved referring phrases. Substitution test run;
   "center of gravity" appears only in these notes, not the draft.
