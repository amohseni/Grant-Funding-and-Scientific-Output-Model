# Section 6 (Peer review), draft v2 (2026-09-06)

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

Peer review enters our model as a second source of evidence about researchers: a noisy
signal of each researcher's capability, with informativeness governed by the noise
parameter $\tau_K$. The value of review, for a Bayesian funder, is the expected value of
its reduction of the funder's uncertainty regarding the capability-resource gap. We can
quantify this value in terms of output: the additional expected output of a funder
allocating with the review signal over the same funder allocating without it.

The direction of review's effect is what one would expect: a more informative signal
brings the funder's allocation closer to the optimal allocation and its output closer to
what complete information would deliver, and with a sharp enough signal the funder
recovers nearly all of the output that a funder relying on records alone leaves
unrealized.\footnote{At the default parameters, the correlation between the funder's
grants and the optimal grants increases from 0.22 at the noisiest signal in our range to
0.95 at the sharpest; under heavy-tailed capability ($\alpha_K = 1.3$), from 0.34 to
0.97. The output shortfall to the complete-information benchmark falls from 17.2 to 1.7
percent of no-funding output.
% TODO: re-derive via verify_all_claims.R before lock.
} Two further effects are less obvious: (1) a funder relying on records alone does not
catch up, and (2) what a review signal is worth varies widely with the field. The first
effect has a structural cause. Capabilities compound across rounds while resources do
not accumulate, so each researcher's output approaches the level their resources
sustain, and output at that level carries less and less information about capability;
in our simulations, further rounds of records leave the funder's allocation no closer
to the optimal allocation.\footnote{The correlation between its grants and the optimal
grants stays between 0.13 and 0.18 across rounds at the default parameters.}

[FIGURE 4 placeholder: allocation correlation with the optimal grants as informativeness
increases, records-only funder alongside; data in sweep_results/D_gap_convergence/;
build in the established pgfplots style once the .rds data is exported to .dat.]

Our main result delineates how much review is worth as a function of the field. Review
is worth most where capability is heavy-tailed and the budget is tight, and in the
heavy-tailed case most of its value comes from a signal of only moderate informativeness.
The reason is twofold. Where capability is heavy-tailed, most of the output a funder can
add comes from getting resources to the few researchers of unusually high capability.
And precisely those researchers are the easiest to identify: they stand far apart from
everyone else, so even a noisy signal separates them. Identifying the most capable
researchers captures most of the value and requires the least information. Sharpening
review beyond a modest level therefore adds little.\footnote{At the default parameters,
review's value is nearly constant for noise below $\tau_K \approx 2.5$ and declines
beyond it; our sweep runs from near-perfect review ($\tau_K = 0.05$) to nearly
uninformative review ($\tau_K = 20$).}

A corollary concerns where review effort goes. In our model, what pays is separating
the exceptional researchers from the rest, and that separation requires little
information; fine distinctions among the many applications of comparable middling merit
add comparatively little output. A review system that invests most of its effort in
exactly those fine distinctions may be misallocating it.

The remaining comparisons run the other way. Where capability is spread more evenly, no
researcher matters much more than another, so review of any informativeness adds little.
And where the budget is ample, optimal and uniform funding produce nearly the same
output, so information about whom to fund has little left to add: review is worth most
where money is tight.

These results bear on the dispute over the efficacy of peer review. Among funded NIH
grants, percentile scores predict subsequent productivity barely better than chance. In
our model, a signal of that accuracy is very noisy review.\footnote{An AUC of 0.54
\citep{Fang2016} corresponds to $\tau_K > 20$ at the default parameters.
% TODO: verify the added-value language for tau at this level (heavy vs moderate tails)
% against D_gap_convergence data before lock.
} Whether such review is worth having depends on the field: where capability is
heavy-tailed, it still adds real value, because even a noisy signal separates the most
capable researchers from the rest; where capability is spread evenly, it adds little.
Accuracy is one input to review's value; the field's capability distribution and the
funder's budget are the others. No single number settles the dispute, because the same
accuracy is worth a great deal in one field and nearly nothing in another.

One condition qualifies all of this: the funder must weight review at its true
informativeness. A funder that overtrusts review, treating a noisy signal as sharp, can
do worse than one that ignores review altogether: with moderately spread capability,
overtrust turns review's value negative; with heavy-tailed capability it reduces the
value without reversing it.\footnote{With true noise $\tau_K = 3$ and the funder's
belief at $\tau_K \leq 1$: at $\alpha_K = 2$, review's measured value turns negative,
consistently across cells though within simulation noise cell by cell; at
$\alpha_K = 1.3$ it remains positive at a reduced level.
% TODO: verify before lock; state the noise caveat honestly if the numbers are kept.
} Overtrust costs more than undertrust forfeits.

Review can inform whom to fund. The next section turns to when to fund. In particular,
we examine how funding should be spread across rounds, and how much that choice
matters.

---

## Notes for Aydin (v2, after your restructure)

1. MAIN RESULT RECENTERED per your direction: convergence demoted to "the direction of
   review's effect is what one would expect," kept only to explain the logic and to
   carry the two non-obvious facts (records never catch up; value varies with the
   field). The main result is now the field delineation, led by the heavy-tailed case
   with your two mechanisms (most addable output sits with the few most capable; they
   are the easiest to identify), and your corollary (fine distinctions among borderline
   applications add little; effort spent there may be misallocated).
2. HONESTY FLAG on the corollary: the model's informativeness parameter is uniform
   across applicants, so the corollary reads the flattening in informativeness through
   the identify-the-whales mechanism; a rank-local version (sharper discernment only
   among mid-ranked applications) is not separately simulated. The text stays within
   "in our model" and "may be"; if you want the corollary stronger, a targeted
   simulation would be needed.
3. CROSS-REFERENCES REMOVED per your rule: no back-references to earlier sections by
   number; each borrowed fact is restated in place ("a funder relying on records alone
   leaves unrealized," "optimal and uniform funding produce nearly the same output").
   The one remaining pointer is the forward "The next section," which needs no memory.
   Applied to this section only; the typeset sections' forward pointers (\S8, \S9 in
   \S4; \S6, \S7 in \S5) stand for now; say the word and I convert them to the
   self-contained style too.
4. Your "long-tailed" mapped to the locked term "heavy-tailed" (one concept, one term;
   flag if you want the lock changed instead).
5. Numbers remain footnoted with context; the three TODO flags stand (D-4 re-derivation;
   added-value language at tau > 20; overtrust cell-level noise caveat).
6. Figure 4 design still open: y = allocation correlation (my lean) vs output shortfall;
   .rds -> .dat export needs a Mac session either way.
7. Two candidate skill rules from your last two messages, logged in state.md, promoted
   on your word: (a) re-introduce central objects at section boundaries; (b) minimize
   cross-references by number; restate borrowed facts so sections are self-contained.
