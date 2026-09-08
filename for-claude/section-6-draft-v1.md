# Section 6 (Peer review), draft v4 (2026-09-08; v4 = full claims verification pass)

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
signal of each researcher's capability. The value of review, for a Bayesian funder, is the expected value of
its reduction of the funder's uncertainty regarding the capability-resource gap. We can
quantify this value in terms of output: the additional expected output of a funder
allocating with the review signal over the same funder allocating without it.

The direction of review's effect is what one would expect: a more informative signal
brings the funder's allocation closer to the optimal allocation and its output closer to
what complete information would deliver, and with a sharp enough signal the funder
recovers most of the output that a funder relying on records alone leaves
unrealized.\footnote{At the default parameters, the correlation between the funder's
grants and the optimal grants increases from 0.22 at the noisiest signal in our range to
0.94 at the sharpest; under heavy-tailed capability ($\alpha_K = 1.3$), from 0.34 to
0.96. The output shortfall to the complete-information benchmark decreases from 10.4 to
1.6 percent of no-funding output; under heavy-tailed capability, from 17.2 to 1.7
percent. The increase in review's value is not strictly monotone at the sharpest
signals: under heavy-tailed capability, its value at $\tau_K = 0.05$ is about two
percent below its value at $\tau_K = 0.3$, a small but statistically solid reversal
(paired across seeds, $z = 4.5$).} Two further effects are less obvious: (1) a funder
relying on records alone does not catch up, and (2) what a review signal is worth
varies widely with the field. The first effect has a structural cause. Capabilities
compound across rounds while resources do not accumulate, so each researcher's output
approaches the level their resources sustain, and output at that level carries less and
less information about capability. Further rounds of records do improve the funder's
allocation, but slowly: after twenty rounds of records the funder's allocation remains
further from the optimal allocation than a funder with the review signal is after
one.\footnote{At the default parameters, the correlation between the records-only
funder's grants and the optimal grants increases from 0.08 in the first round to 0.41
by round twenty; a funder with the default review signal ($\tau_K = 1$) is at 0.81 in
the first round. Under heavy-tailed capability: 0.07 to 0.43, against 0.92. Fifty
seeds.}

[FIGURE 4: tex/fig4-review-value.tex + fig4-value.dat. Share of the records-only
shortfall that review recovers, versus signal noise (log scale, tau 0.05 to 20), three
capability distributions (alpha_K = 1.3, 2, 3.5). Fine-grid data: 17 noise levels x 3
distributions x 50 seeds, run in-session in R, validated bit-identical to the canonical
D_gap_convergence cells. Placed here; carries C3, C4, C8, C10, and the main result.]

Our main result delineates how much review is worth as a function of the field. Review
is worth most where capability is heavy-tailed and the budget is tight, and in the
heavy-tailed case even a signal of moderate informativeness captures most of its value.
The reason is twofold. Where capability is heavy-tailed, most of the output a funder can
add comes from getting resources to the few researchers of unusually high capability.
And precisely those researchers are the easiest to identify: they stand far apart from
everyone else, so even a noisy signal separates them. Identifying the most capable
researchers captures most of the value and requires the least information. Sharpening
review beyond a modest level therefore adds little.\footnote{At three times the default
signal noise, review retains 82 percent of its maximum value under heavy-tailed
capability ($\alpha_K = 1.3$), 44 percent at the default distribution ($\alpha_K = 2$),
and 10 percent where capability is spread more evenly ($\alpha_K = 3.5$). The noise
level at which review's value falls to half its maximum is roughly $\tau_K = 10$,
$2.5$, and $0.75$ in the three cases; the default noise is $\tau_K = 1$, and our sweep
runs from near-perfect review ($\tau_K = 0.05$) to nearly uninformative review
($\tau_K = 20$). Fifty seeds per cell.}

A corollary concerns where review effort goes. In our model, what pays is separating
the exceptional researchers from the rest, and that separation requires little
information; fine distinctions among the many applications of comparable middling merit
add comparatively little output. A review system that invests most of its effort in
exactly those fine distinctions may be misallocating it.

The remaining comparisons run the other way. Where capability is spread more evenly, no
researcher matters much more than another, so review of any informativeness adds
little.\footnote{At $\alpha_K = 3.5$, the entire records-only shortfall is 3.5 percent
of no-funding output (24.6 percent under heavy-tailed capability), and even
near-perfect review recovers only 72 percent of it: review's value at its maximum is
2.5 percent of no-funding output, against 23 percent under heavy-tailed capability.}
And where the budget is ample, optimal and uniform funding produce nearly the same
output, so information about whom to fund has little left to add: review is worth most
where money is tight.

These results bear on the dispute over the efficacy of peer review. Among funded NIH
grants, percentile scores predict subsequent productivity barely better than chance. In
our model, a signal of that accuracy is very noisy review.\footnote{An AUC of 0.54
\citep{Fang2016} corresponds to $\tau_K \gtrsim 20$ at the default parameters, the
noisiest signal in our sweep. The mapping from measured accuracy to signal noise itself
depends on the field: at the same noise, accuracy is higher under heavy-tailed
capability, because the most capable researchers stand further apart.} Whether such
review is worth having depends on the field. At the noisiest signal in our sweep,
review still recovers about a third of the records-only shortfall under heavy-tailed
capability, under a tenth at the default distribution, and essentially none where
capability is spread more evenly.\footnote{31, 9, and 0.1 percent of the records-only
shortfall at $\tau_K = 20$, for $\alpha_K = 1.3$, $2$, and $3.5$; fifty seeds per
cell.} Accuracy is one input to review's value; the field's capability distribution and
the funder's budget are the others. No single number settles the dispute: review of the
same informativeness is worth several times more in one field than in another, and the
same measured accuracy corresponds to different informativeness in different fields.

One condition qualifies all of this: the funder must weight review at its true
informativeness. A funder that overtrusts review, treating a noisy signal as sharp, can
do worse than one that ignores review altogether: with moderately spread capability,
overtrust turns review's value negative; with heavy-tailed capability it reduces the
value without reversing it.\footnote{With true noise $\tau_K = 3$ and the funder's
belief at $\tau_K \leq 1$: at $\alpha_K = 2$, review's measured value is between $-0.5$
and $-0.1$ output units across belief levels, each cell within two standard errors of
zero but all negative, against a calibrated value of $3.4$; at $\alpha_K = 1.3$ it
remains positive at a reduced level ($14.0$ against a calibrated $17.2$). The
asymmetry: at the default distribution, a funder that overtrusts tenfold loses more
than its full calibrated value, while a funder that undertrusts tenfold retains over
half ($4.2$ of $7.6$). Two hundred seeds per cell.} Overtrust costs more than
undertrust forfeits.

Review can inform whom to fund. The next section turns to when to fund. In particular,
we examine how funding should be spread across rounds, and how much that choice
matters.

---

## Notes for Aydin (v4, full verification pass, 2026-09-08)

All numbers now re-derived. Container R reproduced the canonical D-4 cells
bit-identically (max diff 1e-14), so I ran the owed re-derivations here plus three
extensions: a fine tau grid (17 levels), the alpha_K = 3.5 row, and multi-round
trajectories (T = 8 and T = 20). Verification outcomes and the decisions they forced:

1. TWO NUMBER ERRORS CORRECTED in the convergence footnote. (a) Correlations at the
   sharpest signal are 0.94 / 0.96, not 0.95 / 0.97 (those were the sweep maxima,
   which occur at tau = 0.3, not at the sharpest). (b) "17.2 to 1.7 percent" is the
   heavy-tailed case; the defaults are 10.4 to 1.6. The old footnote attributed
   17.2/1.7 to defaults. Both fixed. "Nearly all" demoted to "most" (recovery at the
   sharpest signal is 86 percent at defaults, 93 heavy).
2. THE FLATNESS CLAIM WAS FALSE and is rewritten. "Further rounds of records leave the
   funder no closer" and the footnoted "0.13-0.18 across rounds" do not survive
   measurement: the records-only correlation INCREASES across rounds, slowly (0.08 to
   0.41 over twenty rounds at defaults; 0.07 to 0.43 heavy), while a review-informed
   funder is at 0.81 / 0.92 in round one. The 0.13-0.18 range was across capability
   distributions at round 1, not across rounds. The body now claims slow improvement
   that stays far behind, which is what the data shows; "does not catch up" survives
   at that scope. The Fisher-decay mechanism (C5b) still stands and still explains why
   records are weak evidence; I did not attribute the residual slow improvement to any
   mechanism in the text, because I have not isolated it (candidates: slow accumulation
   of the finite total information; growing capability dispersion making coarse ranks
   more informative). Say the word if you want the isolation run.
3. NEW MEASURED WRINKLE, footnoted: review's value is not strictly monotone in
   informativeness at the sharp end. Under heavy tails, value at tau = 0.05 is ~2
   percent BELOW tau = 0.3 (paired across seeds, z = 4.5; flat at defaults, z = 0.1).
   Small, real, and theoretically curious (a myopic funder does slightly better with a
   slightly noisy signal; possibly noise-as-exploration under compounding). Options:
   keep the footnote (current), drop it as clutter, or chase the mechanism. Your call.
4. ELBOW RESCOPED: tau ~ 2.5 was the defaults' half-value point, not a kink. The
   footnote now reports what the fine grid shows: value retained at 3x default noise
   (82 / 44 / 10 percent for alpha 1.3 / 2 / 3.5) and the half-value noise level per
   field (tau ~ 10 / 2.5 / 0.75). The half-value line is one number per field and may
   be the cleanest quantitative statement of the main result; consider promoting it
   from footnote to body.
5. ALPHA 3.5 NOW MEASURED (was owed): the even-spread case is now grounded in output
   value, not just mechanism gradients. Records-only shortfall is 3.5 percent of
   no-funding output there (vs 11.4 defaults, 24.6 heavy); review recovers at most 72
   percent of that small shortfall, and essentially none of it at high noise.
6. EFFICACY PARAGRAPH RESCOPED for a subtlety the verification exposed: the
   accuracy-to-noise mapping is itself field-dependent (AUC 0.54 maps to tau ~ 20 at
   defaults but to tau well beyond 20 under heavy tails, where separation is easier).
   So "the same accuracy is worth a great deal in one field" was not licensed as
   stated. The paragraph now makes the licensed pair of claims: same informativeness,
   several times the value; same accuracy, different informativeness. The closing
   sentence changed accordingly.
7. OVERTRUST VERIFIED (C13, C14) with the honesty caveat kept: the three negative
   cells at defaults are each within 2 SE of zero but consistently negative; heavy
   tails positive at a reduced level (14.0 vs 17.2 calibrated). The asymmetry sentence
   is grounded: tenfold overtrust loses more than the full calibrated value; tenfold
   undertrust retains over half. Two hundred seeds per cell.
8. FIGURE 4 BUILT (tex/fig4-review-value.tex): share of the records-only shortfall
   review recovers vs noise (log scale), three fields. One frame carries C3, C4, C8,
   C10, and the main result. Ranked alternatives if you want a different y: (ii)
   absolute value as percent of no-funding output (shows stakes, hides the flattening
   comparison); (iii) allocation correlation (process, not value; my old lean, now
   demoted since the section's thesis is about value). No AUC marker on the figure:
   an "AUC 0.54" line would sit at different tau per curve (see note 6), so marking
   one x-position would be wrong for two of the three curves.
9. Provenance: all new runs use the staged model.R (hash 21b0d9a), scripts and outputs
   written back to for-claude/ (run_D4_fine.R, gap_convergence_fine.csv,
   verify_s6_rounds.R + output); bit-identity against canonical cells validated before
   any new cell was trusted. The corollary's uniform-tau caveat (v2 note 2) stands.
