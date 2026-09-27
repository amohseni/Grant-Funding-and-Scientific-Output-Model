# Section 8 (Seed grants and lotteries), draft v1 (2026-09-08)

STATUS: SUPERSEDED BY tex/main.tex \S8 (typeset 2026-09-08 on Aydin's approval;
terminology aligned to "equal division" per his call, instances distinguished by
pool; Appendix B input after Appendix A). Edit main.tex from here on.

Function analysis: section-8-prep.md (approved structure; your four decisions applied:
single-round computation as ground; screened-lottery design with the same-pool split
displayed; Jensen in the appendix with only plain explanation in the body; funded
share stated, not optimized). Claims ledger: section-8-claims.md. Figures built:
fig10-floor-cost, fig11-lottery-prices. Appendix B written: tex/appendix-lottery.tex
(Proposition 2 + proof + scope remark), to be \input at typeset. Numbers policy:
percents and shares over named comparisons; ranges over swept settings.

---

## 8. Seed grants and lotteries

Two funding schemes recur in policy debates because they spread money rather than
target it. A uniform seed scheme grants every researcher an equal share of the budget
regardless of records or reviews; in practice it usually appears as a floor, with a
fraction of the budget spread equally and the remainder targeted. A lottery funds a
subset of researchers chosen by chance; as funders run it, a screen first admits
applicants to a pool of the fundable, and chance picks winners within the pool. Our
model prices what such schemes forgo: the expected output a funder gives up by not
targeting the capability-resource gap. What they save, the costs of running review
and of the proposal competition it induces, lies outside our model and has been
priced by others;\footnote{\citet{GrossBergstrom2019} model grant competition as a
contest in which scientists sink costly effort into proposals, show that at low
paylines this effort can rival the scientific value of the funded research, and show
that a partial lottery above a threshold decouples the waste from the payline.} at
the close we set the two sides of the ledger together.

Both schemes pay, in different amounts, one price: the value of targeting, the
additional expected output of the optimal allocation over spreading the same budget.
How much that is depends on the field. For a funder allocating with records and
review, targeting adds more than uniform funding's entire gain over no funding where
capability is heavy-tailed, review sharp, and the budget tight, and about a quarter
of that gain in the default field at the largest budget we sweep.\footnote{Over two
rounds: the review-informed funder's gain over uniform funding, as a fraction of
uniform funding's gain over no funding, decreases from 1.17 to 0.72 across budget
scales $b = 0.1$ to $1$ with heavy-tailed capability and sharp review
($\tau_K = 0.3$), and from 0.74 to 0.25 in the default field with the default signal;
the two regimes differ in both the capability distribution and the signal's noise.
Two hundred seeds per cell.} The conditions that make review worth little, capability
spread evenly, budgets ample, signals uninformative, are the conditions that make
forgoing targeting cheap.

Consider first the seed floor. Its cost increases faster than the floored share, and
its scale is set by what targeting is worth in the field
(Figure~\ref{fig:floor-cost}): flooring three quarters of a review-informed funder's
budget costs about six percent of no-funding output where capability is heavy-tailed
and review sharp, and about one percent in the default field.\footnote{Across the
swept cells, the cost is approximately the value of targeting forgone on the floored
share: the informed funder's gain over uniform funding, times the floored share,
times a factor between 0.10 and 0.86, with the factor above 0.3 at tight budgets and
heavy tails. Two hundred seeds per cell.} At a full floor the scheme is uniform
funding exactly, so the cost approaches the full value of targeting. The cost is a
price, not a prohibition; whether a floor's purpose justifies it is a judgment the
model does not make.

A lottery adds chance to spreading, and chance itself has a price. A lottery gives
each pool member the equal division's grant on average, but as a gamble: full funding
or nothing. Because a researcher's expected output increases with resources at a
diminishing rate, the gamble produces less expected output than a certain grant of
the same average size, so dividing a pool's money equally outperforms every lottery
over that pool, and a lottery over the whole population produces less than uniform
funding (Proposition~2, Appendix~\ref{app:lottery}). The measured prices are in
Figure~\ref{fig:lottery-prices}, which evaluates the schemes funders run against the
optimal allocation in a single round: an equal division among the top fifth of
researchers ranked by a screen of the review signal's form; the same division over a
doubled pool; the screened lottery, which funds half of that doubled pool at random;
uniform funding; and a lottery over the whole population. Where capability is
heavy-tailed and the budget tight, the equal division among the screen's top fifth captures about three quarters of the
optimal allocation's gain when the screen is sharp and half when it is very noisy; the screened lottery captures less, and most of its shortfall comes
from discarding the screen's ranking within the pool rather than from the gamble
itself;\footnote{At screen noise $\tau = 1$: the division among the top fifth captures
0.75 of the optimal gain, the division among the doubled pool 0.61, and the screened
lottery 0.56. Of the 0.18 the screened lottery gives up relative to the division
among the top fifth, discarding the ranking accounts for 0.14 and the gamble for 0.04. Single
round, a fifth of the population funded, two thousand populations per point.} the
full lottery captures least, less even than uniform funding. Where capability is
spread evenly and the budget ample, the ordering inverts: spreading itself is what
pays, uniform funding captures more than four fifths of the optimal gain, and every
scheme that concentrates grants, whether by screen or by chance, captures less.

When, then, is a lottery defensible on the output side of the ledger? Where what
chance discards is cheap. Replacing the screen's ranking with chance within the pool
costs a few hundredths of the optimal allocation's gain where capability is evenly
spread or the screen noisy, against roughly a sixth where capability is heavy-tailed
and the screen sharp;\footnote{The screened lottery's shortfall to the division
among the top fifth: 0.01 to 0.07 of the optimal gain where capability is evenly spread, or
under heavy tails at the noisiest screen; 0.14 to 0.18 under heavy tails with sharp
to moderately noisy screens.} and the fields where the ranking is cheap to discard
are the fields where review's fine distinctions are worth little. But the comparison
the lottery debate less often poses is spreading: where capability is evenly spread,
the model's answer is not a lottery among a fundable few but thinner grants to many.
The other side of the ledger, what a lottery saves in review costs and proposal
effort, is measured by the contest models cited above; setting the two sides together
prices the choice field by field.

How widely funds should be spread has depended, throughout, on a production function
in which capability and resources are complements. Next we vary the production
technology and ask which of our results survive.

---

## Notes for Aydin (v1)

1. YOUR FOUR DECISIONS applied as stated. The wide-split line is displayed and does
   the decomposition work (note the figure identifies curves by mark, named in the
   caption; in-panel labels are kept only where unambiguous, since the curves
   converge at the noisy end and the family bans leader lines).
2. TERMINOLOGY discipline: "the value of targeting" keeps its \S4 meaning (optimal
   allocation over spreading); the D3 measurement is always qualified ("for a funder
   with records and review," "the informed funder's gain over uniform funding").
   Scheme names introduced in the body before the figure uses them: screen, pool,
   ranked division, equal division, screened lottery, full lottery. I used "ranked
   division"/"equal division" in body prose but "ranked split"/"wide split" as the
   figure's short labels; ONE-CONCEPT-ONE-TERM FLAG: pick one pair and I will align
   body, figure, and appendix (my lean: "ranked division" and "wide division" in
   both, dropping "split").
3. THE DEFENSIBILITY CLAIM is deliberately the comparative one (chance vs the
   screen's ranking, all else fixed): the absolute claim "lotteries are cheap at even
   spread" would be FALSE in the computation (even there, the screened lottery
   captures 0.42-0.54 of the optimal gain vs uniform's 0.83, because concentration
   itself is the mistake). This is why P4 ends with the inversion and P5 adds the
   spreading point. It also bears on the abstract's sentence ("cheap exactly where
   review is worth little"): true of the ranking-vs-chance comparison, not of
   lotteries tout court; flagged for the abstract pass at lock.
4. The Jensen explanation in the body is the plain-language version per your
   instruction; Proposition 2 with the three-line proof and a scope remark is in
   tex/appendix-lottery.tex (renders as Appendix B; prop counter continues from \S4's
   Proposition 1).
5. \S2's ledger frame is restated in place (no backward section reference); the
   Gross and Bergstrom description in the footnote paraphrases \S2's sentence, so the
   two should be checked for verbatim-collision at typeset (currently they differ in
   wording; if you prefer, the \S8 footnote can shrink to the bare citation).
6. GENERATING CONTEXT: the computation is single-round and holds the funded share at
   a fifth (your "state and not optimize"); both stated in the figure caption and
   footnotes. The D3 numbers are two-round with the regime pairing caveat footnoted.
7. Bridge wording assumes \S9 = production technology + concentration ("which of our
   results survive"); adjust when \S9's function analysis lands.
