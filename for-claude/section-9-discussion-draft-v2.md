# Section 9 (Discussion), draft v2 (2026-09-27)

v2 follows Aydin's opening (his first two paragraphs, lightly corrected: "compared",
"and to what it would be in the absence", "determine this value", "bottlenecked",
"the effectiveness") and continues in the same manner: plain, logically ordered,
one idea per sentence, referents named. Order: what we compared and the key
quantity; what information contributes to targeting; what forgoing targeting
loses; what a funder should do; peer review; scope; tests. Title "Discussion".
Claims: section-9-claims.md rows H1-H14 (unchanged grounds). ~1450 words.

---

## 9. Discussion

We have compared the expected impact on research output of targeted funding to what
it would be under equal division of the same budget, and to what it would be in the
absence of any contribution by the funder. Our results isolate one quantity as key to
strategies for optimal grant funding: the value of targeting researchers with the
greatest gap between their capabilities and their resources. Two features of an area
of research determine this value. The first is the funder's budget relative to
researchers' own resources. Where the budget is small, dollars are scarce and where
they go matters; where the budget can relieve most researchers' bottlenecks, the
optimal allocation and equal division produce nearly the same output
(\S\ref{sec:optimal}, Figure~\ref{fig:targeting-value}). The second is how unequally
capability is distributed. Where a few researchers can produce far more than the
rest, getting resources to those few is most of what a funder can add
(\S\ref{sec:review}).

We showed how much various sources of information contribute to the effectiveness of
targeting this gap. Researcher track records are modestly predictive, but they
underdetermine whether a researcher's productivity is bottlenecked by resources or by
capability, and additional rounds of records improve the funder's allocation only
slowly (\S\ref{sec:records}, Figure~\ref{fig:records-rounds}). Peer review supplies
what records cannot: a signal of capability itself. The value of that signal
increases with the two features that set the value of targeting, and with the
signal's informativeness. In a field where capability is heavy-tailed, most of the
signal's value comes from separating the few exceptional researchers from the rest,
and that separation requires little information, so even coarse review recovers most
of the value (\S\ref{sec:review}, Figure~\ref{fig:review-value}). The size of grants
also determines what a funder can learn. On thin grants, a researcher's output is
approximately proportional to their resources and nearly independent of their
capability, so thin grants reveal almost nothing about who is capable; grants
comparable to capability are what make records informative (\S\ref{sec:records},
Figure~\ref{fig:thin-grants}). Finally, a funder that spends later, after observing
output, allocates with more information than one that spends at once. The gain from
spending later is real, but it is smaller than the gain from the review signal at
every horizon and compounding rate we sweep (\S\ref{sec:timing},
Figure~\ref{fig:schedule-vs-signal}). Choosing whom to fund outweighs choosing when.

We also measured what a funder loses by forgoing targeting. Uniform seed grants and
lotteries both give up part of the value of targeting on purpose, so both lose most
output where that value is largest: tight budgets and heavy-tailed capability
(\S\ref{sec:egalitarian}, Figures~\ref{fig:floor-cost} and~\ref{fig:lottery-schemes}).
A lottery over any pool of researchers produces less expected output than an equal
division of the same money over the same pool, because a certain grant produces more
expected output than a chance at a larger grant of the same average size
(Proposition~\ref{prop:lottery}). Where capability is evenly spread, the model's
answer is therefore thin grants to many researchers, not a lottery among a fundable
few. Whether funding should be concentrated on a few researchers or spread across
many depends on the same two features and on the production technology, and the
budget reverses the technology's effect: harder complementarity concentrates grants
when the budget is tight in a heavy-tailed field and spreads them when the budget is
ample (Figure~\ref{fig:concentration}).

These results have a practical consequence for funders. No funding mechanism has a
verdict of its own. Peer review, seed grants, and lotteries each receive opposite
verdicts in different fields, and the two features of the field decide which verdict
applies. A funder choosing among mechanisms should therefore begin by estimating the
two features for its own field. The funder's budget relative to researchers' total
resources is the funder's own arithmetic. The inequality of capability is harder to
observe. The distribution of output is the natural proxy for the distribution of
capability, but where grants are thin, output tracks resources rather than
capability, so the output distribution understates capability inequality
(\S\ref{sec:records}). Two cases show what follows. A small foundation funding a
heavy-tailed field faces the largest value of targeting per dollar. For this
foundation, review directed at identifying the exceptional few is worth most, and
coarse review recovers most of that value; seed grants lose the most output; and a
lottery among the fundable loses about a sixth of what selection by review gains. A
large agency funding an evenly spread field faces a small value of targeting. For
this agency, review of any informativeness adds little; seed grants and lotteries
lose little; and the best use of the budget is thin grants to many researchers.

Three further implications concern peer review. The first concerns where review
effort should go. Since most of review's value comes from separating the exceptional
few from the rest, and fine distinctions among the many applications of comparable
merit add little output, a review system that spends most of its effort ranking the
middle of the distribution spends that effort where the model locates the least
output (\S\ref{sec:review}). The second concerns how to read the evidence on review's
efficacy. The finding that percentile scores barely predict productivity comes from
funded NIH grants with scores in the top fifth \citep{Fang2016}, and the studies that
find modest predictive power also use funded grants \citep{LiAgha2015}. Both kinds of
study measure prediction after review has already made the separation the model
values, inside the band where the model expects review to add least. A weak
association within that band is consistent with review of considerable value over
the whole applicant pool. For the same reason, the mapping in \S\ref{sec:review} from
measured accuracy to signal noise, which treated the measured accuracy as review's
accuracy over the whole pool, is conservative. The third concerns calibration. A
funder that treats review as sharper than it is loses more output than a funder that
treats review as noisier than it is (\S\ref{sec:review}, Figure~\ref{fig:overtrust}).
A funder uncertain about review's accuracy should therefore err toward undertrust.

The model holds several things fixed, and its results should be read within those
limits. There is one funder. Competition among funders, and the process by which past
funding attracts future funding across funders, lie outside the model, though the gap
rule itself directs grants away from researchers whose resources are already large.
Researchers do not respond strategically. Proposal effort, risk choice, and the cost
of review lie outside the model; the contest models cited in \S\ref{sec:literature}
cover that ground. Output is a scalar rate. Projects have no heterogeneity and
exploration has no value, so results about funding unexplored regions of a field have
no counterpart here. The review signal is unbiased noise. Bias and conservatism in
review, which are among the strongest arguments for lotteries, are not modeled. There
is no entry, exit, or career structure. Grants are consumed and inputs are
complements; the timing results depend on complementarity and the information
results do not (Appendix~\ref{app:technology}). The funder is an idealization that
knows the model and updates correctly, so the values of information reported here
are what such a funder extracts, not what a real funder would extract. Finally, the
magnitudes are orderings and comparative statics over the parameters we sweep, and
none is calibrated to a real funder. What a reader should carry to a real field is
the pair of features that organizes the results, not the percentages.

The model makes three predictions that data on funded research can test. First,
capability proxies, whether prior record or review score, should predict output
better the larger a grant is relative to the recipient's other resources, because on
thin grants output is approximately proportional to resources and nearly independent
of capability (\S\ref{sec:records}, Figure~\ref{fig:thin-grants}). A test compares
the predictive power of review scores across grant sizes within one applicant pool.
Second, review's predictive validity, measured over the whole applicant pool rather
than among funded grants, should be higher in fields whose output distributions have
heavier tails (\S\ref{sec:review}). A test compares fields, or one funder's programs,
using the full range of scores where the outcomes of funded and unfunded applicants
can both be observed. Third, where funders run lotteries beside review, the output
lost to the lottery should be larger in heavier-tailed fields
(\S\ref{sec:egalitarian}, Figure~\ref{fig:lottery-schemes}). The lottery programs now
running at several funders supply the comparison. Each prediction is comparative, and
each can fail.

---

## Notes for Aydin (v2)

1. Your opening is kept nearly verbatim; corrections: "We have compared", "and to what
   it would be in the absence of any contribution by the funder", "determine this
   value" (for "this field"), "bottlenecked by resources or by capability", "the
   effectiveness of targeting". Your "exhibit rapidly diminishing returns in
   informativeness" is rendered as "additional rounds of records improve the funder's
   allocation only slowly," which is the claim \S5 verified (correlation with the
   optimal grants 0.08 in round one, 0.41 by round twenty); "rapidly diminishing" is
   not what the trajectory shows, so I did not keep the word.
2. New in v2 relative to v1: the second paragraph now walks the information sources in
   one order (records, review, grant depth, timing); the third paragraph collects
   what forgoing targeting loses, including Proposition 2 and the concentration
   result, so the discussion covers \S8 in full. Cut: the "stakes times capture"
   framing and the closing sentence "the results move together because they measure
   one quantity from different sides."
3. Still pending from v1: CONFIRM Li & Agha 2015 uses funded grants only (Fang
   confirmed from the paper's text); the one-clause \S6 consistency edit (wording in
   v1 note 2).
4. "a lottery among the fundable loses about a sixth of what selection by review
   gains" is \S8's single-round heavy-tailed, tight-budget cell (0.14 to 0.18 of the
   optimal allocation's gain at sharp to moderate review).
5. \S10 conclusion: sketch in section-9-prep.md \S5; draft after your verdict.
