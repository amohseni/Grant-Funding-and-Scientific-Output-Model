# Section 9 (Discussion), draft v1 (2026-09-09)

Per decisions: title "Discussion"; prose only (no table); current-support remark cut;
three tests. Function analysis: section-9-prep.md. Claims: section-9-claims.md
(rows H1-H14). ~1350 words. (section-9-draft-v1.md is the SUPERSEDED
production-technology draft; this file is the discussion.)

---

## 9. Discussion

Our results concern one quantity, the value of targeting: the expected output that
the optimal allocation adds over an equal division of the same budget. Two features
of a field set this value. The first is the funder's budget relative to researchers'
own resources. Where the budget is small, dollars are scarce and where they go
matters; where the budget can relieve most researchers' bottlenecks, the optimal
allocation and equal division produce nearly the same output
(\S\ref{sec:optimal}, Figure~\ref{fig:targeting-value}). The second is how unequally
capability is distributed. Where a few researchers can produce far more than the
rest, getting resources to those few is most of what a funder can add
(\S\ref{sec:review}). These two features set what is at stake. The information
results say how much of the stake a funder obtains. Records alone obtain little of it,
and further rounds of records add little (\S\ref{sec:records},
Figure~\ref{fig:records-rounds}). Review obtains a share that increases with the same
two features and with review's informativeness (\S\ref{sec:review},
Figure~\ref{fig:review-value}). A funder that spends later, after observing output,
obtains more of it, but the addition is smaller than what the review signal adds at
every horizon and compounding rate we sweep (\S\ref{sec:timing},
Figure~\ref{fig:schedule-vs-signal}); in a resource-poor field, a funder whose grants
are deep obtains far more of it than one whose grants are thin
(\S\ref{sec:timing}, Figure~\ref{fig:depth-schedule}). Seed grants and lotteries give
up part of the stake on purpose, and lose most where the stake is largest
(\S\ref{sec:egalitarian}, Figures~\ref{fig:floor-cost} and~\ref{fig:lottery-schemes}).
Complementarity between capability and resources raises the stake and with it every
quantity that depends on it (Appendix~\ref{app:technology}). In every field we
examine, the results move together because they measure one quantity from different
sides.

For a funder, the first consequence is that no mechanism carries a verdict of its
own. Peer review, seed grants, and lotteries each receive opposite verdicts in
different fields, and the two features decide which. A funder choosing among them
should begin by estimating the two features for its own field. The first is the
funder's own arithmetic: its budget against the total resources of the researchers it
might fund. The second is harder to observe. The distribution of output is the natural
proxy for the distribution of capability, with the caveat of \S\ref{sec:records}: on
thin grants output tracks resources, not capability, so a field's output distribution
understates its capability inequality where grants are thin. Two cases show what
follows. A small foundation funding a heavy-tailed field faces the largest stake per
dollar. In this field, review directed at identifying the exceptional few is worth
most, and coarse review recovers most of that value; seed grants lose the most output
here; a lottery among the fundable loses about a sixth of what selection by review
gains. A large agency funding an evenly spread field faces a small stake. In this
field, review of any informativeness adds little; seed grants and lotteries lose
little; and the best use of the budget is thin grants to many researchers, not a
lottery among a few.

Three implications concern peer review. The first is where review effort should go.
Most of review's value comes from separating the few researchers of unusually high
capability from the rest, and that separation requires little information; fine
distinctions among the many applications of comparable merit add little output
(\S\ref{sec:review}). A review system that spends most of its effort ranking the
middle of the distribution spends it where the model locates the least output. The
second is how to read the evidence on review's efficacy. The finding that percentile
scores barely predict productivity comes from funded NIH grants with scores in the top
fifth \citep{Fang2016}, and the studies that find modest predictive power also use
funded grants \citep{LiAgha2015}. These studies measure prediction after review has
made the separation the model values, inside the band where the model expects review
to add least. A weak association there is consistent with review of considerable value
over the whole applicant pool, and the mapping in \S\ref{sec:review} from measured
accuracy to signal noise, which treated the measured accuracy as review's accuracy
over the whole pool, is for that reason conservative. The third is calibration. A
funder that treats review as sharper than it is loses more output than a funder that
treats review as noisier than it is (\S\ref{sec:review}, Figure~\ref{fig:overtrust}).
A funder uncertain about review's accuracy should err toward undertrust.

Two implications concern records and time. Cross-sectional output reflects resources
as well as capability (\S\ref{sec:records}), so allocation by publication record,
whether by a committee or by a bibliometric formula, rewards the well-resourced along
with the capable, and can produce less output than an equal division
(\S\ref{sec:optimal}). What the funder needs is the difference between a researcher's
capability and their resources, and records alone do not supply it. For resource-poor
fields the lever is depth, not haste: grants comparable to capability are what make
records informative, and a funder that spends everything in the first round allocates
before observing any output (\S\ref{sec:records}, \S\ref{sec:timing}). More generally,
how the budget is spread over time matters less than to whom it goes, at every horizon
and compounding rate we sweep (\S\ref{sec:timing}, Figure~\ref{fig:schedule-vs-signal});
the questions of schedule that funders debate are second-order relative to selection.

The model holds several things fixed, and its results should be read within them.
There is one funder; competition among funders, and the process by which past funding
attracts future funding across them, are outside the model, though the gap rule itself
directs grants away from researchers whose resources are already large. Researchers
do not respond strategically: proposal effort, risk choice, and the cost of review lie
outside the model, and the contest models cited in \S\ref{sec:literature} hold that
ground. Output is a scalar rate: projects have no heterogeneity and exploration has no
value, so results about funding unexplored regions of a field have no counterpart
here. The review signal is unbiased noise; bias and conservatism in review, which are
among the strongest arguments for lotteries, are not modeled. There is no entry, exit,
or career structure. Grants are consumed and inputs are complements; the timing
results depend on complementarity and the information results do not
(Appendix~\ref{app:technology}). The funder is an idealization: it knows the model and
updates correctly, and the values of information reported here are what such a funder
extracts, not what a real funder would. Finally, the magnitudes are orderings and
comparative statics over the parameters we sweep; none is calibrated to a funder. What
a reader should carry to a real field is the pair of features that organizes the
results, not the percentages.

The model makes predictions that data on funded research can test. First, capability
proxies, whether prior record or review score, should predict output better the larger
a grant is relative to the recipient's other resources, because on thin grants output
is approximately proportional to resources and nearly independent of capability
(\S\ref{sec:records}, Figure~\ref{fig:thin-grants}). A test compares the predictive
power of review scores across grant sizes within one applicant pool. Second, review's
predictive validity, measured over the whole applicant pool rather than among funded
grants, should be higher in fields whose output distributions have heavier tails
(\S\ref{sec:review}). A test compares fields, or one funder's programs, using the full
range of scores where the outcomes of funded and unfunded applicants can both be
observed. Third, where funders run lotteries beside review, the output lost to the
lottery should be larger in heavier-tailed fields (\S\ref{sec:egalitarian},
Figure~\ref{fig:lottery-schemes}). The lottery programs now running at several funders
supply the comparison. Each prediction is comparative, and each can fail.

---

## Notes for Aydin (v1)

1. CHECK PASSED (Fang): Fang, Bowen & Casadevall 2016 analyze 102,740 FUNDED R01
   grants with percentile scores of 20 or better; AUC 0.54 (their text, fetched
   2026-09-09). CHECK PARTIAL (Li & Agha): the paper is paywalled here; my recall is
   137,215 funded R01 grants, 1980-2008, with a one-standard-deviation better score
   associated with more publications and citations. The 2022 NBER review of peer
   review for funding confirms the general point ("information typically exists only
   on funded projects"). Confirm from the paper before lock.
2. CONSISTENCY EDIT PROPOSED in \S6 (not made): after "In our model, a signal of that
   accuracy is very noisy review.", add: "The accuracy is measured among funded grants
   with the best scores, so review's accuracy over the whole applicant pool is likely
   higher, and the mapping is conservative (\S\ref{sec:discussion})." Say the word.
3. MODAL FORCE, P1: "In every field we examine, the results move together because
   they measure one quantity from different sides" is a structural description backed
   by the co-movements (Fig 2, 5, 10, 11; E4's factor law; Appendix C), not a theorem.
   If you want it weaker: "In every field we examine, the results move together."
4. P2's "about a sixth" is \S8's single-round, heavy-tailed, tight-budget figure
   (0.14-0.18 of the optimal allocation's gain at sharp-to-moderate review). The
   sentence pairs it with a small foundation in a heavy-tailed field, which matches
   the cell. "Coarse review recovers most of that value" is \S6's 82 percent at three
   times the default noise under alpha_K = 1.3.
5. P5's "the gap rule itself directs grants away from researchers whose resources are
   already large" is Proposition 1 (grants decrease in R_i0). The Matthew-effect
   literature (Bol et al. 2018 is in the bibliography) is not cited here; add if you
   want the contrast explicit.
6. The idealized-funder sentence deliberately does not claim the model's values are
   upper bounds for real funders: the forward funder is a heuristic optimizer, so the
   bound is not proven.
7. Cut per your decision: the current-support (resource declaration) remark.
8. \S10 conclusion: not drafted; sketch in section-9-prep.md \S5. Next after your
   verdict on this.
