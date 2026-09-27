# Section 9 (Discussion: implications and limitations), prep (2026-09-09)

Status: function analysis and architecture for Aydin's verdict. No prose drafted.
Promises to pay: the intro's "\S9 discusses implications and limitations"; \S2's claim
that the model "allows us to address questions raised by these literatures"; the
title's "optimal science funding." Companion: \S10 conclusion (sketched at the end).

## 1. What the section is for

Three jobs, in order of importance.

(a) Unify. The paper has eight sections of results and no statement of what holds
them together. The unifying claim is available and is true in the model: every result
is a statement about one quantity, the value of targeting (the expected output the
optimal allocation adds over equal division), factored into two parts.

  - What is at stake: the value of targeting itself. Two features of a field set it:
    the funder's budget relative to researchers' own resources, and how unequal
    capability is. Tight budgets and heavy tails raise it; ample budgets and even
    spread make it vanish (\S4, Corollary 2; Fig 2).
  - What can be captured: the share of that value a funder actually obtains, which
    is set by its information. Records alone capture little and improve slowly (\S5);
    review captures a share that rises with the same two features and with review's
    informativeness (\S6); depth and scheduling raise capture modestly by making
    records informative (\S5, \S7); seed grants and lotteries lower capture on purpose
    (\S8). Complementarity scales the stakes and with them everything downstream
    (Appendix C).

  Stakes times capture. Every figure in the paper is a point on one of these two
  factors. This is the sentence the paper has earned and has not yet said.

(b) Tell a funder what to do with it. The model's practical content is not a
recommendation about any mechanism; it is that the same mechanism gets opposite
verdicts in different fields, and that two features of the field decide which. So the
practical instruction is diagnostic first: estimate the two features, then read off
what follows for review, seed grants, lotteries, and schedule. Two worked cases carry
this (the skill's example-density rule): a small foundation in a heavy-tailed field;
a large agency in an evenly spread one.

(c) State scope without apology. What is held fixed, which results are scoped
(timing to complementarity), what magnitudes mean (orderings and comparative statics,
not transplantable numbers), and what would test the model.

## 2. Architecture (six movements, ~1400 words)

P1 Unification: stakes times capture, as above, with each clause pointing to its
figure. Candidate device: Table 3, "What two features of a field decide." Rows are
the paper's questions (Is targeting worth much? Is review worth investing in? Where
should review effort go? Do seed grants lose much? Is a lottery defensible? Does
schedule matter? Concentrate or spread?). Two columns: tight budget and heavy-tailed
capability; ample budget and evenly spread capability. Cells are direction words with
a figure pointer, no numbers. Every cell is already figure-backed. Risk: qualitative
tables can read as glib; mitigation is that each cell cites its figure, and the table
replaces a paragraph of summary that would otherwise be walking through the sections.
YOUR CALL: table or prose only.

P2 The diagnostic: what a funder should estimate before choosing a mechanism, and
why the two features are estimable in principle (budget relative to the field's
total resources is the funder's own arithmetic; the tail of the capability
distribution is proxied, imperfectly, by the tail of the output distribution,
with the \S5 caveat that thin grants hide capability). Two worked cases.

P3 Peer review. Three implications, each earned by \S6. (i) Effort should go to
separating the exceptional from the rest, not to ranking the middle; review systems
that spend most of their effort on fine distinctions among fundable proposals are
spending it where the model says the least output is. (ii) How to read the efficacy
evidence: studies that find scores barely predictive measure prediction among FUNDED
grants, that is, after the separation the model values has been made and inside the
band where the model expects review to add least; the null is consistent with review
being worth a great deal across the whole applicant pool. [CHECK before drafting:
that Fang, Bowen & Casadevall 2016 and Li & Agha 2015 both analyze funded grants
only. I believe both do; verify from the papers.] NOTE the tension with \S6's current
text, which maps the AUC of 0.54 to "very noisy review" over the whole pool. If the
measurement is truncated, \S6's mapping is conservative (true review is sharper than
tau = 20), and \S6 should say so in one clause. (iii) Miscalibration is asymmetric:
overtrusting review loses more than undertrusting it, so a funder uncertain about
review's accuracy should lean toward undertrust (Fig 6).

P4 Records, depth, and timing. Cross-sectional records reward resources as much as
capability, so allocation by bibliometric record inherits the track-record failure of
\S4; the remedy in the model is information about resources, which funders already
collect in part (current-support declarations) and which the model gives a reason to
weight. [Modal force: "gives a reason to weight," not "shows should be weighted";
the model's resource signal is an idealization of such declarations.] For poor
fields, the lever is depth, not haste: grants comparable to capability are what make
records informative, and spending early does not (\S5, \S7). Choosing whom outweighs
choosing when, so the funding-schedule debates are second-order relative to selection.

P5 Scope. In prose, as precision. Held fixed: a single funder (no competition among
funders, so no Matthew-effect dynamics across funders, though the gap rule itself
cuts against funding by accumulated resources); researchers who do not respond
strategically (proposal effort, risk choice, and review costs are outside; Gross and
Bergstrom's contest models hold that ground); a scalar output rate (no project
heterogeneity, no exploration value, so Avin's landscape results have no counterpart
here); an unbiased review signal (bias and conservatism, the lottery advocates'
strongest points, are not modeled); no entry, exit, or career structure; the two
substantive assumptions of \S1 (grants consumed, inputs complementary), with the
timing results scoped to complementarity (Appendix C) and the information results
not. The funder is an idealized Bayesian who knows the model: the value it extracts
from information is a ceiling for real funders with the same information, and the
shortfalls it suffers are floors. Magnitudes are orderings and comparative statics
over swept parameters, not calibrated to any funder; the two features that organize
the results are the ones a reader should carry to a real field.

P6 Tests. Two or three concrete predictions the model makes and the data that would
check them. (i) Capability proxies (prior record, review score) should predict
output better the larger the grant is relative to the recipient's other resources,
because output on thin grants is approximately proportional to resources and nearly
independent of capability (\S5, Fig 3). (ii) Review's predictive validity, measured
over the whole applicant pool, should be higher in fields with heavier-tailed output
distributions (\S6, and the accuracy-to-noise mapping footnote). (iii) The output
lost to lotteries, where funders run them beside review, should be larger in
heavy-tailed fields (\S8, Fig 11). Each is stated as what to measure and against
what, not as "future work would be valuable."

Bridge to \S10: none needed beyond the section ending on the tests.

## 3. What is worth the reader's attention (ranked)

1. Stakes times capture (P1). Without it the paper is a list of results.
2. Diagnose the field before choosing the mechanism (P2). This is what a funder
   can act on tomorrow.
3. Where review effort should go, and the reading of the efficacy evidence (P3).
   This is the paper's sharpest contact with an existing dispute.
4. Whom over when (P4). Already stated in \S7; here it becomes advice.
5. The lottery verdict as field-conditional (P2/P3). Already in \S8.

## 4. Honesty checks before drafting

- Stakes times capture is a description of the model's structure, not a theorem.
  Draft it as "in every field we examine," with the co-movements that ground it
  (E4's factor law; \S6 and \S8 turning on the same two features; Appendix C).
- Records-only output can fall below equal division (\S4 examples), so "capture" is
  not a share bounded in [0,1]. Say "share" loosely or say "the part of the value a
  funder obtains."
- P3(ii) needs the CHECK above and a consistency edit in \S6.
- P4's current-support remark is an interpretation; keep the modal force low.
- No new results in \S9. Everything cites a body figure, an appendix, or a
  cited paper.
- Voice: same as the body. No "we hope," no "future work," no announced structure.

## 5. \S10 conclusion (sketch)

Four to six sentences. Situate: the allocation problem the intro posed. Restate:
fund the gap; records underdetermine it; review's value and the alternatives' losses
are set by budget scale and capability inequality; whom over when. Do-not-show: no
claim about any real funder's parameters. Final sentence, bold and bounded: the
model's answer to "should we use review, lotteries, or seed grants" is that the
question is malformed until the field's two features are named. Candidate wording to
come with the \S9 draft.

## 6. Decisions for Aydin

1. Table 3 (qualitative, figure-pointered) or prose only.
2. Include P3(ii), the funded-grants reading of the efficacy evidence, pending the
   CHECK, and make the one-clause consistency edit in \S6.
3. P4's current-support remark: keep (low modal force) or cut.
4. P6: two tests or three; which.
5. Section title: "Discussion" or "Implications and limitations."
