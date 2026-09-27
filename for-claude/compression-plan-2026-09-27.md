# Compression plan (2026-09-27): from 37 pages toward 20

No changes made. This is the plan for your review.

## 1. Where the 37 pages go

| Part | Pages now | Words | Notes |
|---|---|---|---|
| S1 Introduction | 1-2 | 877 | |
| S2 Related literature | 3-4 | 626 | |
| S3 Model | 5-7 | 1262 + 81 fn | longest setup |
| S4 Optimal allocation | 8-10 | 917 + 105 fn | two figures |
| S5 Track record | 11 | 482 + 41 fn | one figure |
| S6 Peer review | 12-16 | 1124 + 407 fn | three figures; heaviest footnotes |
| S7 Timing | 17-19 | 1074 + 336 fn | three figures |
| S8 Spreading | 20-23 | 1635 + 407 fn | three figures; longest section |
| S9 Discussion | 24-27 | 1517 | no figures |
| Appendix A | 28-32 | 2179 | full proof, 2 corollaries, 2 examples, 3 remarks |
| Appendix B | 33 | 397 | |
| Appendix C | 33 | 651 | one figure, one table |
| References | 34-37 | | plainnat here; the journal's style will differ |
| Captions (13) | | 1311 | 71 to 166 words each |
| Footnotes (all) | | 1377 | |

Total about 15,400 words including appendices, captions, and footnotes. Body prose alone is 9,500.

Two facts shape the plan. First, layout is doing a lot of the damage: the file is set at one-and-a-half spacing with 4 cm side margins. Single spacing alone gives 33 pages; single spacing with 3 cm margins gives 30 (I compiled both). Second, even at 30 pages the paper is fat in specific places: the captions, the footnotes in S6 to S8, the restatement paragraphs in S9, Appendix A, and the model section. The plan below cuts about a third of the words. With the layout change, that lands at roughly 21 pages including a four-page reference list, or about 17 pages without references. Under 20 with references included needs both the content cuts and the layout change; I recommend doing both, since the journal will reset the layout anyway and the content cuts are improvements on their own.

## 2. Principles for the cuts

1. Every result keeps its figure or proof (your rule). Nothing figure-backed loses its figure; figures get smaller and paired, not removed.
2. Generating context (parameters, number of simulated populations, allocator notes) moves out of captions and footnotes into the simulation-specification appendix that the paper already promises twice ("Appendix [X]"). That appendix is owed anyway; writing it is what lets the captions and footnotes shrink without losing rigor. It carries one parameter table and one table of run settings per figure.
3. Each result is stated once in the body, once in S9 as a one-clause reference, and nowhere else. The abstract, the intro's results paragraphs, and S9's first three paragraphs currently say the same things three times.
4. Mechanism digressions that no later result depends on go: the two-sources-of-growth paragraph in S7, the sharp-end reversal footnote in S6, the fitted exponents.
5. Bridges: one short sentence at the end of a section, none at the start.

## 3. Section-by-section

### S1 Introduction: 877 -> about 600

- Keep P1 (the stakes and the three disputes) and P2 ("Consider the challenge of the grant funder"), trimmed by a sentence each.
- P3 (the model in words): keep, it is the reader's first picture.
- P4 and P5 ("Our results are as follows", "Peer review supplies..."): these two paragraphs, about 330 words, are a full preview of S4 to S8 and duplicate the abstract and S9. Cut to one paragraph of about 150 words: the gap rule in one sentence; records underdetermine the gap; review's value is set by capability inequality, informativeness, and budget; seed grants and lotteries lose exactly where review is worth much; whom over when.
- P6 (two assumptions): keep as you wrote it.
- Roadmap: keep, one line per section.

### S2 Related literature: 626 -> about 500

- Concentration paragraph: cut the Heesen sentence (questionable research practices over a career) unless you want the citation; it does not bear on the model.
- Lottery paragraph: the three implementing funders become one clause ("several funders now run partial lotteries [cites]").
- Formal-models paragraph: keep.
- Fold the final one-sentence paragraph ("Our formulation and analysis...") into the end of the formal-models paragraph.

### S3 Model: 1262 -> about 850

- Production function: keep the definition and the harmonic-mean gloss. The Cobb-Douglas and Leontief sentences shrink to one: "Appendix C repeats the analysis for a family of production functions from Cobb-Douglas to Leontief."
- Accumulation paragraph: keep the update rule and the consumed-grants assumption (one sentence, since S1 now states it).
- Budget paragraph: keep the budget scale and the re-planning sentence; the normalization sentence points to the simulation appendix.
- Signals paragraph: keep.
- Strategies paragraph and Table 1: keep both; the table does the work, so the paragraph drops to the two dimensions and the three baselines.
- Parameter-range paragraph: keep the ranges (they are the evidence base); the Lotka justification becomes one clause with the citation.
- "Randomness enters our model" paragraph (about 200 words on expected-output scoring and common random numbers): move whole to the simulation appendix. In the body, one sentence: "We score each strategy by expected output and compare strategies on identical simulated populations (Appendix D)."
- The closing sentence stays.

### S4 Optimal allocation: 917 -> about 750

- Setup and marginal-value display: keep.
- Proposition: keep.
- "We can understand the rule as follows": this paragraph restates the proposition twice (targets and gaps; small budget, large budget). Cut to about half: the target-and-gap gloss, then the small-versus-large budget sentence, then the frontier sentence that introduces Figure 1. The water-filling and poverty-gap footnote stays, short.
- Intuitive-schemes paragraph: keep; it is the payoff of the section. The two "makes this exact" footnotes become one: "Examples 1 and 2 in Appendix A make both failures exact."
- Budget paragraph (value of targeting vanishes): keep.

### S5 Track record: 482 -> about 400, and consider merging into S6

- The "what output can identify in principle" paragraph: keep; it is the argument.
- The measurement paragraph: keep the shortfall numbers in the footnote only.
- The thin-grants paragraph: keep, but cut its last two sentences ("In our simulations of a resource-poor community... S7 measures this" and "For new and resource-poor fields the implication is direct"), since S7 states the result with its figure.
- Structural option: make S5 the opening subsection of S6 under a single heading, "What the funder can learn: records and peer review." Saves a section boundary, its two bridges, and a heading; the two sections are one argument (records cannot separate capability from resources; review can). I recommend it.

### S6 Peer review: 1124 + 407 fn -> about 800 + 150 fn

- Opening: one sentence (the value of review, defined in output).
- "Direction of review's effect" paragraph: keep the sentence and the correlation numbers in a short footnote. Cut the sharp-end reversal (the 2 percent non-monotonicity and its mechanism, about 120 words of footnote): it is a curiosity that nothing else uses. If you want to keep a trace, one clause: "the increase is not strictly monotone at the sharpest signals."
- "Two further effects" and the records-do-not-catch-up paragraph: keep, tightened; Figure 5 (records-rounds) stays.
- Main result paragraph: keep; it is the paper's central result. Footnote keeps the three retention percentages and the half-value noise levels; drop the sentence restating the sweep range.
- "A consequence concerns where review effort should go": keep, two sentences.
- "In the other fields" paragraph: merge into the main result paragraph as its last two sentences.
- Efficacy paragraph: this and S9's efficacy paragraph overlap. Keep here only the mapping (AUC 0.54 corresponds to the noisiest signal in the sweep; at that noise review still recovers a third of the shortfall under heavy tails and nothing under even spread) in about 80 words plus one footnote; S9 keeps the interpretation (funded-grants truncation, no single number). Delete "Accuracy is one input to review's value; ... No single number settles the dispute" from S6, since S9 says it.
- Overtrust paragraph: keep, footnote to one sentence.

### S7 Timing: 1074 + 336 fn -> about 750 + 120 fn

- Opening paragraph: cut to two sentences.
- Spend-early argument: keep the argument and its failure (this is good); cut the planner-free footnote to its two numbers.
- Simulations paragraph and Figure 8: keep; footnote keeps only "every cell whose schedule leans early has epsilon <= 0.03".
- Grant-size paragraph and Figure 7: keep, shorter; the footnote keeps the three percentages and drops the budget-scaling sentence (moves to the simulation appendix).
- Two-sources-of-growth paragraph and its footnote (about 200 words): cut entirely. It explains why the schedule leans late; the results do not depend on the reader having the mechanism, and Figure 8's caption already says compounding sets the lean. If you want to keep one line: "The lean toward late spending comes from growth in output the funder did not pay for; growth from grant-added output alone would lean the schedule slightly early."
- "What is the schedule choice worth" paragraph and Figure 9: keep. Cut the fitted-exponents footnote. Keep the Cobb-Douglas sentence with its pointer to Appendix C.

### S8 Spreading funds: 1635 + 407 fn -> about 1100 + 180 fn

- Opening: keep the two alternatives and their rationales; cut the Gross-Bergstrom footnote, which repeats S2 word for word (point to S2 instead).
- "Both alternatives lose the same thing" paragraph: keep; footnote keeps the four ratios only.
- Seed grants paragraph and Figure 10: keep; drop the factor-law footnote (0.10 to 0.86) or cut it to one clause.
- Lottery paragraph: the five-scheme description (about 130 words) becomes a compact enumeration in one sentence each, and the results sentences stay. Footnote keeps the 0.75 / 0.61 / 0.56 decomposition, one line.
- "When does a lottery among the fundable lose little" paragraph: merge its first half into the lottery paragraph; keep the closing point (the debate rarely asks whether to concentrate at all; where capability is evenly spread the answer is thin grants to many) as the bridge to the concentration paragraph.
- Concentration paragraph: cut from about 330 to about 180 words. Keep: the question, the family in one sentence with the Appendix C pointer, the three findings (tight budget concentrates; ample budget spreads; even-and-ample is non-monotone), the closing sentence. The generating-context footnote moves to the simulation appendix; the Gini footnote keeps only the endpoints (0.92 to 0.95; 0.35 to 0.21; 0.16, peak 0.196, 0.185).

### S9 Discussion: 1517 -> about 900, and absorb the conclusion

- P1 (key quantity, two features): keep at about 120 words.
- P2 and P3 (information sources; what forgoing targeting loses): these restate S5 to S8 at about 550 words. Cut to one paragraph of about 200 words with one clause per result and its figure reference. The reader has just read the sections.
- P4 (no mechanism is good or bad in itself; estimate the two features; two cases): keep at about 200 words; it is the practical content.
- P5 (peer review: effort, efficacy evidence, calibration): keep at about 180 words. Since S6's interpretive sentences move here, nothing is lost.
- P6 (limitations): keep at about 170 words; the list is already compact.
- P7 (three tests): keep at about 130 words, one sentence of prediction and one of test per item.
- S10 Conclusion: do not write a separate section. Retitle S9 "Discussion and conclusion" and end it with a four-sentence closing paragraph after the tests (the allocation problem the intro posed; fund the gap; the two features decide what review, seed grants, and lotteries are worth; whom over when). A separate conclusion would restate S9, which restates the body. Update the intro roadmap accordingly.

## 4. Captions: 1311 -> about 400

Rule: a caption has three parts and nothing else. (1) One sentence saying what the figure shows, in the form of the result. (2) What is on the axes, if the axis labels do not already say it. (3) Marks, only where curves are not labeled in the plot. Parameters, seed counts, and generating context go to the simulation appendix's per-figure table, which the caption cites once as "(settings: Appendix D)". Numbers that appear in the body's footnotes do not appear again in the caption.

Two examples.

Figure 4 now (125 words): "The value of review by field. Each curve shows the share of the records-only shortfall that review recovers: the additional expected output of a funder allocating with the review signal over the same funder allocating without it, as a fraction of the output the records-only funder leaves unrealized relative to complete information. Signal noise tau_K increases to the right, so informativeness decreases; the default is tau_K = 1. Under heavy-tailed capability (alpha_K = 1.3), review recovers most of the shortfall across a wide range of noise; where capability is spread more evenly (alpha_K = 3.5), review recovers little at any noise. Fifty simulated populations per point; the records-only shortfall is 24.6, 11.4, and 3.5 percent of no-funding output at alpha_K = 1.3, 2, and 3.5."

Figure 4 proposed (38 words): "Review is worth most where capability is heavy-tailed. The share of the records-only shortfall that the review signal recovers, as review noise increases (informativeness decreases), for three capability distributions. Settings: Appendix D."

Figure 11 now (166 words) -> proposed (about 55): "What five funding schemes gain, by field. Each scheme's gain over no funding as a share of the optimal allocation's gain, as review noise increases, in a heavy-tailed field with a tight budget (left) and an evenly spread field with a large budget (right). Filled circles: equal division among the top fifth by review score; open circles: top two fifths; triangles: the partial lottery; dashed: uniform funding; dotted: a lottery over all researchers. Settings: Appendix D."

Apply the same rule to all thirteen. Figures 1 to 3 (analytic) need no settings pointer.

## 5. Figures: pair them

Six figures become three two-panel figures, each at the current width and about 5 cm tall: Figures 2 and 3 (value of targeting; thin grants), Figures 4 and 5 (review value; records do not catch up), Figures 7 and 8 (grant size and schedule; compounding and schedule). Figures 10 and 11 already sit together; Figure 11 is two-panel. This saves about two pages and puts each pair of results side by side where the text discusses them together. Figure 12 stays in Appendix C; Figure 13 stays in S8.

## 6. Footnotes: 1377 -> about 550

Rule: a footnote carries at most the numbers that back the sentence above it, in one or two lines. It carries no mechanism, no curiosity, and no generating context (which moves to Appendix D). Delete outright: the sharp-end reversal (S6), the fitted exponents (S7), the Gross-Bergstrom repeat (S8), the factor-law footnote (S8, or one clause), the concentration generating-context footnote (S8). Shorten to one or two lines: all others.

## 7. Appendices

- Appendix A (gap rule, 2179 words, five pages): essential, since Proposition 1 is the paper's analytic result. Compress to about 1300 words: keep Lemmas 1 to 3 with proofs (Lemma 2's proof can be four lines); keep the Proposition and its four-step proof; keep Corollary 2 (value of targeting vanishes) since S4 cites it; keep Corollary 1 as a statement without proof (it is immediate from the formula); keep Examples 1 and 2 (S4 cites them) with their numbers only, no commentary sentence; delete Remarks 1 to 3 (water-filling and poverty-gap are already in the S4 footnote; the scope remark is one sentence at the top of the appendix).
- Appendix B (lotteries, 397 words): essential, since Proposition 2 backs the S8 lottery claim. Keep as is; delete the closing remark's last sentence.
- Appendix C (production technology, 651 words plus Figure 12 and Table 2): essential for the robustness claims S6 and S7 point to. Keep; cut the family-definition paragraph by a third (the S8 sentence already introduces the family).
- Appendix D (simulation specification): new, about 400 words plus two tables. Budget normalization, the default parameter table, the per-figure settings table (populations per point, rounds, allocator version, any caveat such as the Leontief kink rule), scoring by expected output, common random draws, and the code pointer. This appendix is what makes the caption and footnote cuts rigorous rather than lossy.

## 8. Layout

Single spacing and 3 cm side margins for the working version. The journal will typeset in its own format, so the current one-and-a-half spacing and 4 cm margins buy nothing but pages. If you want a double-spaced version for referees, that is a compile option, not a reason to keep it as the working layout.

## 9. Expected outcome

| | Now | After |
|---|---|---|
| Body words | 9,500 | about 6,700 |
| Footnote words | 1,377 | about 550 |
| Caption words | 1,311 | about 400 |
| Appendix words | 3,227 | about 2,500 + Appendix D 400 |
| Figures in body | 12 (+1 appendix) | 9 blocks (three paired) + 1 appendix |
| Pages, current layout | 37 | about 27 |
| Pages, single spacing, 3 cm margins | 30 | about 21 with references, about 17 without |

## 10. Order of work, if approved

1. Layout change and figure pairing (mechanical; measures the true starting point).
2. Appendix D written from the existing generating-context notes (so nothing is lost before it is cut elsewhere).
3. Captions and footnotes cut to the rules above.
4. Body cuts, section by section, each delivered for your approval with a word count before and after.
5. Appendix A compression.
6. S9 closing paragraph and roadmap update.

Decisions for you: merge S5 into S6 (recommended); absorb the conclusion into S9 (recommended); cut the two-sources-of-growth paragraph entirely or keep one line; cut the sharp-end reversal entirely or keep one clause; the Heesen sentence in S2.
