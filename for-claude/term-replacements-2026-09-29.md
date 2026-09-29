# Coined-term sweep (2026-09-29)

Scope: main.tex, the four appendices, all figure captions, pnas-main.tex, the four PNAS
composite figures. Rule applied: a word used in a sense a reader must guess at, or a label
we coined for a quantity that already has a plain name, is replaced by the plain name.
Counts are occurrences across all files. Every replacement is reversible by search.

## Replaced

| Term | Replacement | Count | Reason |
|---|---|---|---|
| coarse review / coarse signal | even noisy review / even a noisy signal from peer review | 8 | "coarse" left the object (review noise) unnamed |
| sharp review; review is sharp; sharpest; sharp enough signal | reliable review; review is reliable; most reliable; reliable enough signal | 22 | "sharp" for "low noise" is our coinage |
| Sharpening review beyond a modest level | Making review more reliable beyond a modest level | 1 | same |
| treating a noisy signal as sharp / treats review as sharper than it is | treating a noisy signal as reliable / treats review as more reliable than it is | 4 | same |
| thin grants; a thin grant; grants are thin | small grants; a small grant; grants are small | 15 | "thin" is figurative; the text already says "small relative to capability" where the size matters |
| with thin resources | with few resources | 1 | same |
| spread the budget thinly across many researchers | spread the budget in small grants across many researchers | 1 | same |
| resource-poor fields / community; a poor field | fields whose researchers have few resources / a community with few resources; such a field | 9 | "poor" is figurative and its object (baseline resources) was unnamed |
| Resource poverty raises the early share | Scarce baseline resources raise the early share | 2 | same |
| the compounding rate, not resource poverty, sets how late the schedule leans (captions) | the compounding rate, not the scarcity of baseline resources, sets how far spending shifts toward later rounds | 2 | "poverty" and "leans" both replaced |
| holds a below-even share ... deploys the bulk | spends less than an even share ... spends the bulk | 2 | "below-even share" was a coinage; "deploys" replaced by "spends" |
| puts an above-even share in the early rounds | spends more than an even share in the early rounds | 1 | same |
| whose schedule leans early | whose optimal schedule spends more than an even share early | 1 | "leans" |
| the lean toward late spending | the shift of spending toward later rounds | 3 | same |
| would lean the schedule slightly early | would shift the schedule slightly toward earlier rounds | 1 | same |
| larger when spending leans late (Table 2 caption) | larger when spending shifts toward later rounds | 1 | same |
| the informed funder's grants | the review-informed funder's grants | 3 | two labels for one funder; captions already used "review-informed" |
| the exceptional few; the exceptional researchers | the few researchers of unusually high capability | 5 | "exceptional" was an evaluative label for a defined group |
| fundable applications / a lottery among the fundable / a fundable few | applications that pass a review threshold / a lottery among applicants above a review threshold / a few applicants above a review threshold | 7 | "fundable" as a noun was our coinage; kept only in the advocates' phrase "the applicants it deems fundable" |
| money is tight / where it is ample | the budget is small relative to researchers' resources / where the budget is large | 5 | names the feature the paper defines instead of the idiom |
| the natural stand-in | the natural proxy | 2 | standard term |
| resource-starved (Appendix A examples) | with few resources / the researcher with the fewest resources | 2 | figurative |
| What is the schedule choice worth? Less than information | The choice of schedule is worth less than the signal provided by grant peer review | 2 | object of "information" unclear; rhetorical question removed |
| the pair of features that organizes the results | the qualitative features of the dynamics that organize the results | 2 | your wording |
| Each prediction is comparative, and each can fail. | Long version: "The three tests use data that funders already hold: grant sizes and recipients' other funding, review scores and outcomes over whole applicant pools, and the outputs of lottery-funded and review-funded grants where a funder runs both." PNAS: one test per prediction, in one sentence (the grant-size comparison within one applicant pool; the cross-field comparison scored over the whole pool; the lottery-versus-review comparison where a funder runs both). | 2 | names the analyses |

## Considered and kept

| Term | Why kept | Your call? |
|---|---|---|
| overtrust / undertrust (23) | dictionary words; each is glossed at first use ("treating a noisy signal as reliable"; "treats review as noisier than it is") | replace with the gloss everywhere if you prefer |
| tight budget / ample budget (about 50) | ordinary English for budgets; "tight" and "ample" are each defined by the budget-scale ranges in the figures | could become "small budget" / "large budget" throughout |
| no-funding output (23) | defined once as the normalization ("improvement over no funding, in percent") | alternative: "output without funding" |
| records-only funder / records-only shortfall (13) | defined operationally in the text and in the Appendix D settings table | |
| informativeness (13) | standard term for signal quality; paired with "noise" throughout | |
| heavy-tailed / evenly spread (about 90) | standard, and "spread more evenly" is your wording | |
| bottleneck (7) | standard | |
| the gap rule; the frontier; targeting; equal division; partial lottery; selection by review; center of mass | defined objects, one name each | |
| poverty-gap transfer (1) | the FGT literature's own term, cited | |
| does not catch up (2) | plain English | |

## Other edits made in the same pass

- Title, both versions: "A model of optimal science funding: targeting the capability-resource gap" (sentence case as in the file; amsart sets it in capitals on the page).
- PNAS draft moved onto your template (amsart, 11pt, one-and-a-half spacing, 4 cm margins, author block); figures placed at their sections instead of at the end. PNAS content organization unchanged: Significance, untitled introduction, Results, Discussion, Materials and Methods. 14 pages at this layout.
- Three figures spilled into the right margin (their in-plot labels sit outside the axis): Fig 7 axis width 0.82 to 0.66 of the text width, Fig 10 0.82 to 0.76, Fig 11 panels 0.46 to 0.44 each; PNAS Fig 4 likewise. The plots' heights are unchanged.
- One line in Appendix D reworded so the budget formula no longer overruns the margin.

## Second pass (same day): your six requests

| Request | What changed | Where |
|---|---|---|
| "using a Bayesian approach" | replaces "with standard Bayesian methods" | intro, both versions |
| "research output" | inserted at the first mention of output in every paragraph of the abstract, Significance, introduction, and Discussion, in the model definition, in the first paragraph of each Results subsection, and in every caption's first sentence; later mentions in a paragraph stay "output" | both versions, all captions |
| comparative metric | every number that was a percent of no-funding output is now a percent of the research output that funding adds (the funder's expected output minus the output under no funding, funder named at each use). Seed grants at three quarters: 11% (heavy-tailed, reliable review, b = 1), 17% (b = 0.2), 3.5% (default), was 6%/1%. Scheduling gain 0 to 4.5% vs review signal 7 to 27% (30 settings, not 32: the earlier count double-counted two T = 5 cells), was 0 to 0.7% vs 0.8 to 7.5%. Records-only shortfall 27% (40% heavy-tailed) of complete-information funding's gain, was 11.4% (24.6%). Review's maximum value at alpha_K = 3.5: 9% vs 39% heavy-tailed, was 2.5% vs 23%. Appendix C signal value 16/47/51% (CD, gamma = -3, Leontief), was 6.8/29.3/33.7. Table 2 scheduling value up to 6.5% at gamma = -3, was 1.54. Figs 6, 10, 12 (and PNAS Fig 4A, S2, S5) re-plotted from the sweep data on the new denominator; the ordering of the seed-grant curves by budget is now the intuitive one (smallest budget loses the largest share). Shares that were already relative (share of records-only shortfall, share of the optimal allocation's gain, retained fractions, Gini) are unchanged. | S5, S6, S7 text and footnotes; App C, App D; fig6, fig10, fig12; PNAS Results, Methods, fig-p4, SI |
| heuristics vs uniform | "each can perform worse than simply dividing the budget uniformly across all researchers" | intro, both versions |
| track record is a signal | "peer review can provide a further signal of capability" (four places); PNAS Significance "track records alone underdetermine that gap"; Results: "a noisy one, but one whose noise does not come from resources"; small grants: "on small grants output carries little information about capability" (text and Fig 3 / 2A caption), replacing "reveal almost nothing" and "uninformative" | abstract kept ("an additional, noisy signal") |
| "In practically all cases" | replaces "In every field we examine" | conclusion, both versions |
| Figure 1 | two frontiers, R = c_1 K and R = c_2 K with plotted slopes 3/2 (smaller budget b_1) and 2/3 (larger budget b_2); black = funded at both, gray = funded only at b_2, open = unfunded at both; arrows = grants at b_1 (drawing both sets of arrows cluttered the dense cluster; the caption says what the b_2 arrows would show). Labels "frontier at b_1", "frontier at b_2", and one label per group. Data: fig1-both.dat, fig1-b2only.dat, fig1-none.dat (same population as before). | fig1-frontier, fig-p1, both captions, one sentence in each text |
| panel letters | letters moved up (spacing -1.2em to -0.3em) so they clear the y-axis labels | PNAS Figs 2 to 4 |
| b = 1 is the whole field | budget = b n E[R_0] (was 2 b n E[R_0]); every b doubled: default 0.5 to 1; ranges 0.1 to 1 become 0.2 to 2; fig7 labels 1 and 6; fig10 labels 0.2, 1, 2; fig11 and fig13 settings; App D and PNAS Methods; a sentence in App D records that the code's parameter is half of the reported b. Figs 1 and 2 already used the whole-field convention, so they are unchanged. | S3, S7 footnote, App D, figs 7/10/11/13, PNAS Methods, SI |
| SI Appendix | built inside the PNAS document (to become the separate SI PDF at submission): opening paragraph and a claim-to-section table; SI Text 1 literature, 2 model with Table S1, 3 gap rule, 4 lotteries, 5 extended results (each subsection opens with the main-text sentence it supports; Figs S1 to S4), 6 production function (Fig S5, Table S2), 7 simulation specification. Main-text references are now \ref-based. Appendix C retitled "Robustness to the production function" (one term). | tex/pnas/si-*.tex, figS1..S5 |

Numbers verified in-container from the staged sweep summaries (D3_seed_signal, T_run_smooth, D_misspecified_trust, sigma_tierA_*) and a 12-cell rerun of the review-value grid (review_cells.csv, 50 populations per cell, bit-matching the earlier percent-of-no-funding values).
