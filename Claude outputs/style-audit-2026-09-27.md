# Style and clarity audit of main.tex (2026-09-27)

Scope: every paragraph of the abstract, S1-S9, Appendices A-C, and all thirteen captions. No changes made. Each item gives the location, the current text, the problem, and a proposed replacement. Problem codes: READ (a use of "read"), VERDICT, STYLE (style in the way of clarity), JARGON (a term not strictly necessary), REF (the object of a sentence could be mistaken; unpack it), GRAMMAR.

## Part 1. The two words you asked for

READ appears six times; VERDICT three times.

R1. S3, footnote on output units. "they may be read as publications, as findings, or as units of scientific value." -> "they may be interpreted as publications, as findings, or as units of scientific value."

R2. S3, signals paragraph (twice). "gives a noisy reading of each researcher's baseline resources" and "gives a noisy reading of capability" -> "gives a noisy estimate of each researcher's baseline resources" and "gives a noisy estimate of capability".

R3. S9, peer review paragraph. "The second concerns how to read the evidence on review's efficacy." -> "The second concerns how to interpret the evidence on whether review scores predict productivity."

R4. S9, limitations paragraph. "The model holds several things fixed, and its results should be read within those limits." -> "The model holds several things fixed, and its results hold within those limits."

R5. S9, limitations paragraph, last sentence. "What a reader should carry to a real field is the pair of features that organizes the results, not the percentages." (also a "this, not that" construction) -> "What applies to a real field is the pair of features that organizes the results; the percentages do not carry over."

V1. S9, practical-consequence paragraph (three uses). "No funding mechanism has a verdict of its own. Peer review, seed grants, and lotteries each receive opposite verdicts in different fields, and the two features of the field decide which verdict applies." -> "No funding mechanism is good or bad in itself. Peer review, seed grants, and lotteries each add output in some fields and lose output in others, and the two features of the field decide which."

## Part 2. Terms used throughout (decide once, apply everywhere)

G1. JARGON "tranche" (3 in main, 12 in appendices). Proposed: "the round's budget" (and $B$ where already symbolized). Example, S3: "uniform funding, which splits each round's budget equally"; S4 display: "subject to $\sum_i g_i \leq B$, the round's budget".

G2. JARGON "purse" (5, all S3). Proposed: "budget" everywhere; "deeper purses" -> "larger budgets".

G3. JARGON "efficacy" (4: S1, S2, S6, S9). It names the dispute, so one label is right, but "efficacy" is abstract. Proposed: "whether peer review predicts research outcomes" at first mention (S1) and "the predictive-accuracy dispute" or simply "the dispute over peer review" afterward. If you prefer to keep "efficacy", keep it in all four places and gloss it once in S1: "whether peer review is effective, that is, whether review scores predict research outcomes".

G4. JARGON "seeds" (9 in main, 22 in captions and footnotes): "fifty seeds per cell", "two hundred seeds per point". Simulation jargon. Proposed: "fifty simulated populations per point" (and "per cell" -> "per parameter setting" or "per point"). The single-round computation already says "two thousand populations per point", so this unifies the label.

G5. JARGON "depth", "deep grants", "standard depth", "funding depth" (8 in main, 4 in captions) alongside "thin grants". Two metaphors for one dimension (grant size). Proposed: keep "thin grants" if you like it, but replace "deep" with "large": "large grants", "grants of the default size" for "standard depth", "requires large grants, comparable to capability itself" for "requires funding depth". Or use "small grants" / "large grants" throughout and drop "thin".

G6. JARGON "production form" (4 + 1). Proposed: "production function" everywhere (the term S3 defines).

G7. JARGON "comparative static(s)" (S4, S9, and once in Appendix A's corollary title). Proposed: S4: "This dependence on the budget recurs later"; S9: "the numbers we report are comparisons across the parameter values we sweep"; Appendix A corollary title can keep the technical term.

G8. JARGON "underdetermine(s)", "underdetermination" (5: abstract, S1, S5 x3). Philosophers' term. Proposed: "do not reveal" / "does not reveal"; S5's "how much output this underdetermination loses" -> "how much output a funder loses because records do not reveal the gap".

G9. JARGON "harder complementarity" (3 + caption). Proposed: define once in S8 where the family is introduced ("the harder it is to substitute one input for the other, the stronger the complementarity") and then say "stronger complementarity".

G10. JARGON "channel", "grant-fed channel", "free channel" (8, all S7 and its footnote). Proposed rewrite in item T15 below.

G11. STYLE, three times (S1, S2, S6): "the expected value of its reduction of the funder's uncertainty regarding the gap". Abstract and hard to parse; the paper's own operational definition is simpler. Proposed, all three places: "the additional expected output of a funder allocating with the review signal over the same funder allocating without it" (S6 already has this sentence; in S1 and S2 a shorter form: "the expected output that review adds by reducing the funder's uncertainty about the gap").

G12. JARGON "epistemic capability" (4 + 1). "Epistemic" is defined at first use ("what a researcher knows and can do"). Proposed: use "epistemic capability" once, at the definition, and "capability" everywhere after, including the abstract's "epistemic capabilities" -> "capabilities".

G13. JARGON "benchmark" (3, S5-S6): "the complete-information benchmark". Proposed: "the output a funder with complete information would produce" at first use, then "the complete-information output".

G14. JARGON "calibration", "calibrated" (3). S6 footnote "review's full calibrated value" -> "review's full value when the funder weights it correctly"; S9 "The third concerns calibration." -> "The third concerns how much a funder should trust review."

## Part 3. Line-by-line items, in order

### Abstract

A1. GRAMMAR/parallelism. "researchers differ in their epistemic capabilities and in resources" -> "researchers differ in their capabilities and in their resources".

A2. REF. "this signal's value increases with capability inequality and with its informativeness": "its" could be inequality's. -> "with the signal's informativeness".

A3. REF. "even a coarse signal suffices": suffices for what. -> "even a coarse signal captures most of that value".

A4. JARGON (G8). "Track records underdetermine the gap" -> "Track records do not reveal the gap".

### S1 Introduction

I1. STYLE. "In spite of this, basic dimensions of allocation are contested." -> "Basic questions about how to allocate it are contested."

I2. REF. "while other designs find modest but real predictive power": "designs" has no noun. -> "while studies with other designs find modest but real predictive power".

I3. STYLE. "This paper develops a model of that problem and works out its logic: whom to fund, ..., what information about researchers is worth, and what egalitarian alternatives lose." -> "This paper develops a model of that problem and answers four questions: whom to fund, whether to spend early or wait and learn, what information about researchers is worth, and what seed grants and lotteries lose." ("egalitarian alternatives" is a label the reader has not met.)

I4. STYLE. Second person. "Consider the challenge of the grant funder ...: you have a fixed budget and must choose how to allocate it across researchers and research projects to create the most positive impact." -> "A grant funder, whether a federal agency disbursing billions or a private foundation disbursing millions, has a fixed budget and must decide how to divide it among researchers so as to produce the most research."

I5. STYLE. "There are two intuitive schemes. ... Each scheme can be intuitively compelling." Second sentence repeats the first. -> cut "Each scheme can be intuitively compelling."

I6. JARGON. First use of "compounds": "compounds their capability in later rounds". -> "raises their capability in later rounds; we call this compounding." (The word is then defined for the rest of the paper.)

I7. STYLE. "In the ideal condition where the funder has complete information regarding researchers' capabilities and resources" -> "When the funder knows every researcher's capability and resources".

I8. JARGON (G8). "In reality the funder lacks direct access to the gap, and track records alone underdetermine it" -> "A real funder cannot observe the gap, and track records alone do not reveal it".

I9. STYLE/REF (G11). "Its value for a Bayesian funder is the expected value of its reduction of the funder's uncertainty regarding the gap": two "its" with different referents. -> "The value of that signal to the funder is the expected output it adds by reducing the funder's uncertainty about the gap".

I10. STYLE. "The same reasoning gives what the egalitarian alternatives lose: uniform seed grants and lotteries lose the value of the targeting they displace, little where review is worth little, most where review is worth most." -> "The same reasoning says what uniform seed grants and lotteries lose: the value of the targeting they replace, which is small where review is worth little and large where review is worth most."

I11. STYLE. "The budget is a variable throughout, so our results speak across funder scales: a small foundation, whose grants are modest beside researchers' total resources, is the funder for whom targeting matters most per dollar." -> "We vary the budget throughout, so our results apply to funders of every size. Targeting matters most per dollar for a small foundation, whose grants are small relative to researchers' total resources."

I12. STYLE. "We show how strategic timing of when one provides funds matters far less than strategic targeting of whom one provides funds" -> "We show that when to fund matters far less than whom to fund".

I13. STYLE/JARGON. "Two substantive assumptions bound these results. Grants are consumed, a grant's durable residue (skills, publications, standing) being carried as capability; and capability and resources are complements throughout, neither fully substituting for the other." Absolute construction and "durable residue". -> "Two assumptions limit these results. First, grants are used up in the round they are given; whatever lasts from a grant (skills, publications, standing) is counted as capability. Second, capability and resources are complements throughout: neither fully substitutes for the other."

### S2 Related literature

L1. STYLE. "Our contribution sits at the intersection of four literatures" -> "Our model bears on four literatures".

L2. One label. "more dispersal than current practice" -> "more spreading than current practice" (the paper's word elsewhere is "spread").

L3. REF. "a dynamic that formal models show can reward questionable research practices over a career" -> "and formal models show that this self-reinforcement can reward questionable research practices over a career".

L4. STYLE. "with a tight budget, effectiveness requires concentrating funds on the researchers with the largest capability-resource gaps" -> "with a tight budget, output is highest when funds are concentrated on the researchers with the largest capability-resource gaps".

L5. JARGON (G3). "The efficacy of peer review is contested."

L6. STYLE. "Analyses with different designs find real but modest predictive power" -> "Studies with other designs find real but modest predictive power".

L7. STYLE. "a funding natural experiment finds higher-scored proposals outproduce lower-scored ones" -> "a natural experiment in funding finds that higher-scored proposals produce more than lower-scored ones".

L8. STYLE (G11). "Our model instead measures review's value: the expected value of its reduction of the funder's uncertainty regarding the capability-resource gap." -> "Our model instead measures review's value: the expected output that review adds by reducing the funder's uncertainty about the capability-resource gap."

L9. REF. "The existing case measures what a lottery saves": "case" is ambiguous (argument? instance?). -> "The existing formal argument for lotteries measures what a lottery saves".

L10. STYLE. "What has not been measured is what a lottery forgoes. Our model measures it: the value of targeting, which runs from negligible, where budgets are ample, signals uninformative, or capability evenly spread, to substantial, where capability is heavy-tailed and budgets tight (S8)." Long, with an inverted list. -> "What has not been measured is what a lottery gives up: the value of targeting. Our model measures that value. It is negligible where budgets are ample, review is uninformative, or capability is evenly spread, and substantial where capability is heavy-tailed and budgets are tight (S8)."

L11. JARGON. "simulates funding strategies on an epistemic landscape" -> add a gloss: "on an epistemic landscape, a map of research problems of varying value".

L12. JARGON. "treat allocation as a decision problem under heavy-tailed uncertainty about project value, pairing an empirical bibliometric signal with a biased-lottery mechanism" -> "model allocation as a decision under heavy-tailed uncertainty about project value, and combine a citation-based signal with a lottery weighted by that signal".

L13. GRAMMAR. "Our formulation and analysis of the capability-resource gap allows us" -> "allow us".

L14. STYLE. "when lotteries are defensible (S8)" -> "when lotteries lose little (S8)" (matches what S8 measures).

### S3 Our model

M1. STYLE. "a correlation parameter rho allows capable researchers to tend to be well resourced" -> "a correlation parameter rho lets capability and resources be correlated across researchers".

M2. READ (R1).

M3. STYLE/JARGON. "The harmonic mean sits between two boundary cases that fail in opposite directions: Cobb-Douglas, under which a capable researcher can compensate almost fully for vanishing resources, and Leontief, under which no compensation is possible at all." "Fail" is unexplained, and the two names are economics jargon. -> "The harmonic mean sits between two boundary cases, each too extreme for our purposes: Cobb-Douglas production (output is the geometric mean of the two inputs), under which a capable researcher can compensate almost fully for vanishing resources, and Leontief production (output is the smaller of the two inputs), under which no compensation is possible at all."

M4. STYLE. "This is a substantive assumption, not bookkeeping: whatever endures from a grant, the skills, publications, and standing it produces, is carried by the capability update, so funding's lasting effect flows through what researchers become rather than what they keep. Science is in this respect an accumulation process: success begets the means of further success." The last two clauses are flourishes. -> "This is a substantive assumption: whatever lasts from a grant (skills, publications, standing) enters the model as capability, so a grant's lasting effect is the capability it builds."

M5. JARGON (G2). "purse", three times in the budget paragraph and once in the parameter-range paragraph.

M6. READ (R2), twice.

M7. JARGON. "allocate on posterior expected marginal returns" -> "and allocate each dollar to the researcher whose expected output, given the funder's current beliefs, rises most with it".

M8. JARGON (G1). "which splits each round's tranche equally".

M9. STYLE. "Our conclusions rest either on analytic results supported by proofs or on simulation results that are established across the following parameter ranges." -> "Each of our conclusions rests on a proof or on simulations across the following parameter ranges."

M10. STYLE. "our range covers this territory and extends to markedly more equal fields" -> "our range covers these values and extends to much more equal fields".

M11. JARGON. "with common random numbers across strategies, so every comparison between strategies is paired" -> "using the same random draws for every strategy, so that strategies are compared on identical populations and signals".

M12. STYLE. "With the model in place, the funder's problem is exact" -> "is well defined".

### S4 The optimal allocation

O1. STYLE. "Further, we initially restrict our attention to a single round of funding, temporarily putting aside the contribution of resources to the accumulation of researchers' capabilities across rounds." -> "and to a single round, setting aside for now the growth of capability across rounds".

O2. JARGON (G1). "subject to sum g_i <= the round's tranche" -> "<= B, the round's budget".

O3. REF. "when marginal values are equal across all funded researchers and no unfunded researcher's marginal value exceeds theirs" -> "exceeds that common value".

O4. STYLE. Proposition: "where the constant c > 0 is set by spending the budget" -> "where the constant c > 0 is the value at which the grants sum to the budget".

O5. REF. "a productive researcher whose resources approach or exceed their target has a small or negative gap, and an additional dollar adds little output there" -> "adds little to their output".

O6. STYLE. "Funding the under-resourced fails for the mirror-image reason" -> "for the opposite reason".

O7. JARGON (G7). "This comparative static recurs: it sets how much output is lost by overriding the optimal allocation with seed grants and lotteries (S8), and it determines when optimal funding concentrates rather than spreads (S8)." -> "This dependence on the budget recurs later: it sets how much output is lost when seed grants and lotteries override the optimal allocation, and it determines when optimal funding concentrates rather than spreads (both S8)."

### S5 The track record

R6. JARGON. "Cross-sectional output, however carefully recorded, underdetermines what the funder most needs to know." -> "Output observed in one round, however carefully recorded, does not reveal what the funder most needs to know: the gap."

R7. JARGON (G8). "Our model measures how much output this underdetermination loses." -> "Our model measures how much output a funder loses because records do not reveal the gap."

R8. JARGON (G13). "its output falls well short of the complete-information benchmark" -> "its output falls well short of what a funder with complete information would produce".

R9. STYLE. "This shortfall is the span from records alone to complete information. Next we ask how much of that span a review signal recovers, and there display this records-only baseline beside a funder holding the review signal." -> "This shortfall is the most that any further information could recover. In the next section we ask how much of it a review signal recovers, and we show the records-only funder beside the funder with the review signal."

R10. JARGON (G5). "separating capability from resources by observing output requires funding depth, grants comparable to capability itself" -> "requires large grants, comparable to capability itself".

R11. JARGON (G5) and one label. "a funder allocating at standard depth gains almost nothing over uniform funding, while the same funder with deep grants recovers most of the value of discrimination" -> "a funder giving grants of the default size gains almost nothing over uniform funding, while the same funder giving much larger grants recovers most of the value of targeting". ("value of discrimination" is a second label for the value of targeting.)

R12. STYLE. "For nascent and resource-poor fields" -> "For new and resource-poor fields".

R13. REF. "That is what peer review can provide under the right conditions." -> "Peer review can provide such a signal under the right conditions."

### S6 Peer review

V1. STYLE (G11). Opening paragraph, two sentences saying the same thing. -> "The value of review is the additional expected output of a funder allocating with the review signal over the same funder allocating without it."

V2. STYLE. "recovers most of the output that a funder relying on records alone leaves unrealized" -> "recovers most of the output that a funder relying on records alone fails to obtain".

V3. STYLE. Footnote. "a small but statistically solid reversal" -> "a small but statistically significant reversal". And: "The cause is the funder's other estimate: resources are also observed with noise, and a slightly noisy capability signal tempers the allocation against errors in the resource estimate, while a near-perfect one commits the allocation to them." -> "The cause is the funder's estimate of resources, which is also noisy. A slightly noisy capability signal makes the funder hedge against errors in its resource estimate; a near-perfect capability signal makes the funder act on those errors."

V4. STYLE. "Our main result delineates how much review is worth as a function of the field" -> "Our main result is how much review is worth in different fields".

V5. JARGON (G6). "does not depend on our production form" -> "on our choice of production function".

V6. STYLE/REF. "A corollary concerns where review effort goes. ... fine distinctions among the many applications of comparable middling merit add comparatively little output. A review system that invests most of its effort in exactly those fine distinctions may be misallocating it." -> "A consequence concerns where review effort should go. ... fine distinctions among the many applications of comparable merit add little output. A review system that spends most of its effort on those fine distinctions may be misallocating that effort."

V7. REF. "The remaining comparisons run the other way." Which comparisons. -> "In the other fields, review is worth less."

V8. JARGON (G3). "the dispute over the efficacy of peer review".

V9. REF. "One condition qualifies all of this" -> "One condition qualifies these results".

V10. JARGON (G14). Footnote: "loses more than review's full calibrated value" -> "loses more than review's full value when weighted correctly".

### S7 Timing

T1. STYLE. "between rounds two things move" -> "between rounds two things change".

T2. STYLE. "this section asks whether it should not: whether to spend early, evenly, or late, and how much the choice matters" -> "this section asks whether it should spend early, evenly, or late, and how much the choice matters".

T3. STYLE. "Researchers without resources produce nothing; where nothing is produced, nothing is observed and nothing compounds; so money held for later rounds produces neither evidence nor growth in the meantime." -> "A researcher without resources produces nothing, so the funder observes nothing and no capability grows; money held for later rounds therefore produces neither evidence nor growth in the meantime."

T4. STYLE. "should spend at once, to start the engine" -> "should spend at once".

T5. STYLE. "The argument fails on inspection." -> "The argument fails."

T6. STYLE. "exactly as blind as one that spends everything in the last" -> "and knows exactly as little as a funder that spends everything in the last round".

T7. REF. "That favors putting some money early." -> "This complementarity favors putting some money early."

T8. JARGON. Footnote: "In a planner-free check over a grid of two-block schedules executed by a funder re-deciding each round" -> "In a check that fixes the schedule in advance (two blocks of rounds, each with its own share of the budget) and lets the funder re-decide only whom to fund each round".

T9. REF. "Our simulations bear this out and say which way the asymmetry runs: observation comes first, money after." "The asymmetry" has not been named. -> "Our simulations bear this out and say which should come first: observation, then money."

T10. STYLE. "at more than a token rate" -> "at more than a negligible rate". "spending late is rewarded because late money is informed money, and because early output compounds on its own" -> "spending late pays because money spent later is spent with more information, and because early output compounds on its own".

T11. JARGON (G5). "Deep grants change this." -> "Large grants change this." "As grants approach the scale of capability itself" fine.

T12. STYLE. "on a schedule that puts observation before money: a small early share to make output observable, the mass deployed once informed" -> "on a schedule that spends a small share early, so that output can be observed, and the rest once the funder has learned who is capable".

T13. JARGON (G5). Footnote: "at standard depth", "six times that depth", "twelve times" -> "at the default grant size", "six times that size", "twelve times that size".

T14. JARGON (G10). Whole paragraph. "The pattern has a single main cause. Capabilities compound through all research output, funded or not; grants add to output and so add to compounding. Separating the two channels: with only the grant-fed channel active, the optimal schedule leans mildly early, since earlier grants compound longer. With only the free channel active, it leans strongly late. The free channel dominates once it reaches even a fraction of the grant-fed rate, so in any field where research compounds regardless of who is funded, the funder's reason to wait outweighs its reason to hurry." -> "The pattern has a single main cause. Capability grows with all research output, whether or not a grant paid for it. We separate the two sources of growth: growth from the output that grants add, and growth from the output researchers produce anyway. With only the first source active, the optimal schedule leans slightly early, since earlier grants have longer to compound. With only the second source active, the schedule leans strongly late. The second source dominates once its rate reaches even a fraction of the first's, so in any field where research compounds regardless of who is funded, the funder's reason to wait outweighs its reason to hurry." The footnote's "free channel" and "grant-fed channel" follow the same relabeling.

T15. STYLE. "the gain from deliberate scheduling over even installments re-decided each round" -> "the gain from choosing the schedule, over spending equal installments each round,".

T16. STYLE. "and saturates beyond them" -> "and levels off beyond five rounds".

T17. JARGON (G6). "as in our production form" -> "as in our production function".

### S8 Spreading funds

E1. STYLE. Dense sentence. "For a funder allocating with records and review, targeting adds more than uniform funding's entire gain over no funding in a heavy-tailed field with sharp review and a tight budget, and about a quarter of uniform funding's gain in the default field at the largest budget we sweep." -> "For a funder allocating with records and review: in a heavy-tailed field with sharp review and a tight budget, targeting adds more output than uniform funding's entire gain over no funding; in the default field at the largest budget we sweep, targeting adds about a quarter of that gain."

E2. STYLE. Inverted list. "under the same conditions that make review worth little: evenly spread capabilities, large budgets, uninformative peer review." -> "under the same conditions that make review worth little: where capability is evenly spread, the budget is large, or review is uninformative."

E3. STYLE. "behind both lies the older question of whether funding should be concentrated" -> "and both raise the older question of whether funding should be concentrated".

E4. STYLE. "To test this hope we embed" -> "To test this, we embed".

E5. JARGON (G9). "harder complementarity" three times; define once at the family's introduction.

E6. REF. "a pattern consistent with a funder that raises every researcher toward their bottleneck once the budget allows it" -> "a pattern consistent with a funder that, once the budget allows, brings every researcher's resources up toward the level their capability can use".

### S9 Discussion

D1. REF. "Researcher track records are modestly predictive" -> "modestly predictive of future output".

D2. STYLE. "so both lose most output where that value is largest: tight budgets and heavy-tailed capability" -> "where that value is largest, in fields with tight budgets and heavy-tailed capability".

D3. STYLE ("this, not that"). "the model's answer is therefore thin grants to many researchers, not a lottery among a fundable few" -> "the model's answer is thin grants to many researchers; a lottery among a fundable few produces less."

D4. VERDICT (V1).

D5. STYLE. "The funder's budget relative to researchers' total resources is the funder's own arithmetic." -> "The funder can compute the first feature, its budget relative to researchers' total resources, from its own budget and what it knows about the field."

D6. JARGON. "The distribution of output is the natural proxy for the distribution of capability" -> "The natural stand-in for the distribution of capability is the distribution of output".

D7. READ (R3).

D8. REF/STYLE. "Both kinds of study measure prediction after review has already made the separation the model values, inside the band where the model expects review to add least." -> "Both kinds of study measure prediction among applicants that review has already selected, that is, after review has made the separation the model values most, and within the range of applicants where the model expects review to add least."

D9. JARGON. "is conservative" -> "understates review's informativeness".

D10. JARGON (G14). "The third concerns calibration." -> "The third concerns how much a funder should trust review."

D11. READ (R4).

D12. JARGON. "Output is a scalar rate." -> "Output is a single number per round."

D13. JARGON. "Projects have no heterogeneity and exploration has no value" -> "Projects do not differ in kind, and exploring new directions has no value in the model".

D14. STYLE. "The review signal is unbiased noise." -> "The review signal is noisy but unbiased."

D15. STYLE. "There is no entry, exit, or career structure." -> "Researchers do not enter or leave the field, and have no careers."

D16. JARGON (G7). "the magnitudes are orderings and comparative statics over the parameters we sweep" -> "the numbers we report are comparisons across the parameter values we sweep".

D17. READ and "this, not that" (R5).

D18. JARGON. "capability proxies, whether prior record or review score" -> "measures that stand in for capability, such as prior record or review score". "review's predictive validity" -> "how well review scores predict output".

### Appendices

P1. STYLE. Appendix A: "Lemma 1 states the economics of the problem." -> "Lemma 1 summarizes the structure of the problem."

P2. JARGON (G1). "tranche" throughout Appendix A (twelve uses, including the examples: "Track-record funding sends most of the tranche", "Scarcity alone misdirects the tranche"). -> "the round's budget" / "the budget".

P3. Appendix C: "the substitutable end of the family" is fine; "the lean toward late spending" fine.

### Captions

C1. Fig 7. "Observation first, money after." slogan -> "The funder spends a small share first and the rest once it has observed output."

C2. Fig 8. "Compounding, not poverty, sets the schedule's tilt." ("this, not that"; "tilt") -> "The compounding rate, not resource poverty, sets how late the schedule leans." If you want the construction gone entirely: "The compounding rate sets how late the schedule leans; resource poverty changes it only slightly."

C3. Fig 9. "over one that spends even installments and re-decides each round" fine. "an order of magnitude below or more at the default rate and below" -> "at least an order of magnitude below at the default rate and lower rates".

C4. Fig 12. "across the CES family" -> "across the family of production technologies (the CES family)".

C5. Fig 13. "harder complementarity" (G9). "toward perfect complements (Leontief, detached marks at the left edge)" fine.

C6. Figs 4, 5, 6, 8, 10, 12, 13 and the S3 parameter footnote: "seeds" (G4).

C7. Fig 10. "A review-informed funder gives out a fraction of its budget as equal grants to all researchers and targets the rest." fine.

## Part 4. Items I looked at and left alone

The following recurring terms are technical but necessary and consistently used: "capability-resource gap", "gap rule", "funding frontier", "compounding rate" (once defined per I6), "heavy-tailed", "Bayesian funder", "myopic" and "forward-looking", "records-only", "review-informed", "partial lottery", "equal division", "selection by review", "Gini coefficient", "Cobb-Douglas" and "Leontief" (once glossed per M3). "Thin grants" is a metaphor but a clear one; the decision in G5 is whether to keep it.

Count: 6 READ, 3 VERDICT, 14 global terms, about 95 local items.
