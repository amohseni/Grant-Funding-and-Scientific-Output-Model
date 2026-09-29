# Introduction: function, hooks, and two refined versions (2026-09-29)

Target: PNAS. Audiences: the general PNAS reader, funders (including the Templeton
program officers), and metascience researchers. No changes made to main.tex.

## 1. What the introduction has to do

For this paper and these readers, the introduction has five jobs, in this order.

1. Establish in one or two sentences that the problem is large and unsettled. The
   reader should know by the end of the first paragraph that funders face live
   disputes (concentrate or spread; whether review predicts; whether to use
   lotteries) and that these are not settled by the evidence.
2. Say why the evidence has not settled them. This is the sentence the current draft
   lacks and the one that justifies a model: the empirical studies measure how well
   review scores predict outcomes among funded grants, and the concentration
   literature measures returns to funding, but none of them can say what a funder
   should do, because that depends on what the funder is uncertain about and what its
   information is worth. A decision-theoretic model is the instrument for that
   question.
3. State the contribution accurately and in a measured way: a model in which funding
   is a sequential decision under uncertainty about two attributes of researchers,
   analyzed with standard Bayesian decision theory, from which the answers to the
   three disputes are derived and shown to depend on the same two features of a
   field. The claim is a theoretical one (this is how the problem is structured and
   this is what follows), not an empirical one about any funder.
4. Give the results at full strength, briefly, so a reader who stops here knows what
   the paper found.
5. Name the two assumptions that bound the results.

The current draft does 1, 4, and 5 well. It does 3 only implicitly ("develops a model
of that problem") and does not do 2 at all, so the reader is told that disputes exist
and then handed a model without being told why a model is the right tool. The refined
version below adds one paragraph for jobs 2 and 3.

For funders specifically, the constructive framing they need is already the paper's
finding: no mechanism is good or bad in itself; the two features of the field decide.
That sentence belongs in the introduction's results paragraph, because it is what a
program officer will act on.

## 2. Format facts that change the introduction's shape

PNAS research reports are about 4,000 words of main text at the preferred six pages
(twelve maximum), with an untitled introduction, then Results, Discussion, and
Methods, a 250-word abstract, a 120-word Significance Statement, and an SI Appendix
that does not count toward the page limit. [check: against pnas.org author center;
the numbers above are from a secondary summary dated 2026-09-21.]

Two consequences. First, the paper as it stands (8,400 words of body, nine sections,
four appendices) is a different object from a PNAS paper: a PNAS version would carry
the model, the three main results, and the discussion in about 4,000 words, with the
derivations, the full sweeps, and the robustness work in the SI Appendix. Second, the
PNAS introduction is untitled and short (three to five paragraphs), it absorbs the
related-literature material, and it has no roadmap paragraph. So there are two
introductions to write: the one for the current full-length format, and the one for
the PNAS format. Both are below. The current-format one is a refinement of your
draft; the PNAS one is a compression of it.

## 3. Hooks: five first lines

H1, stakes first. "The United States spent close to a trillion dollars on research and
development in 2024, roughly a fifth of it funded by the federal government, and
how best to allocate that money remains unsettled." (The source gives a fifth as federally funded, not as grants; the hook keeps the source's claim.) Concrete, accurate, and it puts
the number in the first line rather than the second. Best for the general PNAS reader.

H2, the decision first. "A research funder has a fixed budget, a population of
researchers it knows only through their past output and the judgment of reviewers,
and one question: who should get the money?" Puts the reader inside the problem
immediately. Best for funders. Risk: the informal "who should get the money" reads as
a flourish to some referees.

H3, the result first. "The best grant a funder can give a researcher depends on two
things, how capable the researcher is and how much they already have to work with,
and the two rules funders most often use each track only one of them." Findings-first,
which is the paper's own rule for introductions. Risk: a reader who has not yet seen
the disputes does not know why this matters.

H4, the model first. "Allocating research funding is a decision under uncertainty:
the funder must divide a fixed budget among researchers whose capabilities and needs
it cannot observe, using only their records and the judgment of reviewers." Accurate
and measured; states the paper's stance in the first line. Slightly dry for a first
line.

H5, the disputes first. "Whether research funding should be concentrated or spread,
whether peer review predicts anything, and whether lotteries should replace it are
three disputes that share one unstated question: what a funder with a fixed budget
and imperfect knowledge of researchers should do." Dense but content-bearing. Risk: a
54-word first sentence.

Recommendation: H1 as the first sentence, H4's content as the second paragraph's
opening (the "why a model" paragraph), and H3's content as the first sentence of the
results paragraph. That order gives stakes, then the instrument, then the finding,
each in the place where the reader is ready for it.

## 4. Refined introduction, current format

Science funding is large and consequential, and how best to allocate it remains unsettled. The United States spent close to a trillion dollars on research and development in 2024, roughly a fifth of it funded by the federal government \citep{NCSES2026, NSB2026}. Whether that money should be concentrated on a few researchers or spread across many is a long-running dispute \citep{Aagaard2020}. How well peer review predicts research outcomes is contested: among funded NIH grants, percentile scores predicted subsequent productivity barely better than chance \citep{Fang2016}, while studies with other designs find modest but real predictive power \citep{Gallo2014, LiAgha2015, ParkLeeKim2015}. And whether review should be replaced, in part, by lotteries has moved from proposal to practice \citep{Roumbanis2019, Heyard2022}.

These disputes concern a single underlying problem: what a funder with a fixed budget and imperfect knowledge of researchers should do. The empirical evidence cannot settle them on its own, because each study measures one input to that decision, how well scores predict outcomes among funded grants or how output responds to grant size, and none can say what the funder should do with what it knows. That question depends on what the funder is uncertain about, on what its information is worth, and on how the answer changes with the budget. It is a question for decision theory. This paper treats the funder's problem as a sequential decision under uncertainty and analyzes it with standard Bayesian methods. The model is deliberately simple, and its contribution is to our theoretical understanding: it shows how the three disputes are structured, derives the optimal allocation, and identifies the two features of a field that decide what peer review is worth, what seed grants and lotteries lose, and whether funds should be concentrated or spread.

Consider the challenge of the grant funder, whether a federal agency disbursing billions or a private foundation disbursing millions: you have a fixed budget and must choose how to allocate it across researchers and research projects to create the greatest positive impact. Two intuitive schemes might be considered. Fund the track record: give to those who have produced the most, since past productivity is evidence of future productivity. Fund the under-resourced: give to those who have the least, since a dollar goes further where dollars are scarce. Both, we will show, are suboptimal, and for the same reason: each accounts for only one of the two quantities the impact of a grant depends on.

We model this decision situation as follows. A funder with a fixed budget allocates grants, over multiple rounds, to researchers who differ in two ways: epistemic capability, what a researcher knows and can do, and resources, what a researcher has to work with. Research output requires both and is limited by the scarcer of the two; the funder can observe output, but not capability or resources themselves. A grant adds to a researcher's resources in the round it is given and, because doing research builds expertise, raises their capability in later rounds. Each round, the funder observes the output that results, updates its beliefs about the researchers, and allocates again.

Our results are as follows. When the funder knows every researcher's capability and resources, the optimal allocation funds the capability-resource gap: each funded researcher receives a grant increasing in their capability and decreasing in their existing resources, and the budget determines how many researchers receive funding at all. Each intuitive scheme fails because it accounts for only one of the two quantities. In reality the funder lacks direct access to the gap, and track records alone underdetermine it: we show that a funder relying on the track record performs far worse than one that estimates the gap.

Peer review supplies an additional signal of capability. The value of that signal to the funder is the expected output it adds by reducing the funder's uncertainty about the gap, and that value is greatest where capability is unequally distributed and review informative; in fields with heavy-tailed distributions of capability, even a coarse signal captures most of this value. The same reasoning says what uniform seed grants and lotteries lose: the value of the targeting they replace, which is small where review is worth little and large where review is worth most. We vary the budget throughout, so our results apply to funders of every size, and the budget is one of the two features that decide every result: targeting matters most per dollar for a small foundation, whose grants are small relative to researchers' total resources, and whether optimal funding should concentrate on a few researchers or spread across many depends on the budget relative to those resources. No funding mechanism is good or bad in itself; the funder's budget and the field's distribution of capability decide which mechanism adds output. Finally, when to fund matters far less than whom to fund: the best split of the budget between early and late rounds shifts with the funder's budget, but the difference is small relative to the difference between funding the right and the wrong researchers.

Two assumptions limit these results. First, grants are used up in the round they are given; whatever lasts from a grant (skills, publications, standing) is in the form of naturally compounding researcher capability generated by the chance to produce more research. Second, capability and resources are complements throughout: neither fully substitutes for the other. We defend each assumption where it is used.

We proceed as follows. [roadmap unchanged]

Changes from your draft, so you can accept or reject each: (a) first sentence fixed ("how best do it it") and the number moved into the first line; (b) new second paragraph doing jobs 2 and 3; (c) "what seed grants and lotteries cost" kept as "lose" (the paper's term after the price-language sweep; your call); (d) the results paragraph gains the funder-facing sentence "No funding mechanism is good or bad in itself..." and the sentence on the budget as one of the two deciding features; (e) "Each intuitive scheme fails because it accounts for only one of the two quantities" now appears once, in the results paragraph, since paragraph three already says it. Everything else is your text.

## 5. Introduction, PNAS format (untitled, about 480 words)

The United States spent close to a trillion dollars on research and development in 2024, roughly a fifth of it funded by the federal government \citep{NCSES2026, NSB2026}, and how best to allocate that money remains unsettled. Whether funding should be concentrated on a few researchers or spread across many is a long-running dispute \citep{Aagaard2020}. How well peer review predicts research outcomes is contested: among funded NIH grants, percentile scores predicted subsequent productivity barely better than chance \citep{Fang2016}, while other studies find modest but real predictive power \citep{LiAgha2015, ParkLeeKim2015}. And several funders now allocate part of their budgets by lottery \citep{Liu2020, Heyard2022}.

These disputes concern one problem: what a funder with a fixed budget and imperfect knowledge of researchers should do. The empirical studies each measure one input to that decision, and none can say what the funder should do with what it knows; that depends on what the funder is uncertain about, what its information is worth, and how the answer changes with the budget. We treat the funder's problem as a sequential decision under uncertainty and analyze it with standard Bayesian methods. Researchers differ in capability, what they know and can do, and in resources, what they have to work with; output requires both and is limited by the scarcer; the funder observes output, and the noisy judgments of reviewers, but neither attribute directly. A grant adds resources in the round it is given and, because doing research builds expertise, raises capability in later rounds.

Three results follow. First, the optimal allocation funds the gap between each researcher's capability and their resources: grants increase with capability and decrease with existing resources, and the budget sets how many researchers are funded at all. The two rules funders most often use, funding the track record and funding the under-resourced, each track one of the two quantities and can each produce less output than dividing the budget equally. Second, track records underdetermine the gap, and peer review supplies what they cannot, a signal of capability. The value of that signal, what seed grants and lotteries lose by giving up targeting, and whether funds should be concentrated or spread all depend on the same two features of a field: how much of researchers' resources the funder's budget can supply, and how unequal capability is. Where the budget is tight and capability heavy-tailed, targeting is worth a great deal, even coarse review captures most of that value, and lotteries lose much of it; where the budget is ample and capability evenly spread, little is at stake and spreading funds widely loses almost nothing. No funding mechanism is good or bad in itself. Third, when to fund matters far less than whom to fund: the gain from choosing the schedule of spending is smaller than the gain from a review signal in every setting we examine.

Two assumptions bound these results: grants are used up in the round they are given, with whatever lasts from them entering as capability, and capability and resources are complements. The model is deliberately simple. Its contribution is to the theoretical understanding of the funder's problem: it shows that the three disputes turn on two features of a field that a funder can estimate, and it says which way each dispute goes in each kind of field.

## 6. Significance Statement (PNAS, 120 words), draft

Research funders decide how to divide fixed budgets among researchers they know only through past output and reviewers' judgments. Whether to concentrate funds or spread them, whether peer review is worth its cost, and whether lotteries should replace it are contested, and the empirical evidence has not settled them. We model funding as a sequential decision under uncertainty and show that the optimal allocation funds the gap between each researcher's capability and their resources, that track records alone cannot reveal that gap, and that what peer review is worth, what lotteries and seed grants lose, and whether to concentrate or spread all turn on two features of a field: the funder's budget relative to researchers' resources, and how unequal capability is. [118 words]

## 7. Decisions for you

1. Hook: H1 (recommended), or another of the five.
2. The new "why a model" paragraph: keep, shorten, or cut.
3. "lose" or "cost" for what seed grants and lotteries give up.
4. Whether to prepare the PNAS-format version now (main text about 4,000 words, SI
   Appendix carrying the derivations, sweeps, and robustness), or to finish the
   full-length version first and compress afterward. If PNAS is the target, I
   recommend writing the PNAS version now, since the compression decides what the
   introduction must carry.
