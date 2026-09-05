# Section 3 (Our model), draft v1 (2026-08-14)

Drafted plain-first. Sources: the model write-up with errata E1-E7 applied (canonical
notation; corrected budget; sweep-base defaults); strategy set compressed to the contrasts
the paper uses; full specification deferred to the appendix. Decision flags for Aydin at the
bottom. No formula beyond the two the section needs (output intensity; capability update).

---

## 3. Our model

A funder faces a population of $n$ researchers over $T$ funding rounds. Researcher $i$ has
epistemic capability $K_i$ and resources $R_i$. Initial capabilities and resources are drawn
from Pareto distributions with tail parameters $\alpha_K$ and $\alpha_R$; smaller tail
parameters produce more unequal fields, and a correlation parameter $\rho$ allows capable
researchers to tend to be well resourced. The funder knows these distributions; it never
observes any researcher's capability or resources directly.

In each round, researcher $i$ produces research output at random, in discrete units, with
expected output
$$\lambda_i = A\,\Lambda(K_i, R_i),$$
where $A$ is a productivity constant and $\Lambda$ is a production function that increases
in both arguments but is limited by the scarcer of the two. Our main analysis takes
$\Lambda$ to be the harmonic mean:
$$\lambda_i = \frac{2A K_i R_i}{K_i + R_i}.$$
The harmonic mean makes capability and resources complements without making either an absolute
requirement: a capable researcher with thin resources still produces, by working more slowly
or borrowing equipment, but no surplus of one input fully substitutes for a shortage of the
other. The harmonic mean sits between two boundary cases that fail in opposite directions:
Cobb-Douglas, under which a capable researcher can compensate almost fully for vanishing
resources, and Leontief, under which no compensation is possible at all. \S9 repeats our
analysis across the full family between these boundaries.

Capabilities increase through research. Between rounds,
$$K_i \leftarrow K_i + \epsilon\,\Lambda(K_i, R_i),$$
where $\epsilon$ is the compounding rate: doing research builds expertise, and the same
bottleneck that limits output limits growth. Growth requires both inputs, and capability
never declines. Resources, by contrast, do not accumulate: resources in a round are baseline
plus grant, $R_i = R_{i0} + g_i$. The baseline $R_{i0}$ recurs every round, representing the
time and support a researcher has independent of the funder, while grants are consumed in
the round they are given. This
is a substantive assumption, not bookkeeping: whatever endures from a grant, the skills,
publications, and standing it produces, is carried by the capability update, so funding's
lasting effect flows through what researchers become rather than what they keep. Science is
in this respect an accumulation process: success begets the means of further success.

The funder's total budget is fixed before any outcomes are observed. We parameterize it by
a budget scale $b$: at $b = 0.5$, the purse equals the community's expected baseline
resources for one round, and the purse scales linearly in $b$. (The exact normalization is
in Appendix [X].) The budget does not grow with the horizon: a longer horizon varies the timing of spending, never its amount. Strategies that
do not plan ahead spend an equal share $B/T$ each round. Forward-looking strategies may
spend on any schedule: each round, they re-plan the entire remaining budget over the
remaining horizon and execute the current round's allocation.

We stay agnostic about the interpretation of output units: they may be read as
publications, as findings, or as units of scientific value. What a researcher has produced
to date is their track record.

The funder learns as it goes. It observes each researcher's track record as it accumulates,
and two further signals are available. A resource signal, observed once at the
start, gives a noisy reading of each researcher's baseline resources (noise $\tau_R$), as
one might obtain from institutional information. A review signal, drawn afresh in each
round that a strategy uses it, gives a noisy reading of capability (noise $\tau_K$), as one
might obtain from proposal evaluation and expert judgment. The review signal's
informativeness is governed by $\tau_K$: the smaller the noise, the more informative the
review. Bayesian strategies update beliefs about every researcher after each round and
allocate on posterior expected marginal returns.

We compare seven funding strategies, which differ along two dimensions: information, with
or without the review signal, and planning, myopic or forward-looking. A myopic strategy
optimizes expected output in the current round; a forward-looking strategy optimizes
expected output over the whole remaining horizon. Three non-Bayesian baselines complete the
set: no funding, against which every result is normalized; uniform funding, which splits
each round's tranche equally; and track-record funding, which allocates in proportion to
the previous round's output.

| Strategy | Information | Planning |
|---|---|---|
| No funding | none | none |
| Uniform | none | none |
| Track record | track record | none |
| Myopic, without review | track record + resource signal | current round |
| Forward, without review | track record + resource signal | remaining horizon |
| Myopic, with review | track record + resource signal + review | current round |
| Forward, with review | track record + resource signal + review | remaining horizon |

Uniform floors and lotteries, which modify these strategies rather than add to them, are
introduced in \S8.

Our conclusions rest on parameter sweeps, not on a single setting. The tail parameters
$\alpha_K$ and $\alpha_R$ range from 1.3 to 3.5: empirical distributions of research
productivity are heavy-tailed, with tail exponents near 2 [cite: Lotka-type evidence;
check], so the range spans fields markedly more unequal than that evidence suggests and
fields markedly more equal. The budget scale $b$ ranges from 0.1 to 1, from a funder whose
purse is small beside the community's own resources to one able to relieve most
researchers' bottlenecks; \S5 and \S7 also examine deeper purses. The compounding rate
$\epsilon$ ranges from near zero, where funding buys only current output, to 0.85, where
this round's funding substantially raises later capability. The review noise $\tau_K$
ranges from 0.05, nearly perfect review, to 20, nearly uninformative review, a range wide
enough to contain the noise implied by empirical estimates of review's predictive accuracy
(\S6). The correlation $\rho$ ranges from $-0.5$ to 0.8, and the horizon $T$ from one
round to ten.

[Placement: caption of the first figure.] Unless a figure states otherwise, unstated
parameters take the values $T = 2$, $n = 50$, $\epsilon = 0.1$, $b = 0.5$,
$\alpha_K = \alpha_R = 2$, $\tau_K = \tau_R = 1$, and $\rho = 0$.

Randomness enters our model through the initial population draws, the signals, and realized
output. Realized output shapes the funder's beliefs and, through them, its allocations.
Capability growth, by contrast, follows expected output: the update above uses $\Lambda$
itself, so a researcher's accumulation depends on their capability and resources, not on
luck in any given round. This is a modeling choice: it treats a round as containing enough
separate acts of research for luck to average out of accumulation, while leaving the
funder's inference fully exposed to noise. For the same reason, we score each strategy by
the sum of expected outputs along its run, which equals realized output averaged over
repetitions of the same allocations; this removes simulation noise from comparisons without
favoring any strategy. Each configuration runs on 200 independently drawn populations, with
common random numbers across strategies, so every comparison between strategies is paired.
Throughout, we report each strategy's improvement over no funding, in percent. The full
specification, parameter table, and code are in Appendix [X].

With the model in place, the funder's problem is exact: choose grants, given beliefs, to
maximize expected total output over the horizon. \S4 solves the case in which the funder
has complete information regarding researchers' capabilities and resources.

---

## Notes for Aydin

1. RESOLVED (Aydin, 2026-08-14): Zollman reference cut (PhilSci archive only); the
   accumulation sentence stays uncited.
2. Track-record funding (the code's naive strategy S3) is included in the strategy set so
   that \S4's refutation of "fund the track record" has a measured counterpart. The
   "fund the under-resourced" option has no implemented counterpart; \S4 refutes it
   analytically. Flag if you want this asymmetry noted there.
3. Strategy names replace the code's S1-S9 labels; the appendix will carry the mapping
   (S1 no funding, S2 uniform, S3 track record, S4/S7 without review, S5/S8 with review,
   S6/S9-S11 seeded).
4. To the appendix: Poisson output, prior structure and posterior sampling (M = 400),
   productivity normalization A = 1/2 under the family convention, k_min = r_min = 1,
   validation gates, budget_ref option, and the exact budget normalization B = 2bnE[R]
   (the 2 preserves the two-round purse convention; moved out of the main text 2026-08-14).
   Nothing in the main text depends on them.
5. The two displays are the only formulas; the gap rule itself debuts in \S4 per your edict.
6. Baseline-vs-initial verified in model.R (2026-08-14, Aydin's query): B_total = 2 b n E_ref
   with E_ref the THEORETICAL Pareto mean of the resource distribution, computed from
   parameters before the population is drawn (lines 1485-1490); R0 fixed at its initial draw
   and recurring per round (R = R0 + g everywhere); no resource accumulation is implemented,
   so the purse cannot grow with anything. No error. "Baseline" kept over "initial" because
   R0 is a recurring per-round endowment, now defined at first use; Aydin may flip the word.
7. 2026-08-14 second pass per Aydin: output-unit agnosticism (output generally; track
   record for what has been produced; agnostic sentence added; table and track-record
   strategy reworded); base-values paragraph replaced by the sweep paragraph (ranges +
   rationales; Lotka-type tail evidence carries [check] pending citation verification) and
   a figure-defaults paragraph. Terminology ripples applied to abstract, intro, and \S2.
8. Randomness paragraph rewritten after Aydin's query; verified in model.R: growth update
   uses the intensity (line 90), never Poisson draws; scoring sums intensities
   (total_expected). The "enough separate acts of research for luck to average out" rationale
   sentence is Claude's proposed gloss of the choice, for Aydin's check.
9. Forward-only check done: every referring phrase has an antecedent; "the same bottleneck
   that limits output limits growth" refers to Lambda, introduced two paragraphs above.
