# Section 4 (The optimal allocation), draft v1 (2026-08-14)

Drafted plain-first from the function analysis (chat, 2026-08-14). One object (the marginal
value of a dollar), one result (the gap rule, Proposition 1, proof to appendix), three
consequences (two refutations; the budget comparative static). Precedents at the
derivation per the section-4 ledger. Figures are placeholders pending regeneration.

---

## 4. The optimal allocation: fund the capability-resource gap

First, consider the simple case where the funder has complete information regarding
researchers' capabilities and resources. Further, we initially restrict our attention to a
single round of funding, temporarily putting aside the contribution of resources to the
accumulation of researchers' capabilities across rounds. The funder's problem is to divide
the round's tranche among researchers so as to maximize expected output:
$$\max_{g_1, \ldots, g_n \geq 0} \; \sum_i \lambda_i(K_i, R_{i0} + g_i)
\quad \text{subject to} \quad \sum_i g_i \leq \text{the round's tranche}.$$

Now, consider the marginal value of a dollar for such a funder. Under the harmonic
production function,
$$\frac{\partial \lambda_i}{\partial g_i} = \frac{2A\,K_i^2}{(K_i + R_{i0} + g_i)^2}.$$
The marginal value of funding researcher $i$ increases with their capability and decreases
with the resources they already hold, including grants. An allocation is optimal exactly
when no dollar can be moved to a researcher who can produce more value with it, that is,
when marginal values are equal across all funded researchers and no unfunded researcher's
marginal value exceeds theirs.

**Proposition 1** (the gap rule). *For every researcher, the optimal grant fills the gap,
if any, between the resources they hold and a target proportional to their capability:*
$$g_i^* = \max(c\,K_i - R_{i0},\, 0),$$
*where the constant $c > 0$ is set by spending the budget.* [Proof: Appendix.
Footnote: $c = \sqrt{2A/\nu} - 1$, with $\nu$ the common marginal value at the optimum.]

We can understand the rule as follows. Each researcher's target is $c$ times their
capability, and their gap is the difference between their target and the resources they
hold. The rule grants every researcher with a positive gap exactly their gap, and grants
nothing to the rest; the constant $c$ is defined as the value at which these grants sum to
the budget. Because $c$ is defined by the budget, the budget always covers the gaps: a
question of which gaps to fund first does not arise. When the budget is small, $c$ is
small: targets are low, few researchers have positive gaps, and their grants are small.
When the budget is large, $c$ is large: targets are high, and funding extends across
nearly the whole population. In both cases the funded researchers, those with positive
gaps ($R_{i0} < c\,K_i$), are those with the least resources relative to their capability,
and these are the researchers for whom a dollar has the greatest marginal
value.\footnote{The optimality condition behind the rule, funding to a common marginal
value, is the water-filling solution familiar from information theory (Cover and Thomas,
2006). If all researchers shared one capability, the rule would reduce to filling each
researcher's gap below one common target, the poverty-gap transfer of development
economics (Foster, Greer, and Thorbecke, 1984); our rule scales the target with
capability.} The condition $R_{i0} < c\,K_i$ draws a frontier through the population, the
funding frontier, separating the funded from the unfunded.

[FIGURE 1 placeholder: optimal grants across a population in (K, R) space, with the funding
frontier; report Fig. 1, regenerate at 200 seeds per provenance rule.]

The gap rule explains why each of the intuitive funding schemes discussed earlier fails.
Funding the track record fails because output reflects resources as well as capability: a
productive researcher whose resources approach or exceed their target has a small or
negative gap, and the marginal dollar buys little there. The failure is not merely a shortfall from the optimum: because
track-record funding channels money toward researchers whose bottleneck is not resources,
it can produce less output than spreading the budget evenly, that is, than not targeting
at all.\footnote{The appendix gives a two-researcher example.} Funding the under-resourced
fails for the mirror-image reason: scarcity is no evidence of capability, and dollars
allocated to researchers of low capability buy little however scarce their resources. It,
too, can produce less than not targeting at all.\footnote{Again, the appendix gives an
example.} Each scheme accounts for one of capability and resources, whereas the marginal
value of a dollar is determined by how a researcher's resources compare with their
capability, and so responds to both.

How much targeting matters depends on the budget. A researcher's expected output
increases with every dollar of resources, but with diminishing marginal returns, and it
never exceeds $2AK_i$: expected output is bounded above by a quantity set by capability,
which it approaches but does not attain. When the budget is large, optimal and uniform
funding therefore differ little in the output they produce: under either, every
researcher's output is near its upper bound, and the difference between the two vanishes
as the budget increases (appendix). When the budget is small, dollars are scarce,
researchers differ in the marginal value of a dollar, and output depends on where each
dollar goes. Targeting matters where money is tight and little where it is ample. This
comparative static recurs: it sets the price of overriding the optimal allocation with
floors and lotteries (\S8), and it determines when optimal funding concentrates rather
than spreads (\S9).

To use the gap rule, however, a funder requires an estimate of the gap between
researchers' capabilities and their resources. We turn now to how such estimation can
work: from the track record alone (\S5), and with peer review added (\S6).

---

## Notes for Aydin

0. 2026-09-05 figurative-language purge: every sentence in the section now passes the
   substitution test against the formalism. Verify in particular: (a) the bound sentence
   in the budget paragraph (increases with every dollar, diminishing returns, bounded
   above by 2AK_i, approached not attained); (b) "at capacity" replaced with "whose
   resources approach or exceed their target" because your ceiling correction applies to
   it too (no researcher is ever strictly at capacity); veto if you want your wording
   back; (c) Proposition 1 now "fills the gap, if any" (one concept one term; "shortfall"
   dropped; the proof doc still says shortfall in two remarks, to align at LaTeX
   conversion); (d) the abstract's "rises" became "increases" under your "standardly
   throughout" (the abstract was locked; veto if it should not have been touched);
   (e) \S3's "Capability grows through research" is your approved text and was left;
   change to "increases" on your word.

1. NUMBERS OUT (2026-08-14, your call): the section is now qualitative only; every claim
   in it is derivable from the gap rule. The dropped simulation numbers (track-record
   23.5 vs uniform 26.1 at defaults; 30.1 vs 29.8 at heavy tails; targeting value 12
   percent at b = 0.1 to 5 percent at b = 1) are banked in section-4-notes.md for the
   quantitative sections, where the parameter context exists to make them meaningful.
2. Both "can produce less than not targeting at all" claims are now ANALYTIC: exact
   two-researcher examples added to the proof document (Examples 1 and 2, verified in
   exact rational arithmetic), one per option. The budget static is also now analytic:
   Corollary 2 (proof doc) shows the optimal-minus-uniform output difference vanishes as
   the budget grows, absolutely and relative to uniform's own gain. The old "optimal
   allocation converges toward uniform funding" sentence was WRONG as stated (grant
   shares converge toward capability-proportional, not equal) and is gone; the true
   claim, outputs converge, is what Corollary 2 proves.
3. Your closing sentence applied with one terminology alignment flagged for your veto:
   "the gap between researchers' funds and their capabilities" became "the gap between
   researchers' capabilities and their resources" (locked term is resources; order
   matches "capability-resource gap"). Second-sentence candidates, ranked:
   (a) "We turn now to how such estimation can work: from the track record alone (\S5),
   and with peer review added (\S6)." [in the draft; earns the roadmap]
   (b) "We turn now to how such estimation might work." [your original]
   (c) "We turn now to what the funder can learn about the gap." Your word choice.
4. 2026-09-05 flags from your line-edit pass:
   (a) "schemes" adopted and PROPAGATED to the intro in main.tex (three spots: "two
   intuitive schemes," "Each scheme can be intuitively compelling," "Each intuitive
   scheme fails"); one concept, one term. State.md terminology lock updated.
   (b) Your "who's at capacity": applied as "who is at capacity." Two flags: "capacity"
   sits one letter from "capability" and is not otherwise used in the paper (an earlier
   consistency pass kept "capabilities" over your one "capacities"); if you want zero
   collision risk, an alternative is "a productive researcher whose resources already
   meet their target has a small or negative gap." Your call; his wording stands in the
   draft.
   (c) Your closing sentence "determined by their difference" is not applied as written,
   flagged for truth: the marginal value 2AK^2/(K+R)^2 is determined by resources
   relative to capability (the ratio R/K), not by the difference K - R; two researchers
   with the same difference at different scales have different marginal values. (The
   GRANT is a scaled difference, cK - R; the marginal value is not.) In the draft:
   "...whereas the marginal value of a dollar is determined by how a researcher's
   resources compare with their capability, and so responds to both." Ranked
   alternatives: (i) as drafted; (ii) "...is determined by the two together";
   (iii) your verbatim, if you judge "difference" acceptable as a loose gloss.
5. FIGURE 2 (targeting value vs budget) removed from this section with the numbers;
   banked for \S8. FIGURE 1 stays but is now an analytic illustration: optimal grants in
   (K, R) space computed directly from Proposition 1 on a sample population, no
   simulation needed, so it needs no seed provenance.
6. "Scarcity is no evidence of capability" holds at rho = 0 and understates the case at
   positive rho (where scarcity is evidence of LOW capability). The stronger conditional
   is available; the weak form is never false in the direction that would rescue the
   option.
7. Forward-only check redone after the rewrite. "The round's tranche" is defined in S3;
   for forward strategies S7's re-planning applies and the proposition governs the
   within-round division either way (proof doc, Remark 3).
