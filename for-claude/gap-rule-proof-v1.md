# The gap rule: proposition and proof, draft v1 (2026-08-14)

Standalone derivation of Proposition 1, destined for the appendix. Written to be
self-contained: the setup states every assumption the proof uses, and the proof invokes
nothing beyond elementary calculus and convexity. Nothing is assumed true in advance; the
gap rule is derived, not verified. A numerical check against a direct optimizer is
recorded in the notes.

---

## A. The within-round allocation problem

### A.1 Primitives

A population of $n$ researchers. Researcher $i$ has epistemic capability $K_i$ and
baseline resources $R_{i0}$, both fixed within the round. A grant $g_i \geq 0$ adds to
researcher $i$'s resources for the round, so that their expected output is
$$\lambda_i(g_i) = A\,\Lambda(K_i,\, R_{i0} + g_i), \qquad
\Lambda(K, R) = \frac{2KR}{K + R},$$
with $A > 0$ the common productivity scale and $\Lambda$ the harmonic production function
of \S3. Write
$$\varphi_i(g) = \frac{2A\,K_i\,(R_{i0} + g)}{K_i + R_{i0} + g}$$
for researcher $i$'s expected output as a function of their grant alone.

### A.2 Assumptions

Throughout: $A > 0$; $K_i > 0$ and $R_{i0} \geq 0$ for every $i$; and the round's tranche
$B > 0$. (Under the model's continuous distributions for $K$ and $R$, the conditions on
$K_i$ and $R_{i0}$ hold with probability one.) Nothing else about the model is used: the
proof does not depend on the distributions of capability and resources, on the noise in
realized output, or on the dynamics across rounds. In particular, the result holds for
every realized population and every tranche, not merely on average.

### A.3 The problem

The funder observes $(K_i, R_{i0})_{i=1}^n$ and divides the tranche so as to maximize
expected output:
$$(\mathrm{P}) \qquad \max_{g_1, \ldots, g_n \geq 0} \; \sum_{i=1}^n \varphi_i(g_i)
\quad \text{subject to} \quad \sum_{i=1}^n g_i \leq B.$$

## B. Three lemmas

**Lemma 1** (the return to a grant). *For each $i$, on $g \in [0, \infty)$:*

*(i) $\varphi_i$ is strictly increasing, with*
$$\varphi_i'(g) = \frac{2A\,K_i^2}{(K_i + R_{i0} + g)^2} > 0;$$

*(ii) $\varphi_i$ is strictly concave:*
$$\varphi_i''(g) = -\,\frac{4A\,K_i^2}{(K_i + R_{i0} + g)^3} < 0;$$

*(iii) $\varphi_i'(g) < 2A$ whenever $R_{i0} + g > 0$, and $\varphi_i(g) \uparrow 2A K_i$
as $g \to \infty$.*

*Proof.* Write $S = K_i + R_{i0} + g$. Then $\varphi_i(g) = 2AK_i(S - K_i)/S =
2AK_i - 2AK_i^2/S$, and (i) and (ii) follow by differentiating in $g$ (note $dS/dg = 1$).
For (iii): $\varphi_i'(g) = 2A\,K_i^2/S^2 < 2A$ exactly when $S > K_i$, that is, when
$R_{i0} + g > 0$; and $\varphi_i(g) = 2AK_i - 2AK_i^2/S \to 2AK_i$ as $S \to \infty$.
$\square$

Lemma 1 states the economics of the problem. The marginal value of a dollar to researcher
$i$, $\varphi_i'(g)$, increases with their capability and decreases with the resources
they already hold, grants included. Each dollar granted lowers the value of the next
(strict concavity). And expected output increases with every dollar but never exceeds
$2AK_i$: it is bounded above by a quantity set by capability, which it approaches as
$g \to \infty$ but does not attain.

**Lemma 2** (existence, uniqueness, exhaustion). *Problem $(\mathrm{P})$ has exactly one
solution $g^* = (g_1^*, \ldots, g_n^*)$, and it exhausts the tranche:
$\sum_i g_i^* = B$.*

*Proof.* The feasible set $F = \{g \geq 0 : \sum_i g_i \leq B\}$ is compact and convex,
and the objective $f(g) = \sum_i \varphi_i(g_i)$ is continuous, so a maximizer exists.
$f$ is strictly concave on $F$: for distinct feasible $g \neq h$ and $t \in (0,1)$,
concavity of each $\varphi_i$ gives $f(tg + (1-t)h) \geq t f(g) + (1-t) f(h)$, and the
inequality is strict because $g_j \neq h_j$ for some $j$ and $\varphi_j$ is strictly
concave. If two distinct maximizers existed, their midpoint would be feasible (convexity
of $F$) and strictly better, a contradiction; so the maximizer is unique. Finally, if
$\sum_i g_i^* < B$, then increasing any one $g_i^*$ slightly stays feasible and raises
$f$ (Lemma 1(i)), contradicting optimality; so $\sum_i g_i^* = B$. $\square$

**Lemma 3** (equal marginal value). *A feasible allocation $g$ with $\sum_i g_i = B$
solves $(\mathrm{P})$ if and only if there is a $\nu > 0$ such that*
$$\varphi_i'(g_i) = \nu \;\text{ for every } i \text{ with } g_i > 0,
\qquad \varphi_i'(0) \leq \nu \;\text{ for every } i \text{ with } g_i = 0.$$

*Proof.* **Necessity.** Let $g^*$ solve $(\mathrm{P})$. Since $B > 0$ is exhausted
(Lemma 2), some researcher is funded: $g_i^* > 0$ for some $i$. Fix any such $i$ and any
$j \neq i$, and for $0 < \delta \leq g_i^*$ consider the transfer that replaces $g_i^*$
by $g_i^* - \delta$ and $g_j^*$ by $g_j^* + \delta$, leaving the rest unchanged. The
result is feasible, so optimality gives
$$\varphi_j(g_j^* + \delta) - \varphi_j(g_j^*) \;\leq\; \varphi_i(g_i^*) -
\varphi_i(g_i^* - \delta).$$
Dividing by $\delta$ and letting $\delta \downarrow 0$ yields $\varphi_j'(g_j^*) \leq
\varphi_i'(g_i^*)$. If $j$ is also funded, the same argument with $i$ and $j$ exchanged
gives the reverse inequality, so all funded researchers share a common marginal value;
call it $\nu = \varphi_i'(g_i^*) > 0$ (positivity by Lemma 1(i)). For unfunded $j$, the
displayed inequality reads $\varphi_j'(0) \leq \nu$.

**Sufficiency.** Let $g$ satisfy the condition with multiplier $\nu$, and let $h$ be any
feasible allocation. Concavity of $\varphi_i$ gives, for every $i$,
$$\varphi_i(h_i) \;\leq\; \varphi_i(g_i) + \varphi_i'(g_i)\,(h_i - g_i).$$
For funded $i$, $\varphi_i'(g_i) = \nu$, so $\varphi_i'(g_i)(h_i - g_i) =
\nu\,(h_i - g_i)$. For unfunded $i$, $g_i = 0$ and $h_i - g_i = h_i \geq 0$, so
$\varphi_i'(0)\,(h_i - g_i) \leq \nu\,(h_i - g_i)$. Summing over $i$:
$$f(h) \;\leq\; f(g) + \nu \sum_i (h_i - g_i) \;=\; f(g) + \nu\Big(\sum_i h_i - B\Big)
\;\leq\; f(g),$$
since $h$ is feasible and $\nu > 0$. So $g$ solves $(\mathrm{P})$. $\square$

Lemma 3 is the formal content of the sentence in \S4: an allocation is optimal exactly
when no dollar can be moved to a researcher who can produce more value with it.

## C. The proposition

**Proposition 1** (the gap rule). *The unique solution of $(\mathrm{P})$ is*
$$g_i^* = \max(c\,K_i - R_{i0},\, 0), \qquad i = 1, \ldots, n,$$
*where $c > 0$ is the unique constant satisfying the budget identity*
$$\sum_{i=1}^n \max(c\,K_i - R_{i0},\, 0) = B.$$
*Equivalently, $c = \sqrt{2A/\nu} - 1$, where $\nu \in (0, 2A)$ is the common marginal
value of Lemma 3.*

*Proof.* Let $g^*$ be the unique solution (Lemma 2) and $\nu$ the multiplier that Lemma 3
associates with it. The proof proceeds in four steps.

**Step 1: $\nu \in (0, 2A)$, so $c = \sqrt{2A/\nu} - 1$ is well defined and positive.**
Positivity of $\nu$ is part of Lemma 3. For the upper bound: some researcher $i$ is
funded, and for them $R_{i0} + g_i^* \geq g_i^* > 0$, so Lemma 1(iii) gives
$\nu = \varphi_i'(g_i^*) < 2A$. Then $2A/\nu > 1$, so $c = \sqrt{2A/\nu} - 1 > 0$.

**Step 2: funded researchers receive exactly their gap.** Let $g_i^* > 0$. By Lemma 3,
$$\frac{2A\,K_i^2}{(K_i + R_{i0} + g_i^*)^2} = \nu.$$
Solving, and taking the positive square root (both sides of the next display are
positive),
$$K_i + R_{i0} + g_i^* = K_i\,\sqrt{2A/\nu} = (1 + c)\,K_i,$$
so $g_i^* = c\,K_i - R_{i0}$. Since $g_i^* > 0$, this equals
$\max(c\,K_i - R_{i0},\, 0)$.

**Step 3: unfunded researchers have no gap.** Let $g_i^* = 0$. By Lemma 3,
$$\frac{2A\,K_i^2}{(K_i + R_{i0})^2} = \varphi_i'(0) \leq \nu,$$
so $K_i + R_{i0} \geq K_i\sqrt{2A/\nu} = (1 + c)\,K_i$, that is, $R_{i0} \geq c\,K_i$.
Then $\max(c\,K_i - R_{i0},\, 0) = 0 = g_i^*$.

Steps 2 and 3 establish $g_i^* = \max(c\,K_i - R_{i0}, 0)$ for every $i$, and the budget
identity is exhaustion (Lemma 2). It remains to show $c$ is the only constant satisfying
that identity.

**Step 4: the budget identity pins $c$ uniquely.** Define
$$G(c) = \sum_{i=1}^n \max(c\,K_i - R_{i0},\, 0), \qquad c \geq 0.$$
$G$ is continuous and nondecreasing, with $G(0) = 0$ and $G(c) \to \infty$ as
$c \to \infty$. Moreover $G$ is strictly increasing wherever it is positive: if
$G(c) > 0$, then $c\,K_j > R_{j0}$ for some $j$, and for any $c' > c$,
$G(c') \geq G(c) + K_j\,(c' - c) > G(c)$. Now suppose $c_1 < c_2$ both satisfy
$G(c) = B$. Since $B > 0$, $G(c_1) > 0$, so $G(c_2) > G(c_1) = B$, a contradiction.
$\square$

## D. Consequences

**Corollary 1** (funding frontier and budget comparative statics). *(i) Researcher $i$
is funded, $g_i^* > 0$, exactly when $K_i > R_{i0}/c$. (ii) As a function of the tranche
$B$, the constant $c(B)$ is continuous and strictly increasing; hence the funded set is
nondecreasing in $B$, every funded researcher's grant $c(B)\,K_i - R_{i0}$ is strictly
increasing in $B$, and the common marginal value $\nu(B) = 2A/(1 + c(B))^2$ is strictly
decreasing in $B$.*

*Proof.* (i) is immediate from the formula: $\max(cK_i - R_{i0}, 0) > 0$ iff
$cK_i > R_{i0}$. For (ii): $G$ is continuous, strictly increasing where positive, and
unbounded, so $c(B) = G^{-1}(B)$ is well defined, continuous, and strictly increasing on
$B > 0$ (Step 4). The remaining claims follow by monotonicity of $c \mapsto cK_i -
R_{i0}$ and $c \mapsto 2A/(1+c)^2$. $\square$

**Corollary 2** (the value of targeting vanishes in the budget). *Let $f^*(B)$ denote the
optimal value of $(\mathrm{P})$ and $f^u(B) = \sum_i \varphi_i(B/n)$ the value of uniform
funding. Then $0 \leq f^*(B) - f^u(B)$ for every $B$, and $f^*(B) - f^u(B) \to 0$ as
$B \to \infty$.*

*Proof.* The first inequality holds because uniform funding is feasible. For the second:
by Lemma 1(iii), $\varphi_i(g) < 2AK_i$ for every $g$, so $f^*(B) \leq 2A\sum_i K_i$; and
$f^u(B) = \sum_i \varphi_i(B/n) \to 2A\sum_i K_i$ as $B \to \infty$. Hence
$0 \leq f^*(B) - f^u(B) \leq 2A\sum_i K_i - f^u(B) \to 0$. $\square$

As the budget increases, any allocation that spreads funds widely brings every
researcher's output near its upper bound $2AK_i$, so allocations come to differ little in
total output. Since uniform funding's own gain over no funding does not vanish (it
approaches $2A\sum_i K_i - \sum_i \varphi_i(0) > 0$), the same conclusion holds for the
value of targeting measured relative to that gain.

**Example 1** (funding the track record can produce less than uniform funding). Two
researchers, $A = 1$, tranche $B = 2$. Researcher 1: $K_1 = 8$, $R_{10} = 8$ (productive
and well resourced; expected output $\varphi_1(0) = 8$). Researcher 2: $K_2 = 4$,
$R_{20} = 1$ (capable but resource-starved; $\varphi_2(0) = 8/5$). Allocating in
proportion to expected output gives $g = (5/3,\, 1/3)$ and total output
$570/53 \approx 10.75$. Uniform funding, $g = (1, 1)$, gives $568/51 \approx 11.14$. The
gap rule ($c = 3/4$) gives $g^* = (0, 2)$ and $80/7 \approx 11.43$. Track-record funding
sends most of the tranche to the researcher with the small gap and produces less output
than uniform funding.

**Example 2** (funding the under-resourced can produce less than uniform funding). Two
researchers, $A = 1$, tranche $B = 2$. Researcher 1: $K_1 = 10$, $R_{10} = 5$.
Researcher 2: $K_2 = 1/2$, $R_{20} = 0$ (the scarcest resources in the population).
Directing the tranche to the most resource-starved researcher gives total output
$112/15 \approx 7.47$. Uniform funding gives $49/6 \approx 8.17$. The gap rule
($c = 2/3$) gives $g^* = (5/3,\, 1/3)$ and $42/5 = 8.4$. Scarcity alone misdirects the
tranche toward low capability; note that the optimum still grants the scarce researcher
a little, because their marginal value at zero is high, but that marginal value is
determined by capability, not by scarcity.

**Remark 1** (water-filling). At the optimum, every funded researcher's total resources
are filled to the capability-scaled level $R_{i0} + g_i^* = c\,K_i$, equivalently
$K_i + R_{i0} + g_i^* = (1+c)K_i$. Funding to a common marginal value by filling
shortfalls from a level is the water-filling solution of information theory (Cover and
Thomas, 2006); here the fill level scales with capability rather than being common to
all.

**Remark 2** (the capability-independent special case). If all researchers share one
capability, $K_i = K$, the fill level $cK$ is common, and the rule reduces to filling
each researcher's resource shortfall from a single target: the poverty-gap transfer of
Foster, Greer, and Thorbecke (1984).

**Remark 3** (scope across rounds). Proposition 1 governs the division of a given tranche
within a round, taking that round's $(K_i, R_{i0})$ as fixed. Across rounds, capability
compounds and the constant $c$ and the funded set are recomputed from the new population
state; nothing in the proof depends on how the population state arose. When the funder
also chooses how to split the budget across rounds (\S7), the proposition governs the
within-round division of whatever tranche that choice assigns.

---

## Notes for Aydin

1. The proof is self-contained: no KKT machinery is invoked. Lemma 3 does the work KKT
   would do, with necessity by the transfer argument (which is also the sentence in \S4
   made formal) and sufficiency by the concavity gradient inequality. I judged this
   cleaner and more checkable than citing Karush-Kuhn-Tucker for a problem this simple;
   say the word if you prefer the standard citation instead.
2. The primary characterization of $c$ is the budget identity, which involves only
   observables of the problem $(K_i, R_{i0}, B)$; the $\nu$-form
   $c = \sqrt{2A/\nu} - 1$ is derived, matching the main text's footnote.
3. Two boundary details the meticulous setup surfaced: (a) $\nu < 2A$ is needed for
   $c > 0$ and is not automatic; Step 1 proves it from the existence of a funded
   researcher. (b) A researcher with $R_{i0} = 0$ is funded at every budget (their
   $\varphi_i'(0) = 2A$ exceeds every admissible $\nu$), which is why the frontier
   statement is $K_i > R_{i0}/c$ rather than a condition that could exclude them.
4. Corollary 1(ii) is the formal basis for \S4's sentence "a larger budget raises $c$,
   moving the frontier down and funding more researchers more deeply," plus the
   marginal-value decline used in \S8-\S9.
5. Numerical verification (this session, verify_gap_rule.py): 81 trials, n = 50,
   crossing tail parameters 1.3/2.0/3.5, budget scales 0.1/0.5/1.0, and A = 0.5/1/2, with
   the R_{i0} = 0 boundary exercised. The gap-rule allocation (bisection on the budget
   identity) against scipy's SLSQP direct optimizer: the formula's objective was never
   below the optimizer's in any trial (max shortfall 0.0; where they differed at all, the
   formula won); max allocation deviation 6.5e-7 of the budget (optimizer tolerance); the
   budget identity held to 6e-16; funded researchers' marginal values matched
   nu = 2A/(1+c)^2 to 9e-16 and every unfunded researcher's marginal value was <= nu.
   Funded counts ranged 2 to 50 of 50, so the frontier was active. Corollary 1(ii)
   spot-checked: c(B) strictly increasing across b = 0.05 to 2 on a fixed population.
6. Conversion to LaTeX (amsart theorem environments) once you approve the content.
