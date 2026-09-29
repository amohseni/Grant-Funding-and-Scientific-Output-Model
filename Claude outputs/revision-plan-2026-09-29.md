# Revision plan in response to the referee-style review (2026-09-29)

Verdict on the review: the technical core is right and should be acted on. Four of its
claims I verified in-session before writing this plan:

- The two-agent counterexample reproduces exactly: K = (10, 1), B = 2; with R = (10, 1)
  the optimal allocation beats uniform by 0.054, with R = (1, 10) by 0.633. The K-R
  association matters, and our sweeps hold it at rho = 0.
- The second derivative in Lemma 1 reads -4A K^2/S^3 and should read -2A K^2/S^3. This
  is my error from this morning's removal of the factor 2 (the old -4A went with 2A);
  the referee saw v6, which had it. Fixed in the first work package below.
- Corollary 2 proves nonnegativity and vanishing as B grows, not monotone decrease. In
  300 random Pareto populations the ratio plotted in Fig. S2 (targeting's value as a
  share of uniform funding's gain) decreased monotonically in every one, on b from 0.02
  to 3; the absolute difference did not (it starts at zero). So the monotone claim is a
  numerical finding about the ratio, with a stated domain, and the theorem is the limit.
- The alpha_K sweep at fixed Pareto scale moves E[K] from 4.33 to 2 to 1.40. The
  confound is real. The code has the scale parameter (k_min) and the correlation
  parameter (rho_kr), so the normalized and correlated reruns are feasible here.

Two facts about the code that the paper must state and now does not: the review signal
is S = K + N(0, tau_K^2), drawn afresh each round it is used, and the resource signal
R_0 + N(0, tau_R^2), observed once; output is Poisson with the stated rate; posteriors
are importance-sampled with M = 200 particles; the forward-looking funder is a
certainty-equivalent receding-horizon planner (it re-plans the remaining budget each
round conditioning on expected output, and does not value the information its own
grants generate), not a Bayes-adaptive dynamic program.

## Disposition of the ten major points

| # | Referee point | Verdict | Action | Feasible here |
|---|---|---|---|---|
| 1 | "Two features" is too strong; rho and alpha_R matter | Right | Run the review-value and targeting-value results over rho in {-0.5, 0, 0.5, 0.8} and alpha_R in {1.3, 2, 3.5}; test the referee's hypothesis that the deeper quantity is the dispersion of the gaps (cK - R) or of R/K; rewrite the thesis to whatever the data support. Best case: the second feature becomes "how unequally the capability-resource gap is distributed", which is the paper's own idea; capability inequality is then one driver of it. Worst case: "two features, holding the K-R association fixed" | Yes (minutes per grid) |
| 2 | alpha_K sweep confounds inequality with mean capability; fixed additive noise is not fixed informativeness | Right | Rerun Fig. 2C with E[K] held at 2 across alpha (k_min = 2(alpha-1)/alpha); report review informativeness as the rank correlation between signal and capability and rerun at matched correlation; add an n sweep (25, 50, 100, 200) | Yes |
| 3 | Monotone decrease of targeting's value unproved | Right | Main text and SI: "vanishes as the budget grows (Corollary 2) and, in every population we sampled, decreases in it"; define the plotted quantity (ratio) where the claim is made; SI Text 3 no longer attributes monotonicity to the corollary | Yes (text) |
| 4 | phi'' coefficient | Right (my error today) | -2A | Yes |
| 5 | Review is stipulated to observe K; likelihoods unstated; want S = aK + bR + eta and persistent error | Half right | State every likelihood and the planner exactly (must). Add one robustness run, review contaminated by resources, S = K + beta R + eta with the funder knowing beta, for beta in {0, 0.25, 0.5}: shows how much of review's value survives when review partly sees resources. Persistent reviewer error: state as a limitation, not run (it needs a new state variable in the filter) | Likelihoods yes; contaminated review needs a small change to the signal draw and likelihood in model.R, feasible; persistent error no |
| 6 | Lottery theorem narrower than the policy language | Right | Make the divisibility and diminishing-returns conditions explicit at every use; add a three-line example to SI Text 2 (a minimum viable grant size under which a lottery over full awards beats equal division); add indivisibility to the limitations | Yes |
| 7 | Researcher, not researcher-project pair | Right about scope | One sentence in the model and one in limitations: this is a model of investigator-level funding; project-specific quality is the natural extension | Yes |
| 8 | Demote timing; define the planner | Partly | Keep timing in the main text (it is one of the paper's four questions) but move it after spreading, as the extension it is; define the planner (certainty-equivalent receding horizon); "observation, then money" becomes conditional on complementarity | Yes |
| 9 | Empirical bridge; alpha_K justified by productivity tails contradicts the thesis | Right | SI table mapping each parameter to an observable and a plausible range, with honest blanks; reword the tail justification ("productivity tails, which by our own argument understate capability inequality where grants are small; our range therefore extends beyond them"); no calibration in this revision | Yes |
| 10 | Generalize the gap rule to constant-returns technologies | Right, and valuable | Proposition 1 for F = A K h(R/K) with h increasing and strictly concave: equal marginal value gives a common R/K among the funded, so g* = max(cK - R, 0) with c = (h')^{-1}(nu/A). Harmonic and every CES member are instances. Corollary 2 (vanishing) needs h bounded, which holds for the harmonic mean and Leontief and fails for Cobb-Douglas: a real finding, and it explains why the concentration results differ at the Cobb-Douglas end | Yes (proof and text; verify with theorem-proving pass) |

Other points, all accepted: engage Sorzano and Pueche-Granados in one Discussion sentence (their heavy tail is in observed impact; ours decomposes output into latent capability and resources); "Output gives only a lower bound" becomes a statement about expected output; error bars on Fig. 2C (rerun with 200 populations); Fig. 2B gains an output-shortfall measure beside the correlation (needs a T = 20 rerun that records output by round); rerun the production-function robustness with the current allocator so the "earlier routine" sentence goes; the modal edits ("Within the model, no funding mechanism is best in every field"; "in the model, the best use of the budget is..."); mechanical items (grant number, duplicate period, bibliography check notes) are yours.

Declined: moving timing to the SI (your call, I lean keep); adding project quality Q_it as a model component; persistent reviewer error; an empirical calibration. Each is a separate paper's worth of work or a modeling decision, and none blocks the current claims once they are stated at the right strength.

## Work packages, in order

A. Corrections that need no new runs (one session). Lemma 1; Corollary 2 wording in main text, Fig. S2 caption, SI Text 3; expected-output lower bound; exact likelihoods, particle count, and planner definition in Methods and SI Materials and Methods; lottery conditions and the minimum-scale example; investigator-level scope; Sorzano sentence; modal edits; timing subsection moved after spreading, intro's four questions reordered to match. Long version updated in parallel.

B. Generalized gap rule (one session, theorem-proving pass). New Proposition 1 with the CRS family; harmonic as the instance used in the simulations; Corollary 2 with the boundedness condition and the Cobb-Douglas remark; Examples unchanged. Main text states the general rule and says which technologies it covers.

C. Sensitivity runs (compute here; one to two sessions). C1 rho and alpha_R grids for review value and targeting value. C2 mean-normalized alpha_K sweep and matched-informativeness sweep for Fig. 2C, with standard errors, plus n. C3 gap-dispersion analysis: for each cell, the dispersion of cK - R at the cell's budget against the value of targeting and of review, to test the referee's hypothesis. C4 review contaminated by resources. C5 T = 20 rerun recording output shortfall by round (Fig. 2B). C6 production-function robustness rerun on the current allocator (background; the longest). Deliverable: a ledger of every number that changes, then the text.

D. Rewrite of the thesis and Discussion in light of C (one session). The "two features" sentence in Significance, abstract, intro, Discussion, conclusion rewritten to what C supports; the Discussion recap shortened, its space given to what is observable, what needs concavity and divisibility, and what would falsify the model; the SI observables table.

Order: A and B first (they are right regardless of the runs), C next, D last. I will not touch the thesis wording until C is in.
