# Sweep handoff addendum: Packages D, E, F

*2026-08-05, follows docs/SWEEP_HANDOFF_2026-08-05.md (Packages A-C). Revised same day: Package F
(grant persistence) REMOVED by Aydin's decision, and E-1 redesigned to include the product-reduction
test; see the notes in those sections. Same global conventions (section 0 there): commit-first,
validation gates after any model.R edit, 200 seeds unless stated, paired contrasts, percent-of-S1
normalization, RUN_INFO.txt per output dir, results appended to docs/DIAGNOSTICS_RESULTS_2026-08.md
with CONFIRMED / REFUTED / MIXED per prediction. Run after A-C complete.*

*Also confirm the bootstrap-suite verification items are queued (they may already be): the honest
fixed-schedule sweep at the exploration-depth corner; the CE self-consistency check there; 64 to 200
seed re-runs of any resource_regime / exploration cells destined for main-text figures. Spec:
for-claude drafts/grant-funding/notes-bootstrap-integration.md section 4; summarized at the end of
this file in case that folder is unavailable.*

---

## Package D: information integrity (the paper already promises these)

### D-1. Model-implied AUC: locating real peer review on the precision axis (standalone script, no model runs)

**Why.** The paper's empirical anchor is Fang et al.'s AUC of 0.54 for NIH percentile scores, read
as evidence that real review is a very noisy signal sitting above the precision elbow. That reading
is currently qualitative. Make it quantitative: compute the model's implied AUC as a function of
signal noise, then invert.

**How.** Pure Monte Carlo, no simulation engine: draw K_i ~ Pareto(k_min=1, shape=alpha_K),
sigma_i = K_i + N(0, tau_k^2); compute the AUC of sigma ranking for identifying top-q researchers
by K. Grid: alpha_K in {1.3, 2, 3.5} x tau_k in {0.05, 0.1, 0.3, 0.5, 1, 1.5, 2.5, 5, 10, 20} x
q in {0.1, 0.2}; 10^6 draws per cell (seconds). Output: AUC(tau_k) curves per (alpha_K, q), and the
implied tau_k solving AUC = 0.54 per curve, with a bracketing interval.

**Prediction P-D1:** the implied tau_k at AUC = 0.54 lies well above the signal-value half-decay
point (tau_k ~ 2.5 at heavy tails), so the "NIH sits above the elbow; sharpening review still pays"
sentence survives quantification. If instead the implied tau_k lands below the elbow, footnote 13's
argument reverses and the paper must say so.

Output: `sweep_results/D_auc_calibration/` (CSV + one figure: AUC vs tau_k with the 0.54 line and
the elbow marked).

### D-2. Misspecified trust: the funder's believed tau_k vs the true tau_k (small code change)

**Why.** Every signal result assumes the funder weights the signal at its true precision. The
realistic failure is overtrust of noisy review. The peer-review section needs its caveat curve, and
this is the most referee-attractive missing experiment in the information story.

**How.** Decouple generation from inference: `run_simulation_T` gains `tau_k_true = tau_k` and
`tau_k_belief = tau_k`; `draw_signals` uses the true value, `loglik_grant_signal` (and any other
likelihood touchpoint) uses the belief. Grep every tau_k use and classify generative vs inferential
before editing; validation gates as usual (defaults reproduce bit-identically).

Grid (T=2, smooth, strategies c(1,2,4,5,7,8), 200 seeds): tau_k_true in {0.3, 1, 3} x
tau_k_belief in {0.1, 0.3, 1, 3, 10} x alpha_K in {1.3, 2}. 30 cells. Readouts: S5-S4 and S8-S7
(the value of using the signal, given the belief), as % of S1.

**Prediction P-D2 (asymmetry):** overtrust of a noisy signal (belief 0.1-0.3, truth 3) drives the
signal's value negative (worse than pubs-only); undertrust of a sharp signal (belief 3-10, truth
0.3) merely forfeits part of the value, staying nonnegative. If overtrust never goes negative, the
caveat softens to "misweighting only attenuates," which is itself worth reporting.

Output: `sweep_results/D_misspecified_trust/`.

### D-3. Resource-signal ablation (config-only, the report itself flags this as not yet run)

**Why.** Report Appendix A: "measuring the resource signal's own value requires a sweep over its
presence and absence, which is not yet run." Close the report's own open item.

**How.** No code change: `use_resource_signal` already exists. Grid (T=2, smooth, 200 seeds):
use_resource_signal in {TRUE, FALSE} x rho in {-0.5, 0, 0.8} x alpha_K in {1.3, 2}, strategies
c(1,4,5,7,8). 12 cells. Readout: each strategy's output with vs without sigma_R, and the S8-S7
signal value with vs without.

**Prediction P-D3:** the resource signal's own contribution is small at rho >= 0 (redundancy: pubs
plus the talent signal identify R via 1/R = A'/lambda - 1/K), largest at rho = -0.5 where the
talented are resource-poor and output alone under-identifies. Small everywhere = the redundancy
claim as stated; nontrivial at negative rho = add one sentence of scope.

Output: `sweep_results/D_resource_ablation/`.

### D-4. Allocation convergence to the gap rule (small detail sweep)

**Why.** The thesis welds Part I to Part IV: the realized (posterior-based) allocation approaches
the full-information rule g* = cK - R as information sharpens. Footnote 3's r = 0.95-1.00 is a
point claim at base settings; the paper wants the curve.

**How.** Sweep tau_k in {0.05, 0.3, 1, 3, 10, 20} x alpha_K in {1.3, 2} at T=2, strategy S5 (and
S4 as the pubs-only baseline), 50 seeds with per-researcher grants captured (detail=TRUE or
equivalent). Per cell: corr(g_realized, max(cK - R, 0)) with c computed from the same budget, and
the output gap to the full-information allocator (an oracle run with K, R passed directly, which
`docs/` notes as Tier-4 machinery; if no oracle mode exists, add the trivial one: posteriors
replaced by point masses at truth).

**Prediction P-D4:** corr rises monotonically toward the 0.95-1.00 range as tau_k falls; the
output gap to oracle shrinks correspondingly; the pubs-only baseline plateaus below it. This is
the "cost of not knowing the gap" figure stated in the allocation itself.

Output: `sweep_results/D_gap_convergence/`.

---

## Package E: headline insurance (cheap; directions only)

### E-1. Heterogeneous productivity: the scaled rule, and the product-reduction test (moderate code change)

**Why, two questions.** (i) "Every researcher has the same production function" is the limitations
list's sharpest entry. The closed form survives heterogeneity in A: with
lambda_i = A_i Lambda(K_i, R_i), water-filling gives g*_i = c_i K_i - R_i with
c_i = sqrt(2 A_i / nu) - 1, a productivity-scaled gap rule. (ii) Aydin's reduction hypothesis:
variation in (A, K) collapses into the distribution of the product T = A K (which stays heavy-tailed
under power-law or lognormal factors), so a one-dimensional-talent model matched on the product
reproduces the results. The reduction is NOT exact algebraically: A and K enter through two
different combinations (output ceiling 2AK as funding grows; saturation scale K alone, since
lambda = AK, half the ceiling, at R = K), and at fixed product the optimal grant
sqrt(2T/nu)*sqrt(K) - K - R still varies with the split. Whether the differences are material for
the paper's quantities is an empirical property of the model. Either outcome is useful: if the
reduction approximately holds, homogeneous-A is justified as more than convenience (one sentence in
the setup); if it fails, the appendix documents which quantities separate, and the empirical
discussion gains the point that A and K are separately identified only by dose-response (output vs
funding depth), never by cross-sectional output.

**How, two parts.**

- **E-1a (rule tracking; A observable).** Draw A_i ~ Lognormal with E[A] matched to the homogeneous
  value, sdlog in {0, 0.25, 0.5}; A_i is KNOWN to the funder (an institutional multiplier: it enters
  the likelihood and both planners as data). Thread productivity as scalar-or-vector (R
  broadcasting; audit the lambda_rate call sites). Grid: sdlog x signal_value's four corners
  (alpha_K in {1.3, 3.5} x tau_k in {0.3, 10}), T=2, 100 seeds. Readouts:
  corr(g_realized, max(c_i K_i - R_i, 0)) with c_i per the scaled rule; signs of S8-S7 and S9-S7.
  **Prediction P-E1a:** grants track the scaled rule at r > 0.9; all headline directions survive.

- **E-1b (reduction test; A latent).** Fix a target distribution for log T (matching the baseline
  model's log-talent: all variance in K). Construct populations that split the SAME total variance
  of log T between two independent latent factors: Var(log A) = (1 - s) V, Var(log K) = s V, with
  s in {1.0 (baseline), 0.7, 0.4}, means matched so E[log T] is constant. A_i is latent with known
  prior; the funder's posterior atoms become (K_s, A_s) pairs drawn from the two priors (M
  unchanged); the likelihood evaluates lambda = A * Lambda(K, R); the grant signal remains on K
  only (state this design choice in the memo: the signal observes the knowledge factor, so as s
  falls, review sees a shrinking share of what matters, and that is part of what "reduction" means).
  Grid: s x tail-heaviness of T in {heavy, moderate} x tau_k in {0.3, 1}, T=2, 100 seeds, common
  random numbers in T across s cells. Readouts vs s: total optimal output; Gini of grants; signal
  value S8-S7; rank correlation of grant vectors across s (do the same researchers get funded?).
  **Prediction P-E1b (registered, Claude's):** quantities shift with s beyond noise, with signal
  value falling as s falls; if all shifts stay within ~10-15% of the s=1 baseline, report the
  reduction as approximately valid and the paper keeps one-dimensional talent; larger shifts go to
  the appendix with the separation documented.

### E-2. Headline surface at scale (config-only)

**Why.** The forward-advantage surface (T x epsilon) was run at n=50, and the finite-n null was
checked at T=2 only; the budget-independence claim was checked at T=2 and at low epsilon only.

**How.** Two config sweeps, 100 seeds: (a) horizon_growth at n=200 (25 cells); (b) T in {2,3,5} x
b in {0.1, 0.5, 1} at epsilon = 0.85 (9 cells). Readouts: PG and b_idx_S8.

**Prediction P-E2:** the surface shape and b_idx pattern are unchanged at n=200; at high epsilon
the forward advantage scales roughly with the purse but its sign and T-monotonicity are b-invariant.

---

## Package F: REMOVED (decision by Aydin, 2026-08-05)

Grant persistence will not be modeled or swept. Rationale, recorded so no future session
re-proposes it as a gap: grants do not buy durable capacity outside endowments; they are consumed
(salaries, postdocs, running costs), and what persists after a grant ends is skills, publications,
and reputation, which is exactly what K is. The paid-growth channel epsilon * Lambda(K, R0 + g)
already carries the grant's durable residue. The draft therefore states the consumable-grant
assumption in the setup as a substantive modeling claim, not a concession; the limitations section
gets one line naming the exceptions (endowments; capital-intensive fields where durable equipment
is the fixed input, which is the Leontief end of the CES family and so folds into the
substitutability discussion). No code change, no sweep.

---

## Bootstrap verification items (from notes-bootstrap-integration.md, in case that file is unreachable)

1. Honest fixed-schedule sweep at the depth corner (b=3): grid over two-block schedules (share x of
   the purse in rounds 1..k, remainder after), each executed with the myopic within-round allocator
   under true dynamics, no CE planner. Certifies seed-and-harvest planner-free if no early-mass
   schedule beats even-or-late. Add one poverty-epsilon cell.
2. CE self-consistency at depth: does S8 maximize its own CE objective there (honest mispricing) or
   score below S5 on it (bug)? One cell.
3. Re-run at 200 seeds any resource_regime / exploration cells whose numbers or figures enter the
   main text (64-seed provenance otherwise noted in captions).

## Priority order if time-limited

D-1 and D-3 first (hours, no or trivial code change, both promised by the current text). Then D-2
and D-4 (the caveat curve and the thesis figure). Then E-2 (config-only), E-1a, E-1b (the largest
code change in this addendum; the posterior extension to (K, A) atoms is its own validation gate:
at sdlog = 0 and s = 1 respectively, E-1a and E-1b must reproduce the homogeneous model
bit-identically).
