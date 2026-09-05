# Errata and reconciliation: report + model write-up

*2026-08-05. Applies to `Optimal_Funding_Strategies_for_Scientific_Output.pdf` (the report) and
`Optimal_Funding_Strategies_for_Scientific_Output_1.pdf` (the model write-up). The LaTeX sources are
not in the model repo or the for-claude repo; apply these at the source (Overleaf or wherever they
live). Each item: what is wrong, and the exact replacement. Items marked [Kevin] carry a TBD that
must survive into the draft he edits.*

---

## E1. The γ collision (both documents)

The report uses γ for factor productivity (λ = γKR/(K+R)); the model write-up uses γ for the CES
substitution exponent and A for productivity. Same symbol, two meanings, documents read together.

**Canonical, from here on:**

| object | paper symbol | code name | relation |
|---|---|---|---|
| factor productivity | A | `gamma` | A = gamma/2 under the family normalization Λ(x,x) = x |
| CES exponent | γ | `ces_gamma` | admissible γ ≤ 0; harmonic γ = −1 |
| elasticity of substitution | σ | derived | σ = 1/(1−γ) |

The write-up's convention wins (it is the economics convention and the family needs both symbols).
The report's every "γ" as productivity becomes "A" with the factor-of-2 mapping below. The code
keeps `gamma` for productivity internally; `ces_gamma` is the exponent; the mapping is documented
once in the computational appendix and never again.

## E2. Report equations restated under the family convention

With λ = A·Λ_γ(K,R) and Λ_{−1}(K,R) = 2KR/(K+R):

- Report Part I: "λ = γKR/(K+R)" → "λ = A·Λ(K,R), where Λ(K,R) = 2KR/(K+R) is the harmonic mean"
  (or introduce the family first and instantiate; the write-up §2 already does this correctly).
- Report marginal return: ∂λ/∂g = γK²/(K+R+g)² → **2A·K²/(K+R+g)²**.
- Report footnote 2 and the gap rule: c = √(γ/ν) − 1 → **c = √(2A/ν) − 1**. The rule g* = cK − R is
  unchanged. (Check: 2AK²/(K+R+g)² = ν gives K+R+g = K√(2A/ν).)
- Knowledge dynamics K ← K + ελ: consistent under either convention since λ is λ; no change beyond
  the symbol.

Nothing numerical changes anywhere: at the sweeps' gamma = 1, A = 1/2, and every reported quantity
is a strategy contrast unaffected by the relabeling.

## E3. One term for K (both documents)

The report says "talent"; the write-up says "quality/knowledge." One concept, one term.
**Recommendation: talent**, everywhere including the dynamics ("talent compounds"), because the
report already uses it consistently and "talent–resource gap" is the phrase the paper is remembered
by. Word-level call is Aydin's; the alternative ("knowledge", matching the K symbol) costs the
memorable phrase. Whichever wins, the write-up's "Quality/knowledge K" bullet and every later
"knowledge" referring to K are renamed to match.

## E4. Pareto shape symbols (write-up)

Write-up uses α (knowledge shape) and β (resource shape); the report and code use α_K, α_R
(`k_shape`, `r_shape`). Canonical: **α_K, α_R**. Replace α → α_K, β → α_R in the write-up's §8
functional forms and parameter table.

## E5. The ρ collision (write-up)

The write-up uses ρ for grant persistence (accumulating-baseline option) while the report and code
use ρ for the K–R correlation (`rho_kr`). Canonical: **ρ = the K–R (talent–resource) correlation.**
The persistence parameter, if the accumulating baseline is kept at all (see E8), is renamed φ.

## E6. Eq (10) contradicts the implemented model (write-up; substantive)

Eq (10) states a per-round constraint Σ_i g_i^(t) ≤ B with "budget is not transferable across
rounds; total spend is at most TB." The implemented forward strategies allocate the whole remaining
budget across the horizon; the report's back-loading schedules (round-1 share 0.08 to 0.37 at T=5)
are exactly budget transfer across rounds. Eq (10) as written describes only the non-forward
strategies, and "total spend at most TB" describes nothing in the code: the total purse is fixed and
does NOT grow with T.

**Replacement for the Funding section:**

> The funder has a total budget B_total = 2·b·n·E[R], fixed across the horizon. (The constant 2
> preserves the two-round model's purse; total funds do not grow with T, so the horizon varies the
> timing of spending, never its amount.) Non-forward strategies (S1–S6) spend an equal tranche
> B_total/T each round: Σ_i g_i^(t) ≤ B_total/T. Forward strategies (S7–S9) may allocate any
> nonnegative schedule satisfying Σ_t Σ_i g_i^(t) ≤ B_total, implemented as receding-horizon
> control: each round the full remaining budget is re-planned over the remaining horizon and only
> the current round's allocation is executed. A run option budget_ref ∈ {R, K} selects the
> reference mean in the purse (E[R], the historical default; E[K] decouples the purse from baseline
> resources for resource-poverty experiments).

This matters for interpretation, not just bookkeeping: holding the purse fixed in T is what makes
the horizon results pure timing results.

## E7. Default-parameter table (write-up) does not match the data

The write-up's table (T=10, n=100, B=50, τ=1, x_seed=0.5, ...) matches no reported run. The
report's Appendix C base is what every number in both documents was generated from. **Replace the
write-up's defaults column with:** T=2, n=50, ε=0.1, b=0.5 (so B_total = 2·b·n·E[R] = 100 at
E[R]=2), α_K=α_R=2, k_min=r_min=1, τ_K=τ_R=1, ρ=0, A=1/2 (code gamma=1), M=400, x_seed=0.25,
n_pre=0, budget_ref=R. Two code-signature footnotes: the function default is M=200 and every sweep
passes M=400; the function default x_seed=0.25 while the write-up currently says 1/2.

## E8. Spec describes unimplemented machinery as if available (write-up) [Kevin]

- **Lookahead depth h.** The write-up presents h ∈ {1, ..., T−t+1} with "intermediate h
  interpolates and is the practical setting for large T." Only the endpoints exist in code: h=1
  (myopic) and h = T−t+1 (forward). Replacement: "Only the endpoint policies are implemented and
  reported: h = 1 (myopic) and the full remaining horizon (forward). Intermediate depths are an
  extension. [TBD — Kevin: cut the h-interpolation claim or we implement it.]"
- **Baseline options 2 and 3** (redrawn; accumulating with persistence). Only the constant baseline
  exists in code. Replacement: mark both as specified extensions, not run options. [TBD — Kevin:
  cut or implement before submission.]

## E9. The robustness sentence asserts runs that have not happened (write-up; substantive) [Kevin]

"Every result reported in the main text at γ = −1 is replicated in the appendix at γ ∈ {0, −3, −∞}"
is false today: no CES sweep has been run, and the CES family itself is not yet in the code (the
current model.R hardcodes the harmonic form). **Replacement:** "Every result reported in the main
text at γ = −1 is to be replicated at γ ∈ {0, −3, −∞}. [TBD — Kevin: the sweep is specified in
docs/SWEEP_HANDOFF_2026-08-05.md, Package B, and has not yet been run; no claim in the main text
may rest on it until the results memo lands.]" The same correction applies wherever the report or
draft implies the production-form robustness is done.

## E10. The T-scaling overclaim (report)

"...roughly linearly in ε and quadratically in T, because a round-1 boost pays out in every later
round and itself compounds." The quadratic claim has no fitted exponent behind it. **Replacement:**
"increasing in both: approximately linear in ε at fixed T, and superlinear in T at fixed ε over the
range tested (T ≤ 5), with the fitted exponent reported in the appendix." The exponent comes free
from existing data; the sweep handoff's item A5 computes it. Keep the mechanism clause; drop the
exponent claim until the number exists.

## E11. Withdrawal of a previously flagged inconsistency (Claude's error)

Session 2 flagged a contradiction between the report's main text V.2 and the p.4 margin note on
what drives back-loading ("W4"). On re-reading, both say the same thing: free force C dominating
paid force B. The p.4 margin note reads "free knowledge C over paid B"; V.2's text agrees. **There
is no textual inconsistency; the flag was a misreading and is withdrawn.** What remains true and
open is evidential: both forces scale with ε, so no existing sweep discriminates the C-attribution
from a B-attribution. That is exactly what Package C of the sweep handoff tests. The attribution
text needs no edit now; it needs the sweep's verdict before it is drafted into the paper.

## E12. Seed-floor section (report) [held]

Part VI's framing ("uniform seed floors reduce output") and its counter-reading ("floors are nearly
free") are both unlicensed until Package A returns. No replacement text yet; the section is frozen,
and whichever claim survives D3/D4 gets written at that point. Tracked in
drafts/grant-funding/notes-seed-floor.md.

---

## Order of application

E1–E5 are mechanical renamings; apply in one pass. E6–E7 are one-paragraph and one-table
replacements. E8–E10 insert TBD flags that must remain visible to Kevin. E11 requires no edit.
E12 waits on the sweeps. None of E1–E12 blocks the sweep handoff, which already uses the canonical
naming on the code side.
