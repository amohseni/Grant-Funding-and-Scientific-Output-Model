# Section 6 claims ledger (2026-09-06)

Every claim in \S6 v3, formulated precisely, with its modal force, its ground, and its
verification status. "Mechanism test" = verified in this session by direct computation
(verify_s6_mechanisms.py); "owed" = re-derivation from the model repo's sweep data via
verify_all_claims.R or a targeted R run (Mac session). Nothing locks until its row is
green.

| # | Claim (as formulated in the text) | Modal force | Ground | Status |
|---|---|---|---|---|
| C1 | Review is a noisy signal of capability; informativeness governed by tau_K | definitional (\S3 spec) | model definition | OK |
| C2 | Value of review = expected value of its reduction of the funder's uncertainty regarding the capability-resource gap; quantified as expected output with review minus without | definitional | definition + scoring convention (\S3) | OK |
| C3 | A more informative signal brings allocation closer to the optimal allocation, output closer to the complete-information benchmark | "in our simulations, across the sweep" (monotone over tau range) | D-4 sweep | OWED: verify_all_claims.R |
| C4 | With a sharp signal the funder recovers nearly all of the records-only shortfall (17.2 -> 1.7% of no-funding output; corr 0.22->0.95 and 0.34->0.97) | measured at stated parameters | D-4 | OWED: verify_all_claims.R |
| C5a | Records-only allocation correlation stays 0.13-0.18 across rounds at defaults | measured, scoped to simulated horizons ("does not catch up," NOT "never") | D-4 baseline | OWED: verify_all_claims.R |
| C5b | Structural cause: capabilities compound, resources do not accumulate, so output approaches the resource-limited level and carries less and less information about capability | analytic (limit + derivative) | K grows unboundedly (growth -> 2*eps*A*R per round); lambda -> 2AR; dlambda/dK -> 0; per-round Fisher info about K decays 7.9e-2 -> 6.1e-7 over 200 rounds; cumulative info converges (0.44 over 20k rounds) | VERIFIED (mechanism test, this session) |
| C6 | Where capability is heavy-tailed, most of the output a funder can add comes from resourcing the few most capable | "most," heavy-tailed case | share of optimal-allocation output gain from top 10% by K at b=0.1: 0.92 (alpha 1.3), 0.86 (2.0), 0.72 (3.5); 400 draws each | VERIFIED (mechanism test); model-run confirmation owed with C3 |
| C7 | The most capable are easiest to identify; even a noisy signal separates them | comparative ("easiest"; "even a noisy signal") | signal s = K + N(0,tau): top-10% recovery at tau=3: 0.68 (heavy) vs 0.21 (light); middle-pair ordering at tau=3: 0.54 (barely above chance); 2000 draws per cell | VERIFIED (mechanism test) |
| C8 | Review's value nearly constant for noise below tau ~ 2.5, declines beyond; sharpening past a modest level adds little | measured at defaults | D-1 elbow | OWED: verify_all_claims.R |
| C9 | Corollary: fine distinctions among middling applications add comparatively little; heavy investment there MAY be misallocated | "in our model" + "may" | follows from C6 + C7 + C8; C7's middle-pair result (P = 0.50-0.54 correct ordering at realistic noise) is the direct mechanism; uniform-tau caveat stands (rank-local discernment not separately simulated) | MECHANISM VERIFIED; caveat retained |
| C10 | Where capability is spread evenly, review of any informativeness adds little | qualitative comparative | C6/C7 gradients point the right way (0.72 share; 0.21 recovery); output-value magnitude owed | PARTIAL; magnitude owed |
| C11 | Where the budget is ample, information about whom to fund has little left to add | analytic limit + measured | Corollary 2 (proven, appendix) + targeting-value curve | OK (analytic); numbers owed with \S4 set |
| C12 | AUC 0.54 corresponds to tau_K > 20 at defaults | calibration mapping | D-1 | OWED: reproduce mapping (needs the D-1 procedure) |
| C13 | Overtrust can turn review's value negative at moderate tails; attenuates at heavy tails | "can"; cell-noise caveat in footnote | D-2 (-0.53/-0.43/-0.14, each |z| <~ 1.6, consistent across cells) | OWED: verify + keep honesty caveat |
| C14 | Overtrust costs more than undertrust forfeits | asymmetry, measured | D-2 | OWED: verify_all_claims.R |

Notes.
- C5 reformulation (2026-09-06): Aydin's sanity check ("over infinite rounds capabilities
  should become certain since resources grow and output converges to 2AK") fails against
  the model spec in an instructive way: resources do NOT grow across rounds (baseline
  recurs; grants are consumed), while capability DOES compound. So output converges to
  the resource-limited level 2AR, not the capability ceiling 2AK, and identification of
  capability from output degrades rather than improves. "Never catches up" was
  nevertheless too strong as a body claim (the simulations cover finite horizons); the
  text now states the scoped claim plus the structural cause, and the cumulative-Fisher
  computation shows even the infinite-horizon version has an analytic basis (finite
  total information), available if we later want it as a proposition.
- Mechanism tests are this session's direct computations on the model's primitives, not
  runs of the full simulation; they verify mechanisms and directions. Magnitudes quoted
  from the sweeps (C3, C4, C5a, C8, C12, C13, C14) still require the R re-derivation
  before lock.
