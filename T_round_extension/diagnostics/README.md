# diagnostics/ — how the bug was found and confirmed (provenance)

These scripts are the trail that uncovered, localized, and root-caused the SIR-resampling bug in
the forward planner (full story in `../STATE_OF_PLAY.md` §5). Kept for provenance and in case the
investigation needs to be reproduced or extended. Run from the **project root**.

**Run them against the *fixed* `app.R`** to see the corrected (stable) behavior; the numbers quoted
below are what they produced on the *buggy* code, to illustrate the symptom.

| File | What it does | What it found (on buggy code) |
|------|--------------|-------------------------------|
| `pg_focused.R` | Focused `fwd_vs_myo_PG(T×ε)` at n=50, more seeds. First close look at the horizon result. | Non-monotonic PG — the first hint something was off. |
| `ce_tax_vs_T.R` | Decomposes the CE valuation "tax" per round vs T (point-estimate vs scenario-integrated continuation value). Intended to test whether the high-T decline was a CE artifact. | CE tax ~0 (CE is faithful) **but** revealed PG flips sign erratically with n — the real clue. |
| `granularity_check.R` | `fwd_vs_myo_PG(T=5)` vs greedy granularity `n_steps` ∈ {50,100,200,400} across n ∈ {40,50,60}. | The erratic n-dependence *shrinks* with finer granularity but PG *diverges* negative — signature of a discretization-amplified bug, not a real effect. |
| `convergence_check.R` | Convergence of PG vs `n_steps` up to 800, for ε ∈ {0.3, 0.85}; also reports `S7−S4` (pubs-only forward). | PG diverges (to −7); **`S7−S4` is stable** — localized the bug to the *grant-signal* forward path. |

The decisive bug-vs-real test (the forward greedy scoring below myopic on its *own* CE objective at
fine granularity) and the fix confirmation (disabling SIR → convergence) were run inline during the
investigation; the logic is documented in `../STATE_OF_PLAY.md` §5.

**Bottom line:** the bug was a stochastic resampling step (`ce_reweight_posterior`'s SIR branch)
that corrupted the planner's finite-difference marginals. Fixed by making the reweight
deterministic. On the fixed code these scripts show granularity-*stable*, positive forward gains.
