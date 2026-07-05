# Prompt for Claude Chat — statistical analysis of Sweep B (optimal alpha)

Copy everything below the line into Claude Chat. Attach the seven CSVs from
`sweep_results/` (paths in §3) so the model can actually read the data.

---

I'm analyzing results from a science-funding simulation. A bug was just fixed
and the headline sweep ("Sweep B") was re-run. I want your help doing the
statistical analysis on the new results — testing hypotheses, characterizing
the regime structure, fitting parametric models to `alpha*(tau_k, epsilon)`,
and quantifying the before/after change.

## 1. Model (minimum context)

A funder has a total budget of `2B` to split across two rounds of grants to
`n=50` researchers. Per-researcher publication rate is

    lambda(K, R) = gamma * K * R / (K + R)

where `K` is researcher knowledge and `R` is researcher resources (their own
+ any grant). Between rounds, knowledge updates:

    K2 = K + epsilon * K * R / (K + R)

so `epsilon` is the knowledge-compounding rate. The funder observes noisy
signals: `tau_k` is the standard deviation of the grant-quality signal
(higher = noisier). Posteriors are particle-based and refresh between rounds
with round-1 publications.

Nine allocation strategies are run head-to-head per trial; we mostly care
about three forward-Bayes variants:

  - S7  Forward(pubs)
  - S8  Forward(pubs + grant signal)        ← main forward strategy
  - S9  Forward(pubs + seed)

and their myopic counterparts:

  - S4  Myopic(pubs)
  - S5  Myopic(pubs + grant signal)         ← main myopic strategy
  - S6  Myopic(pubs + seed)

Plus three baselines: S1 No-funding, S2 Uniform-seed, S3 Naive (prop. to
pubs).

**The decision variable of interest** is

    alpha = sum(g1) / (sum(g1) + sum(g2))

i.e. the fraction of total budget spent in round 1. Only **forward** planners
optimize over alpha (they have two rounds of lookahead); myopic planners
always spend `B` per round → `alpha = 0.5` by construction.

## 2. The bug and the fix

**Before the fix**: `run_forward_bayes()` spent `B` in round 1 and `B` in
round 2, so `alpha = 0.5` was forced for forward strategies — making the
hypothesis test untestable. The pre-fix `alpha_S{7,8,9}` is exactly 0.5 in
every trial (sd = 0 across 1000 trials).

**After the fix**: the forward planner has total budget `2B` and chooses
alpha by argmax over a grid `seq(0.05, 0.95, 0.05)` of the
posterior-expected two-round output `V(alpha)`. Round-1 greedy is run for
each candidate alpha with a consistent round-2 lookahead, then the chosen
alpha is realized under the true state. Same change applied to S7, S8, S9.
Myopic strategies untouched.

## 3. Files (please load these)

All paths are inside `sweep_results/` of the simulation repo. I'll attach
them as CSVs.

### Sweep B (the headline test — re-run with seeds 1:50 after the fix)

  - **`optimal_alpha_v2_raw.csv`** — 1000 rows (20 cells × 50 seeds).
    One row per `(tau_k, epsilon, seed)`. Columns:
      - Design: `tau_k`, `epsilon`, `seed`
      - Per-strategy expected output: `out_S1` … `out_S9`
      - Per-strategy alpha: `alpha_S1` … `alpha_S9`
        - `alpha_S1` is NA (no-funding has no allocation)
        - `alpha_S{2,3,4,5,6}` is always 0.5 (myopic/naive/uniform — by
          construction)
        - `alpha_S{7,8,9}` is the planner's chosen alpha — this is the
          dependent variable of interest
      - Derived comparisons (also present in summary as `_mean`/`_se`):
        - `fwd_vs_myo_P  = out_S7 - out_S4`
        - `fwd_vs_myo_PG = out_S8 - out_S5`   ← main forward-vs-myopic test
        - `fwd_vs_myo_PS = out_S9 - out_S6`
        - `signal_myo  = out_S5 - out_S4`     (myopic value of grant signal)
        - `signal_fwd  = out_S8 - out_S7`     (forward value of grant signal)
        - `seed_myo    = out_S6 - out_S4`
        - `seed_fwd    = out_S9 - out_S7`
        - `optimal_vs_naive = out_S8 - out_S3`

  - **`optimal_alpha_v2_summary.csv`** — 20 rows (cell-level aggregate).
    Columns: `tau_k`, `epsilon`, `n_trials`, then `_mean` and `_se` (SE of
    the mean across the 50 seeds) for every numeric column in raw.

### Sweep B pre-fix (for before/after comparison)

  - **`optimal_alpha_v1_PREFIX_raw.csv`** — 1000 rows. Same schema; same
    grid; produced with the buggy code where every `alpha_S{7,8,9} = 0.5`.
  - **`optimal_alpha_v1_PREFIX_summary.csv`** — 20 rows.

### Sweep A spot-check (forward_vs_myopic_main, 9 cells × 15 seeds, post-fix)

Sampled from a 5×5 (tau_k, epsilon) grid — corners + middle.

  - **`forward_vs_myopic_main_spotcheck_raw.csv`** — 135 rows. Columns:
    `tau_k`, `epsilon`, `seed`, `out_S5`, `out_S8`, `alpha_S8`.
  - **`forward_vs_myopic_main_spotcheck_summary.csv`** — 9 rows: cell means
    plus `fwd_vs_myo_PG_se`.

### Sweep A pre-fix summary

  - **`forward_vs_myopic_main_v1_PREFIX_summary.csv`** — 25 rows from the
    full 5×5 pre-fix Sweep A. Same schema as the Sweep B summary.

### Grid coordinates

  - **Sweep B**:  `tau_k ∈ {0.1, 1, 5, 20}`,
                  `epsilon ∈ {0.05, 0.20, 0.40, 0.65, 0.95}`
  - **Sweep A pre-fix**: `tau_k ∈ {0.1, 0.3, 1, 3, 10}`,
                         `epsilon ∈ {0.05, 0.15, 0.30, 0.55, 0.85}`
  - **Sweep A spot-check**: subset `tau_k ∈ {0.1, 1, 10}`,
                            `epsilon ∈ {0.05, 0.30, 0.85}`

Common base params (held fixed across all sweeps): `n=50`, `b=0.5`
(dimensionless budget; `B = b * n * E[R]`), `gamma=1`, `tau_r=1`,
`k_min=1, k_shape=2`, `r_min=1, r_shape=2`, `rho_kr=0`, `x_seed=0.5`,
`M=500` (importance-sampling particles), `delta=1` (greedy step).

## 4. What's already known

Headline numbers I've extracted from `optimal_alpha_v2_summary.csv`:

  - `alpha_S8` ranges **[0.404, 0.564]** across the 20 cells.
  - alpha < 0.5 in **11/20 cells (55%)**. The original hypothesis required
    ≥90%, so the cells clause **FAILS**.
  - `alpha_S8` is strictly decreasing in `epsilon` at every fixed `tau_k`
    (per-row range 0.10–0.13). The monotonicity clause **HOLDS**.
  - Mean alpha at `epsilon=0.05`: **0.538** (slight front-loading).
    Mean alpha at `epsilon=0.95`: **0.420** (clear backloading).
  - Strongest backloading at the high-epsilon / high-tau_k corner
    (α=0.404 at eps=0.95, tau_k=20). Strongest front-loading at the
    low-epsilon / low-tau_k corner (α=0.564 at eps=0.05, tau_k=0.1).

Sweep A spot-check finding: pre-fix `fwd_vs_myo_PG_mean ≈ -0.33` uniformly
across all 25 cells. Post-fix mean narrows to **−0.146**; forward beats
myopic at high epsilon (+0.117 at tau_k=0.1, eps=0.85) but loses at low
epsilon. None of the per-cell shifts cross 2 SE at 15 seeds × M=500.

## 5. What I'd like you to do

Please do the statistical analysis I haven't done yet. In order of
importance:

  1. **Quantify the alpha(tau_k, epsilon) surface.** Fit a parametric model
     to `alpha_S8` as a function of `(tau_k, epsilon)` using the
     **raw** (run-level) data, not the cell means, so SEs are honest.
     Reasonable starting candidates: linear in `epsilon` and `log(tau_k)`
     with an interaction; or a logistic transform of alpha if you think
     the unit interval matters at this range. Report coefficient
     estimates, SEs, R², and a residual diagnostic.

  2. **Test the monotonicity claim statistically.** I asserted alpha is
     monotonically decreasing in epsilon at every tau_k. Do a per-tau_k
     trend test (e.g. Jonckheere–Terpstra against the ordered alternative,
     or just a slope test from the regression) and report p-values.

  3. **Locate the alpha = 0.5 contour.** For each tau_k, find the epsilon
     at which `E[alpha_S8] = 0.5` (a 1D root of the fitted model, with a
     CI). This is the "front-load vs backload" boundary and is more
     interesting than the 55%-of-cells summary.

  4. **Before/after comparison.** Using v1 (pre-fix) and v2 (post-fix)
     raw data:
     - Are the `out_S8` distributions different? (They should be: the
       post-fix planner has a real choice variable so total output
       should differ at non-trivial cells.) A paired test by seed within
       cell is the right move.
     - Is the gap `out_S8 - out_S5` significantly different post-fix
       vs. pre-fix? Where is the effect concentrated?

  5. **Sweep A re-verdict.** Given the 9-cell spot-check, what's the
     posterior probability that forward beats myopic at the high-epsilon
     end of the grid? With 15 seeds × M=500, the per-cell SE is large;
     a Bayesian one-sample test or even just a confidence interval on
     the gain at each cell would be useful. Tell me how many seeds I'd
     need to detect a 0.2-unit gain at α=0.05 / power=0.8.

  6. **Anything surprising in the raw data** that I might have missed:
     pathologies in particular cells, signs that the alpha-grid is too
     coarse (e.g. all chosen alphas land exactly on grid points by
     definition, but distribution over grid points across seeds may be
     informative), heteroskedasticity, outliers, etc.

  7. A short note on whether the hypothesis is fairly described as
     "FAIL" given the literal threshold but qualitatively supported, or
     whether you'd phrase the verdict differently for a paper.

I'd like your answer as a written report with the key numbers and figures
inline, plus the code blocks you used so I can re-run them. Use Python or R
— pick whichever lets you do the regression / contour-finding most cleanly.

Don't restate the model setup or repeat what I've already told you;
go straight to the analysis.
