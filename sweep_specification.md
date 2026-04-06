# Parameter Sweep Specification
# Grant Funding & Scientific Output Model
# ============================================================
# This document defines two sweep designs — Full and Minimal —
# for implementation as a standalone R script by Claude Code.
# The simulation entry point is run_simulation() in app.R.
# ============================================================


## SOURCE FILE

Source lines 1–818 of app.R only. Lines 819 onward are Shiny UI/server
and must be excluded. Strip these patterns from the sourced lines:
  - Lines matching: ^options\(
  - Lines matching: ^suppressPackageStartupMessages
  - Lines matching: ^\s*library\(
  - Lines matching: ^\}\)$   (the closing bracket of the library block)


## BUDGET CONVERSION

Budget is specified in the sweep as a fraction of the pool's total
baseline resources. Use the Pareto MEDIAN (not mean) as the per-researcher
baseline, because the mean diverges as alpha → 1.

  Pareto median = x_min * 2^(1/alpha)

For a given condition:
  R_median  <- r_min * 2^(1 / r_shape)
  B_raw     <- budget_fraction * n * R_median

Pass B_raw as the B argument to run_simulation().


## FIXED DEFAULTS

These are fixed at the values below in all conditions unless a sweep
explicitly varies them.

  k_min           = 1.0
  k_shape         = 2.0    # alpha for K
  r_min           = 1.0
  r_shape         = 2.0    # alpha for R (same as k_shape unless alpha is swept)
  rho_kr          = 0.0
  gamma           = 1.0
  epsilon         = 0.1
  delta           = 1.0    # greedy step size; never vary
  tau_r           = 1.0
  tau_k           = 1.0
  use_resource_signal = TRUE
  n_pre_rounds    = 0
  x_seed          = 0.5
  M               = 1500   # importance samples (see note on Minimal sweep)
  n               = 100


## STRATEGY INDEX

1  No funding
2  Uniform seed
3  Naive (prop. to pubs)
4  Myopic (pubs)
5  Myopic (pubs + grant)
6  Myopic (pubs + seed)
7  Forward (pubs)            — expensive; see performance note
8  Forward (pubs + grant)    — expensive
9  Forward (pubs + seed)     — expensive

Strategies 7–9 each require O(B/delta * n * M) evaluations. At defaults
this is ~7.5M per strategy per run. Profile before including them in
large sweeps. In the Minimal sweep they are excluded.


## SIGNAL NOISE ENCODING

tau_k = Inf means the grant signal is maximally uninformative.
tau_r = Inf means the resource signal is maximally uninformative.

IMPORTANT: Do NOT pass Inf directly to run_simulation. In R,
dnorm(x, sd=Inf, log=TRUE) returns -Inf for all finite x, which causes
the importance weight normalization to produce NaN (via -Inf - (-Inf)).

Instead, map Inf to 1e6 inside run_one() before calling run_simulation:
  tau_k_eff <- ifelse(is.infinite(tau_k), 1e6, tau_k)
  tau_r_eff <- ifelse(is.infinite(tau_r), 1e6, tau_r)

At tau = 1e6 the signal log-likelihood is approximately constant across
all prior samples and cancels in normalization, correctly implementing
"signal off" without numerical issues. Store the original Inf value in
the output columns, not the substituted 1e6.


## ESS NOTE

At tau < 0.25, importance weights may collapse (effective sample size
near 1), degrading Bayesian strategy results. This does not cause errors
but silently corrupts inference. Flag any condition with tau_k < 0.25 or
tau_r < 0.25 in the output with a column ess_risk = TRUE. Do not exclude
these conditions; just mark them for inspection.


## ============================================================
## FULL SWEEP
## ============================================================

Run 30 seeds per condition. All other parameters at defaults unless
the sweep explicitly varies them.

Strategies to run: 1–9 (but see performance note; consider running
7–9 as a separate job after profiling).


### PRIMARY SWEEPS
### (one parameter varied at a time; all others at defaults)

--- n (pool size) ---
Values: 10, 25, 50, 100, 250, 500
Note: n = 1000 excluded — forward strategies become computationally
prohibitive. For n > 250, consider running strategies 1–6 only unless
forward strategies have been profiled at those sizes.

--- alpha (Pareto shape, applied to both K and R) ---
Values: 1.5, 2.0, 2.5, 3.0
Note: alpha = 1.0 excluded from budget-normalized runs because
E[R] is undefined at alpha = 1. Run alpha = 1.0 as a separate
fixed-B condition (B_raw = 50) and do not apply the budget
conversion for that case.

--- r_min / k_min (population composition ratio) ---
Values: 0.25, 0.33, 0.5, 1.0, 2.0, 3.0, 4.0
Implementation: fix k_min = 1.0; set r_min = ratio value.
Note: ratio < 1 → most researchers resource-bottlenecked (grant
funding effective). Ratio > 1 → knowledge-bottlenecked (grant
funding less effective). Empirically, grant applicant pools are
likely resource-bottlenecked; ratio < 1 is the primary regime.

--- budget_fraction (B as fraction of n * Pareto_median(R)) ---
Values: 0.01, 0.05, 0.10, 0.20, 0.40, 0.50, 1.0
Convert to B_raw using formula above before passing to run_simulation.

--- tau_k (grant signal noise SD) ---
Values: 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, Inf
Note: values 0.1 and 0.2 flagged ess_risk = TRUE in output.

--- tau_r (resource signal noise SD) ---
Values: 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, Inf
Note: values 0.1 and 0.2 flagged ess_risk = TRUE in output.

--- x_seed (seed fraction of budget) ---
Values: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
Note: only affects strategies 6 and 9. Run this sweep with strategies
1–6 and 9 (strategy 9 must be included for the seed result to appear).


### INTERACTION SWEEPS
### (two parameters co-varied; all others at defaults)
### These are the key cross-sections for the paper's main claims.

--- tau_k × r_min/k_min ---
tau_k:         0.5, 1.0, 2.0, Inf
r_min/k_min:   0.25, 0.5, 1.0, 2.0, 4.0
Conditions:    4 × 5 = 20

--- tau_k × budget_fraction ---
tau_k:          0.5, 1.0, 2.0, Inf
budget_fraction: 0.05, 0.10, 0.20, 0.40
Conditions:     4 × 4 = 16

--- x_seed × r_min/k_min (seed strategy effectiveness) ---
x_seed:        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
r_min/k_min:   0.25, 0.5, 1.0, 2.0, 4.0
Conditions:    10 × 5 = 50
Strategies:    1–6, 9 (forward seed included; 7–8 optional)


### SENSITIVITY / ROBUSTNESS SWEEPS
### Fixed at defaults for main results; vary for appendix only.

--- gamma (output scaling) ---
Values: 0.5, 1.0, 2.0, 5.0

--- epsilon (knowledge growth rate) ---
Values: 0.0, 0.05, 0.10, 0.20
Note: epsilon = 0 disables knowledge dynamics entirely; useful
as a baseline to show whether dynamics matter at all.

--- rho_kr (K–R correlation) ---
Values: -0.5, 0.0, 0.5
Note: negative correlation means high-K researchers tend to be
under-resourced — the primary motivation for grant funding.

--- n_pre_rounds (naive funding rounds before main experiment) ---
Values: 0, 1, 3


## ============================================================
## MICRO SWEEP
## ============================================================

Goal: complete in under 10 minutes on a single machine. One-at-a-time
design: each parameter swept across three values while all others held
at micro-sweep defaults. No interactions. No sensitivity parameters.
For qualitative orientation only — directional effects, not estimates.

Settings:
  M              = 200    (adequate for tau >= 0.5; qualitative use only)
  seeds          = 7      (sufficient for directional conclusions)
  strategies     = per-block (see below; not all blocks need all strategies)

Micro-sweep defaults (used when a parameter is not being swept):
  n                = 50
  alpha            = 2.0
  r_min/k_min      = 1.0   (r_min = 1.0, k_min = 1.0)
  budget_fraction  = 0.10
  tau_k            = 1.0
  tau_r            = 1.0
  x_seed           = 0.5
  (all other parameters at global fixed defaults)

The default condition (all parameters at micro-sweep defaults) is shared
across all blocks — run it once and reuse.


### PARAMETER SWEEPS

--- n (pool size) ---
Values:     10, 50, 100
Strategies: 1, 2, 4, 5, 6
Purpose:    Does pool size change which strategy wins?

--- alpha (Pareto shape) ---
Values:     1.5, 2.0, 3.0
Strategies: 1, 2, 4, 5, 6
Purpose:    Does output concentration change the value of targeting?
Note:       alpha = 1.5 added vs. originally proposed {2, 3} — the
            high-inequality regime is where targeting matters most.

--- r_min / k_min (population composition) ---
Values:     0.5, 1.0, 2.0
Strategies: 1, 2, 4, 5, 6
Purpose:    Resource-bottlenecked vs. knowledge-bottlenecked populations.
Implementation: fix k_min = 1.0; set r_min = ratio.

--- budget_fraction ---
Values:     0.01, 0.10, 1.0
Strategies: 1, 2, 4, 5
Purpose:    Does budget abundance/scarcity change which strategy wins?
Note:       At budget_fraction = 0.01, B_raw < delta for n ∈ {10, 50}
            (degenerate = TRUE; Bayesian results invalid). At n = 50
            (the default), budget_fraction = 0.01 is degenerate. Retain
            for boundary illustration; exclude degenerate rows from plots.

--- tau_k (grant signal noise) ---
Values:     0.1, 1.0, 10.0
Strategies: 1, 4, 5
Purpose:    Core question: when does the grant signal (S5) beat
            pubs-only (S4)? These three strategies are sufficient.
Note:       tau_k = 0.1 flagged ess_risk = TRUE.

--- tau_r (resource signal noise) ---
Values:     0.1, 1.0, 10.0
Strategies: 1, 4, 5
Purpose:    Value of the resource signal relative to pubs-only baseline.
Note:       tau_r = 0.1 flagged ess_risk = TRUE.

--- x_seed (seed fraction) ---
Values:     0.1, 0.5, 1.0
Strategies: 1, 2, 4, 6
Purpose:    Output-vs-seed-fraction curve. S6 is the only strategy
            affected by x_seed; compare against no-funding (S1),
            uniform (S2), and myopic-pubs baseline (S4).


### MICRO SWEEP TOTALS

  Default condition:    1  (strategies 1, 2, 4, 5, 6)
  n sweep:              2 new conditions
  alpha sweep:          2 new conditions
  r_min/k_min sweep:    2 new conditions
  budget_fraction:      2 new conditions
  tau_k sweep:          2 new conditions
  tau_r sweep:          2 new conditions
  x_seed sweep:         2 new conditions
  Total unique:        15 conditions × 7 seeds = 105 runs

Estimated runtime at ~2–3s per run (n=50, M=200, trimmed strategies):
  105 × 2.5s ≈ 4 minutes single-threaded.
  With parallel::mclapply on 4 cores: ~1 minute.
If actual per-run time exceeds 6s, reduce seeds to 5 (75 runs, ~7 min).


## ============================================================
## MINIMAL SWEEP
## ============================================================

Goal: complete in approximately one hour on a single machine.
Captures the core behavioral structure of the model — signal value,
population composition, budget pressure — with everything else fixed.

Settings that differ from Full sweep:
  M              = 500    (reduced from 1500; adequate for tau >= 0.5)
  n              = 50     (reduced from 100)
  seeds          = 20     (reduced from 30)
  strategies     = 1:6    (forward strategies 7–9 excluded)

Note on M = 500: adequate for the tau values used here (minimum 0.5).
If any Bayesian strategy results look erratic, increase M to 1000 and rerun.


### BLOCK A — Signal value (core question)
Vary tau_k and tau_r jointly. All other parameters at defaults.
Primary purpose: show when and how much grant review improves allocation
over publication records alone.

  tau_k:    0.5, 1.0, 2.0, Inf
  tau_r:    0.5, 1.0, Inf
  Conditions: 4 × 3 = 12


### BLOCK B — Population composition
Vary r_min/k_min. Fix tau_k = 1.0, tau_r = 1.0.
Primary purpose: show how the resource-to-knowledge balance modulates
the value of grant review.

  r_min/k_min:  0.25, 0.5, 1.0, 2.0, 4.0
  Conditions: 5


### BLOCK C — Budget pressure
Vary budget_fraction. Fix tau_k = 1.0, tau_r = 1.0.
Primary purpose: show whether the advantage of better allocation shrinks
at very high or very low budgets.

  budget_fraction: 0.01, 0.05, 0.10, 0.20, 0.40
  Conditions: 5

Note: at budget_fraction = 0.01 with n = 50, B_raw ≈ 0.71 < delta = 1.
This condition will be flagged degenerate = TRUE; Bayesian strategy
results are invalid for that row. It is retained to show the degenerate
boundary but should be excluded from strategy comparison plots.


### BLOCK D — Seed fraction
Vary x_seed. Fix tau_k = 1.0, tau_r = 1.0.
Primary purpose: show the output-vs-seed-fraction curve.
Must include strategy 6 (Myopic pubs + seed) in strategies run.
Strategy 9 (Forward pubs + seed) is excluded here (too slow);
note its absence in output.

  x_seed:  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
  Conditions: 10


### MINIMAL SWEEP TOTALS

  Block A:    12 conditions
  Block B:     5 conditions
  Block C:     5 conditions
  Block D:    10 conditions
  Total:      32 conditions × 20 seeds = 640 runs

Estimated runtime at ~5s per run (n=50, M=500, strategies 1–6):
  640 × 5s ≈ 53 minutes single-threaded.
  With parallel::mclapply across conditions: ~10–15 min on 4 cores.
If actual per-run time exceeds 10s, reduce seeds to 10 or drop Block D.


## ============================================================
## OUTPUT FORMAT
## ============================================================

Save results as a flat data frame. One row per (condition × seed).
Save as both RDS (for R analysis) and CSV (for sharing).
Save incrementally: write after each condition completes so partial
runs are not lost. Recommended: one RDS per condition, merged at end.

All conditions — Micro, Minimal, and Full — write to the same output
directory. The sweep_type and block columns distinguish them. When a
later sweep runs a condition that overlaps with an earlier one, skip it
if the output file already exists (the incremental save logic handles this).

Required columns:

  sweep_type        — "micro", "minimal", or "full"
  block             — block label (e.g. "default", "n", "tau_k")
  seed              — integer seed used
  n, alpha, r_min, k_min, r_min_k_min_ratio
  budget_fraction, B_raw
  tau_k, tau_r
  gamma, epsilon, rho_kr, n_pre_rounds, x_seed
  M
  ess_risk          — TRUE if tau_k < 0.25 or tau_r < 0.25
  degenerate        — TRUE if B_raw < delta (Bayesian strategy results invalid)
  alpha_budget_sub  — alpha used in median formula (= alpha unless alpha=1.0,
                      in which case 1.5)

  Per-strategy output (for each strategy s actually run):
    total_expected_s{s}   — sum of lam1 + lam2 across researchers
    total_output_s{s}     — sum of p1 + p2 (realized Poisson draws)

  Pre-funding population summary (same for all strategies; store once per row):
    mean_K_init           — mean(K_at_start)
    mean_R_init           — mean(R0_at_start)
    mean_D_init           — mean((K-R)/(K+R)) at start
    mean_S_init           — mean(((K-R)/(K+R))^2) at start
    cor_KR_init           — cor(K_at_start, R0_at_start)

  Per-strategy post-funding population summary (for each strategy s run):
    mean_K2_s{s}          — mean(strat$K2)
    mean_R2_s{s}          — mean(strat$R2)
    mean_D2_s{s}          — mean((K2-R2)/(K2+R2))
    mean_S2_s{s}          — mean(((K2-R2)/(K2+R2))^2)
    cor_K2R2_s{s}         — cor(strat$K2, strat$R2)

  Note: these summaries are cheap — they are computed from vectors already
  returned by run_simulation(). No additional simulation cost.

  Derived columns (compute at analysis time, not at run time):
    gain_grant_myopic      — total_expected_s5 - total_expected_s4
    gain_grant_forward     — total_expected_s8 - total_expected_s7 (if run)
    gain_seed_myopic       — total_expected_s6 - total_expected_s4
    gain_seed_vs_grant     — total_expected_s6 - total_expected_s5
    gain_resource_signal   — computed cross-condition: for each row, subtract
                             total_expected_s4 in the tau_r=10 reference row
                             (same other parameters). Note in analysis code,
                             not stored in sweep output.
    best_strategy_idx      — which.max across available total_expected columns
    best_strategy_expected — max across available total_expected columns


## DEFAULT CONDITION (special case)

The default condition (all micro-sweep defaults, see above) is the single
most important condition. It receives special treatment:

  M         = 500   (not 200)
  seeds     = 30    (not 7)
  strategies = 1:9  (all strategies including forward)

This is the basis for Plot 1 (strategy ranking). It is the only condition
in the micro sweep that runs forward strategies (7–9). All other micro
conditions use M=200, 7 seeds, and trimmed strategy sets.

In addition, for the default condition only, save the full researcher-level
data as a separate file: default_condition_researchers.rds. This is a list
with one entry per seed, each containing a data frame with columns:
  researcher_id, K_init, R_init, K2_s{s}, R2_s{s}, g1_s{s}, g2_s{s}
for all strategies s in 1:9. This file supports the K-R scatter plots.


## ============================================================
## IMPLEMENTATION NOTES FOR CLAUDE CODE
## ============================================================

1. Structure the script as:
     source_model()       — strips and sources lines 1–818 of app.R
     build_conditions()   — returns data frame of all parameter combinations
     run_one()            — wraps run_simulation(), handles budget conversion,
                            extracts output columns, tags ess_risk
     run_sweep()          — loops or mclapply over conditions × seeds,
                            saves incrementally

2. Use parallel::mclapply for the outer loop (one call per condition).
   Set mc.cores = parallel::detectCores() - 1.

3. Budget conversion inside run_one():
     R_median <- r_min * 2^(1 / r_shape)
     B_raw    <- budget_fraction * n * R_median
   After computing B_raw, check for degeneracy:
     degenerate <- B_raw < delta
   If degenerate = TRUE, still run the condition (strategies 1–3 are
   unaffected) but flag it in the output with a column degenerate = TRUE.
   Bayesian strategy results (4–9) will show zero allocation and should
   be treated as invalid for those strategies.

4. tau = Inf: do NOT pass Inf directly. See SIGNAL NOISE ENCODING above.
   Map to 1e6 inside run_one() before passing to run_simulation. Store
   original Inf in output columns.

5. alpha = 1.0 special case: the budget median formula is undefined at
   alpha = 1. For alpha = 1.0 conditions only, substitute alpha = 1.5
   in the median formula to compute B_raw:
     R_median_sub <- r_min * 2^(1 / 1.5)
     B_raw        <- budget_fraction * n * R_median_sub
   Record alpha_budget_sub = 1.5 in the output row so this substitution
   is transparent.

6. Incremental save pattern:
     result_file <- file.path(out_dir, sprintf("cond_%05d.rds", cond_id))
     if (!file.exists(result_file)) {
       res <- run_one(...)
       saveRDS(res, result_file)
     }
   This allows resuming interrupted runs without recomputing.

7. Log wall time for the first 5 conditions and print an estimated
   total runtime before proceeding with the full sweep.


## ============================================================
## PLOT INVENTORY
## ============================================================

Plot 1 — Strategy ranking at defaults
  Question:   Which of the 9 strategies produces the most output, and
              by how much?
  Chart type: Dot plot or bar chart; strategies on x-axis ordered by
              total_expected; error bars = 1 SD across seeds.
  Data:       Default condition only. Requires strategies 1–9.
  Columns:    total_expected_s1 through total_expected_s9.
  Source:     Default condition (M=500, 30 seeds).

Plot 2 — Per-parameter effect (7 plots, one per swept parameter)
  Question:   How does strategy rank order change as parameter X varies?
  Chart type: Line plot; x-axis = parameter value; y-axis =
              total_expected; one line per strategy; mean across seeds.
  Data:       Each parameter's sweep block.
  Columns:    total_expected_s{s} for strategies run in that block.
  Note:       The tau_k plot is the most important — it directly shows
              the value of grant review. The r_min/k_min plot shows
              when that value is highest.

Plot 3 — Value of the grant signal
  Question:   How many units of output does adding grant review add,
              and how does this depend on signal quality?
  Chart type: Line plot; x-axis = tau_k; y-axis = gain_grant_myopic
              (and gain_grant_forward if forward strategies were run).
  Data:       tau_k sweep block.
  Columns:    gain_grant_myopic = total_expected_s5 - total_expected_s4.
  Note:       Micro sweep covers myopic only. Forward strategies require
              Minimal sweep.

Plot 4 — Value of the resource signal
  Question:   What does the resource signal (tau_r) add independently?
  Chart type: Line plot; x-axis = tau_r; y-axis = gain relative to
              tau_r = 10 reference (computed cross-condition at analysis
              time).
  Data:       tau_r sweep block.
  Columns:    total_expected_s4, total_expected_s5 across tau_r values.
  Note:       Reference row is tau_r = 10.0 (near-uninformative). Subtract
              at analysis time; not stored in sweep output.

Plot 5 — Grants vs. seeds vs. neither
  Question:   Is it better to use information budget on reviewing
              proposals (grant signal) or to redistribute money uniformly
              (seed)? How does this depend on population composition?
  Chart type: Grouped bar or line; strategies S4, S5, S6 compared;
              optionally faceted by r_min/k_min.
  Data:       Default condition + r_min/k_min sweep block.
  Columns:    total_expected_s4, total_expected_s5, total_expected_s6;
              gain_seed_vs_grant = total_expected_s6 - total_expected_s5.

Plot 6 — K-R distribution shift
  Question:   Which strategies move researchers toward balanced (K≈R)
              vs. concentrating resources on already-balanced researchers?
  Chart type: Scatter of (R, K) pre- and post-funding, colored by
              bottleneck direction D = (K-R)/(K+R); one panel per
              strategy, or summary showing mean_D2 per strategy.
  Data:       Default condition researcher-level file for scatter.
              All conditions for summary statistics.
  Columns (scatter):   default_condition_researchers.rds
  Columns (summary):   mean_D_init, mean_D2_s{s}, mean_S_init,
                       mean_S2_s{s}, cor_KR_init, cor_K2R2_s{s}.
