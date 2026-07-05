# Handoff prompt — re-analyze, validate, and write the paper

Paste everything **below the line** into a fresh Claude conversation, and attach the content listed
in §0. This is self-contained: the model, data schema, findings, and known caveats are all included,
so Claude can validate the results and draft the paper without prior context.

---

## 0. Attach these (content specification)

**If using Claude chat / a Project (no code execution):**
- **The 16 sweep-summary CSVs** — `T_round_extension/data/csv/*.csv` (one row per parameter cell;
  columns are `<metric>_mean` / `<metric>_se`). These ARE the data to analyze.
- **`RESULTS.md`** — my current statistical digest of the findings (to be validated/critiqued).
- **`T_round_extension/DATA_DICTIONARY.md`** — the codebook: every column, every sweep's grid, the
  base parameters.
- **`figures/` PDFs** (8 figures) — the current figure set.
- *(optional but recommended)* **`T_round_extension/STATE_OF_PLAY.md`** — development history incl.
  the one bug that was found/fixed and the methods caveats.

**If using Claude Code (code execution + repo):** point it at the repository instead; it can read
the `.rds` files (`readRDS`), re-run analyses in R, and even re-run simulations (`Rscript
reproduce.R`). The `_raw.rds` (per-trial) and `_rawlong.rds` (per-round schedule) files enable
deeper re-analysis than the summary CSVs.

---

## 1. Your task

You are helping the author finish a paper on this model. Three jobs, in order:

**(A) Re-analyze and validate.** Independently verify the claims in `RESULTS.md` against the attached
data. For each claim: recompute the effect size and z (= `mean / se`) from the summary columns, check
it is actually significant and the direction/magnitude stated, and check internal consistency across
sweeps. **Flag anything overstated, unsupported, or fragile.** Be a skeptical referee, not a
cheerleader. Note where a claim rests on a single cell, a near-zero effect, or a caveated measure.

**(B) Strengthen the analysis.** Propose and (if you can run code) perform additional analyses that
would sharpen the paper: parametric fits (e.g. `signal_fwd` as a function of the knowledge-tail
exponent `k_shape` and signal noise `tau_k`; the ε-slope vs the ε→0 intercept of `fwd_vs_myo_PG`),
regime-boundary characterizations, and any robustness the data supports.

**(C) Draft the paper.** Produce a full draft: abstract, introduction, model, results, discussion,
and a methods/appendix section. Lead with the strongest, most robust claims; relegate fragile ones
to robustness checks or drop them. Use the figures. Write for a meta-science / science-of-science /
economics-of-science audience.

---

## 2. The model (self-contained)

A science **funder** allocates a fixed budget across **n researchers** over **T funding rounds** to
maximize total scientific output. Each researcher *i* has latent **knowledge** Kᵢ and **resources**
Rᵢ, both Pareto-distributed (shape `k_shape` = α_K for K, `r_shape` = α_R for R; smaller shape =
heavier tail = more inequality), optionally correlated via a Gaussian copula (`rho_kr` = ρ_c).

- **Production (output rate):** `λ(K, R) = γ · K·R / (K + R)` — output needs both inputs and is
  bottlenecked by the scarcer one. Total output = Σ over rounds and researchers of λ.
- **Grants are non-persistent:** round-*t* resources are `R₀ + g_t` (the grant applies to that round
  only).
- **Knowledge compounds:** after each round, `K ← K + ε · K·R/(K+R)` (`epsilon` = ε, the growth
  rate) — funding permanently raises a researcher's knowledge and thus future output.
- **Observation:** the funder sees noisy publications `p_t ~ Poisson(λ)` each round, plus two
  one-time signals: a **resource signal** σ_R (noise SD `tau_r`) and a **grant / peer-review signal**
  σ_K (noise SD `tau_k`). It infers each researcher's (K, R) by Bayesian importance sampling.
- **Budget:** total `B_total = 2·b·n·E[R]` (E[R] = mean of the resource distribution), split evenly
  across the T rounds for non-forward strategies; forward strategies plan the whole budget freely.

**The 9 strategies** (the analysis compares these):

| # | Name | Grant signal? | Plans ahead? | Seed floor? |
|---|------|---------------|--------------|-------------|
| S1 | No funding | – | – | – |
| S2 | Uniform (flat every round) | – | – | – |
| S3 | Naive (∝ observed pubs) | – | – | – |
| S4 | Myopic (pubs) | no | no | – |
| S5 | Myopic (pubs + grant) | **yes** | no | – |
| S6 | Myopic (pubs + seed) | no | no | round-1 floor |
| S7 | Forward (pubs) | no | **yes** | – |
| S8 | Forward (pubs + grant) | **yes** | **yes** | – |
| S9 | Forward (pubs + seed) | no | **yes** | round-1 floor |

**Forward** planners are a receding-horizon **certainty-equivalent (CE)** Bayesian planner: they plan
the whole remaining horizon (anticipating knowledge compounding and, one step ahead, the information
their grants reveal), execute the current round, then re-plan next round with the newly observed
pubs. **Myopic** planners optimize only the current round. The reported output is *expected* output
(Σλ, deterministic given the drawn state), so contrasts are low-variance.

## 3. The data (what's in the CSVs)

16 sweeps, each a CSV with **one row per parameter cell**, **200 trials/cell**. Columns:
the varied parameter(s), `n_trials`, and for every metric `X` two columns `X_mean` and `X_se`
(SE = sd/√200). Significance ≈ `|mean| > 2·se`.

**Per-strategy metrics** (j = 1..9): `out_Sj` (total expected output), `alpha_Sj` (round-1 share of
spend). Forward-only: `b_idx_S{7,8,9}` (schedule center-of-mass; 0.5 = even across rounds, >0.5 =
front-loaded), `gini_g1_S{5,7,8,9}` (round-1 funding concentration). `rho_s` (realized Spearman K–R
correlation).

**Contrasts** (the workhorses): `signal_fwd` = S8−S7 (**value of the grant signal**); `signal_myo` =
S5−S4; `fwd_vs_myo_PG` = S8−S5 (**value of planning ahead** — the horizon headline); `fwd_vs_myo_P` =
S7−S4; `fwd_vs_myo_PS` = S9−S6; `seed_fwd` = S9−S7 (**value of a seed floor**); `seed_myo` = S6−S4;
`optimal_vs_naive` = S8−S3.

**The 16 sweeps** (full grids in `DATA_DICTIONARY.md` §3): `horizon_growth` (T×ε — the headline),
`horizon_noise`, `horizon_scale`, `signal_precision`, `signal_value` (α_K×τ_K), `funder_scale`,
`seed_value`, `alpha_regime`, `pre_rounds`, `regime_map`, `correlation` (ρ_c×α_K), `pop_size`,
`resource_noise`, `tail_map` (α_K×α_R), `info_value` (T×τ_K at ε≈0), `horizon_long` (T=5..10).
Unvaried parameters take base defaults (`DATA_DICTIONARY.md` §4): T=2, n=50, ε=0.1, b=0.5, α_K=α_R=2,
τ_K=τ_R=1, ρ_c=0, γ=1, M=400 posterior samples.

## 4. Findings to validate (condensed from `RESULTS.md`)

1. **Forward beats myopic — a (horizon × growth) phase diagram** (`horizon_growth`): `fwd_vs_myo_PG`
   is positive everywhere and grows with ε; mechanism is front-loading (`b_idx_S8` rises with ε).
2. **Signal value = inequality × precision** (`signal_value`, `signal_precision`): `signal_fwd` from
   ~21 (heavy tail + sharp signal) to ~0 (light tail + noisy).
3. **Knowledge inequality dominates; resource inequality secondary** (`tail_map`).
4. **Signal value survives K–R correlation** (`correlation`): large even at ρ_s = 0.78.
5. **Uniform seeding never helps** (`seed_value`): `seed_fwd` < 0 everywhere.
6. **Decomposition:** forward's edge = compounding (∝ε) + a small information component that peaks at
   T=3 (`info_value`, ε≈0).
7. **Long horizons** (`horizon_long`): advantage grows monotonically with T, modest, no saturation.
8. Robustness nulls: resource-signal noise irrelevant; no finite-n artifact; budget scale doesn't
   drive the forward gain.

## 5. Known caveats to scrutinize (do not gloss over these)

- **Greedy granularity.** Allocation is a greedy fill in `δ = B/n_steps` steps. The sweeps use
  `n_steps=50`, which is granularity-stable for **T ≤ 5** but too coarse beyond it — it inflated the
  high-T advantage ~4× until `horizon_long` was re-run at `n_steps=400`. If you extend any horizon
  result past T=5, granularity matters. (Verify the T≤5 sweeps don't secretly depend on it.)
- **Certainty-equivalence approximation.** Forward is CE-approximate, not a full stochastic program.
  The CE-vs-scenario valuation gap is < 0.4% across T=2–5 and doesn't grow with T (so results are not
  CE artifacts) — but this is a methods point to state, and worth sanity-checking the argument.
- **Budget confound in `tail_map`.** The α_R axis co-varies with the budget (B ∝ E[R]), so heavier
  resource tail ⇒ more inequality *and* more money. Don't over-claim a pure "resource inequality"
  effect from it.
- **`pre_rounds` scale artifact.** `signal_fwd` rising with pre-round data is partly because naive
  pre-round funding inflates the absolute scale of K, R, λ (and hence absolute contrasts) — interpret
  relatively.
- **`info_value` T=3 peak** is a feature of the planner's one-step information anticipation, not a
  general "information grows with horizon" law. Characterize it honestly.
- A subtle planner bug (stochastic resampling corrupting the forward planner's finite-difference
  marginals) was found and fixed; all results here are post-fix. Pre-fix data (`sweep_results/legacy/`)
  must not be used.

## 6. Output

Deliver: (A) a validation memo (claim-by-claim: holds / overstated / unsupported, with the recomputed
numbers); (B) any additional analyses/fits; (C) the paper draft. Ask for the raw/rawlong data or code
access if a validation needs more than the summary cells.
