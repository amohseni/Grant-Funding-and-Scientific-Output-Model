# Robustness-check handoff — swapping the model's dynamics

*Minimal, self-contained context for a fresh session whose job is to test whether the results are
robust to the model's functional forms — e.g. replacing the current harmonic (Michaelis–Menten)
production/knowledge dynamics with Cobb–Douglas, CES, geometric, etc. Read top-to-bottom once.*

---

## 1. What the project is

A simulation study of how a **science funder should allocate a fixed budget across researchers and over
T rounds** to maximize scientific output, and specifically whether a **forward-looking Bayesian planner**
(learns each researcher's latent state from noisy publications, and anticipates how today's grants raise
future capacity) beats **myopic** and **naive** rules. Pure base-R engine in `model.R`; Shiny front-end
in `app.R`; parameter sweeps in `sweep_T.R` + `sweep.R`; results in `RESULTS.md`.

**Headline results you are stress-testing:**
1. Forward beats myopic, and by more as **knowledge compounds faster (ε↑)** and the **horizon lengthens
   (T↑)** — the advantage is a *compounding* phenomenon (≈0 when ε≈0). Grows to ~+3.4 pubs at ε=0.85, T=10.
2. **Mechanism = back-loading:** the forward planner spends lean early, heavy late (budget center-of-mass
   `b_idx > 0.5`), because output saturates in resources so early resources are wasted on low-knowledge
   researchers. Back-loading intensity correlates 0.97 with the advantage.
3. Many secondary "advantages" (heavy tails, budget/funder scale, noise) were **discretization artifacts**
   of a coarse allocator and vanish under exact optimization.

The robustness question: **do 1 and 2 survive if the production and/or knowledge-growth functions change?**

---

## 2. What the model is (one screen)

- **Agents:** n researchers, each with latent **knowledge K > 0** and **resources R > 0**, drawn from
  Pareto (heavy-tailed) priors, optionally K–R correlated (`draw_initial_population`, `model.R:48`).
- **Output:** each round, researcher i publishes `p ~ Poisson(λ_i)`, where **λ = production(K, R)**.
- **Funding:** a grant `g` adds to resources **this round only** (`R + g`); it is *not* persistent.
- **Knowledge dynamics:** between rounds, **K grows** as a function of (K, R) — this is the *compounding*
  channel and the only thing that couples rounds beyond the shared budget.
- **Funder's information:** never sees (K,R). Infers a **posterior** over each (K,R) from the publication
  history + two noisy signals (on R and on K), via importance sampling with **M=400 atoms**
  (`build_posteriors_hist`, `model.R:890`). All allocation optimizes the **posterior-expected** output.
- **Budget:** `B_total = 2·b·n·E[R]`; non-forward strategies spend an equal tranche `B_total/T` per round;
  forward strategies allocate the whole remaining budget across the horizon.
- **9 strategies** (`out_S1..S9`): S1 none, S2 uniform, S3 naive∝pubs, **S4/5/6 myopic** (pubs /
  +grant-signal / +seed-floor), **S7/8/9 forward** (pubs / +grant / +seed). Headline contrast:
  `fwd_vs_myo_PG = S8 − S5`.
- **Two allocators** behind a switch `allocator = "greedy" | "smooth"` (see §5).

---

## 3. The dynamics live in TWO one-line functions (your robustness levers)

Everything funnels through these (`model.R:65` and `:69`):

```r
lambda_rate      <- function(K, R, gamma)   gamma * (K * R) / (K + R)          # PRODUCTION  (Poisson mean)
update_knowledge <- function(K, R, epsilon) K + epsilon * K * R / (K + R)       # KNOWLEDGE GROWTH (compounding)
```

Both are the **harmonic mean / Michaelis–Menten** form: output needs *both* inputs (→0 if either →0) and
**saturates** in each (doubling R with K fixed less-than-doubles output). Note the current model uses this
same harmonic form for *both* production and knowledge growth. ("Geometric" dynamics are **not** currently
implemented — only harmonic; treat geometric/Cobb–Douglas/CES as variants to add.)

**Why this is the whole game for robustness:** these are the *only* two places the functional form is
defined. `lambda_rate` has **28 call sites**, `update_knowledge` **17** — but they all go through these two
functions, so **editing these two functions changes the form everywhere at once**: the true data-generating
process, the Poisson likelihood (`loglik_pubs`, so posteriors stay consistent automatically), the CE
information term (`ce_reweight_posterior`), and every allocator objective. **One exception — see §5.**

---

## 4. Assumptions worth knowing before you perturb them

- **Both inputs essential + saturating** (harmonic). Cobb–Douglas keeps "both essential" but drops
  saturation (unbounded returns). Linear/additive drops "both essential." These are the interesting axes.
- **Knowledge growth ∝ current output** (`update_knowledge = K + ε·λ/γ`, since `λ/γ = KR/(K+R)`). I.e. you
  learn in proportion to what you produce. Alternatives: geometric (`K·(1+ε)`), Cobb–Douglas learning
  (`K + ε·K^a R^b`), capped/logistic.
- **Grants non-persistent** (resources reset to R₀ each round; only K carries over). This is what makes
  *timing* matter and drives the back-loading result — keep it in mind when interpreting a form change.
- **Objective = posterior-EXPECTED output**, averaged over the M posterior atoms — **not** λ at the
  posterior mean. Any allocator you touch must respect this (`post_lambda_round_t`, `model.R:923`, is
  ground truth: `sum(post$w * lambda_rate(K, R0+g))`).
- **Heavy-tailed (Pareto) K,R priors**; **Poisson** publication noise; **Gaussian** signal noise.
- Budget normalized to mean-E[R]; per-round tranche = B_total/T for non-forward strategies.

---

## 5. THE allocator caveat (read this before changing anything)

There are two allocators, and they differ in how they depend on the functional form:

| Allocator | Marginals via | Form-dependence |
|---|---|---|
| **greedy** (`greedy_round_t`, `plan_forward_ce`) | **finite differences** of the objective (`post_marginal_round_t`, `fwd_marginal`) | **fully form-agnostic** — reads `lambda_rate`/`update_knowledge` only through the objective. Just works after you edit §3. |
| **smooth** — forward (`plan_forward_smooth` → `researcher_grad`) | **central differences** of `fwd_researcher_value` | **also form-agnostic** (differentiates the objective numerically). |
| **smooth** — myopic (`waterfill_round_t`, `model.R`) | **hand-derived analytic KKT** marginal `Σ_m w_m·γK²/(K+R+g)² = ν` | **HARD-CODED for the harmonic form.** Wrong for any other λ. This is the one exception. |

**Consequence:** `waterfill_round_t` bakes in `∂λ/∂g = γK²/(K+R+g)²`, valid only for `λ = γKR/(K+R)`. It's
used by the myopic strategies (S4/5/6) and by `plan_forward_smooth`'s single-round (`H==1`) branch. If you
change the production function and run `allocator="smooth"`, the myopic allocations will be **silently
wrong**.

**Recommendation for robustness sweeps: use `allocator="greedy"` with a high `n_steps` (≥400).** It is
form-agnostic and already validated. Use a *high* n_steps because coarse `n_steps=50` produces a
**discretization artifact** that inflates the forward-vs-myopic contrast at isolated parameter values (this
is exactly what the "smooth allocator" work fixed for the harmonic case — see
`docs/SMOOTH_ALLOCATOR_HANDOFF.md`). At `n_steps≈400-800` greedy converges to the exact optimum.

*(Only if you need the exact/fast smooth path for a new form: re-derive the two lines in `waterfill_round_t`
— `num <- gamma*Wmat*Amat^2` (the `∂λ/∂g` numerator) and the `/d^2`, `/d^3` powers (its g-derivatives) —
for your new λ, or replace `waterfill_round_t` with a generic 1-D concave solver that finite-differences
the objective. The forward smooth planner needs no change.)*

---

## 6. How to run a robustness check (concrete recipe)

**Step 1 — add a form switch (cleanest) or edit in place.** Recommended: parameterize the two functions so
you can sweep the form. Example nesting the current model in a **CES family** (ρ = −1 recovers harmonic,
ρ → 0 is Cobb–Douglas):

```r
# production: CES(rho). rho=-1 -> current harmonic KR/(K+R)*2; rho->0 -> Cobb-Douglas sqrt(K*R); rho=1 -> linear
lambda_rate <- function(K, R, gamma, rho = -1, alpha = 0.5) {
  if (abs(rho) < 1e-8) gamma * (K^alpha * R^(1-alpha))          # Cobb-Douglas limit
  else gamma * (alpha*K^rho + (1-alpha)*R^rho)^(1/rho)
}
# (to reproduce the current model exactly, keep lambda = gamma*K*R/(K+R): that is CES with rho=-1, alpha=0.5,
#  up to a factor of 2 — decide whether to match levels or just the shape.)

update_knowledge <- function(K, R, epsilon, ...) K + epsilon * (lambda_rate(K, R, 1, ...))  # growth ∝ output
# or geometric: K * (1 + epsilon);  or Cobb-Douglas learning: K + epsilon * K^a * R^b
```

Thread the extra params (`rho`, `alpha`, …) through `run_simulation_T(...)` the same way `gamma`/`epsilon`
already flow, OR set them as globals/defaults if you only test one form at a time. **Keep a way to recover
the harmonic baseline** (a regression anchor).

**Step 2 — sanity check** at the baseline form: `Rscript tests/test_T2_reduction.R` should still pass, and
a harmonic run should reproduce `RESULTS.md`. Then flip the form.

**Step 3 — run the headline axes** with the form-agnostic greedy allocator at high granularity:

```r
source("model.R")
r <- run_simulation_T(seed = 1, T_rounds = 5, n = 50, epsilon = 0.85, b = 0.5,
                      M = 400, n_steps = 400, allocator = "greedy", strategies = 1:9)
r$strategies[[8]]$total_expected - r$strategies[[5]]$total_expected   # forward advantage (PG)
r$strategies[[8]]$b_idx                                               # schedule shape (>0.5 = back-load)
```

For a full picture, run a small sweep over the two headline knobs (T ∈ 2..5, ε ∈ {0.05..0.85}) — mirror
`SWEEP_CONFIGS_T$horizon_growth` in `sweep_T.R`, or just loop `run_simulation_T` over a grid and average
over ~100 seeds (use `parallel::mclapply`). `sweep_T.R`'s `run_sweep_T`/`summarize_sweep` machinery
produces the same `_raw/_rawlong/_summary` outputs the analysis tooling expects.

**Step 4 — compare to the harmonic baseline.** The two claims to test:
- **Does the forward advantage still grow with ε×T?** (plot `fwd_vs_myo_PG` vs T, one line per ε).
- **Does the planner still back-load?** (`b_idx_S8 > 0.5`, rising with ε and T; and in the `rawlong` data
  the per-round share `alpha_t` should ramp up over rounds).

If both survive across CES ρ (and under a geometric knowledge law), the results are robust to the dynamics.
Watch for the *predicted* boundary: with **Cobb–Douglas / non-saturating** production, resources are never
"wasted," so the back-loading rationale weakens — a genuinely informative place for the result to change.

---

## 7. Gotchas / lessons already banked

- **Coarse `n_steps` fabricates effects.** Always use `n_steps ≥ 400` with greedy for robustness numbers;
  `n_steps=50` (the old default) inflated contrasts at isolated cells. (Full story:
  `docs/SMOOTH_ALLOCATOR_HANDOFF.md`, `tests/diag_b03_resonance.R`.)
- **Posterior-expected, not at-the-mean** (§4). Easy to get wrong when writing a new marginal.
- **Posteriors auto-stay-consistent:** since `loglik_pubs` calls `lambda_rate`, changing the production
  form also changes the funder's inference model — usually what you want (funder knows the true form). If
  you want a *mis-specified* funder, decouple the likelihood's `lambda_rate` from the true one.
- **The smooth myopic water-fill is harmonic-only** (§5) — don't trust `allocator="smooth"` under a new form
  unless you re-derive it.
- **Keep the harmonic baseline runnable** as a regression anchor; `tests/test_T2_reduction.R` guards the
  T=2 reduction.

---

## 8. Key files & commands

| Purpose | Where |
|---|---|
| Model engine (dynamics, planners, 9 strategies) | `model.R` — `lambda_rate` :65, `update_knowledge` :69, `run_simulation_T` :~1205, `allocate_round` :~1067 |
| Greedy allocator (form-agnostic) | `greedy_round_t`, `plan_forward_ce` in `model.R` |
| Smooth allocator (myopic part = harmonic-specific) | `waterfill_round_t`, `plan_forward_smooth` in `model.R` |
| Sweeps | `sweep_T.R` (`SWEEP_CONFIGS_T`, `run_sweep_T`, `summarize_sweep`), `sweep.R` |
| Baseline results to beat/compare | `RESULTS.md`; data `sweep_results/T_run_smooth/`; figures `figures/smooth_preview/` |
| Allocator/artifact background | `docs/SMOOTH_ALLOCATOR_HANDOFF.md`, `docs/SMOOTH_ALLOCATOR_WRITING_KIT.md` |
| Column dictionary | `T_round_extension/DATA_DICTIONARY.md` |
| Regression anchor | `Rscript tests/test_T2_reduction.R` |

**Up and running in 3 lines:**
```r
source("model.R")                                                    # loads engine
run_simulation_T(seed=1, T_rounds=5, epsilon=0.85, allocator="greedy", n_steps=400)   # one run
# then edit lambda_rate / update_knowledge (model.R:65,69), re-source, re-run, compare.
```
