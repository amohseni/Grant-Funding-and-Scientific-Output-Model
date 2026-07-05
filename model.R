# model.R
# ============================================================
# Grant-funding / scientific-output model — consolidated ENGINE
# ------------------------------------------------------------
# The single source of the simulation model. Pure base R (no Shiny /
# ggplot). Sourced by app.R (the Shiny UI), by sweep.R / sweep_T.R, and
# by the test / diagnostic scripts.
#
# Contents:
#   * Model primitives, Gaussian-copula population, importance-sampling
#     posteriors, and the certainty-equivalent (CE) forward machinery.
#   * The 2-round strategies + run_simulation() — the v5 model, kept as
#     the reference that the T-round path reproduces bit-identically at T=2.
#   * The T-round generalization (run_simulation_T + receding-horizon CE
#     planner) — the DEFAULT model; v5 is its T=2 special case.
#
# Lineage: consolidates the model code formerly split between app.R (v5,
# 2-round) and simulate_T.R (T-round). See T_round_extension/ for the
# analysis package and its STATE_OF_PLAY.md for development history.
# ============================================================

STRATEGY_NAMES <- c(
  "No funding",
  "Uniform seed",
  "Naive (prop. to pubs)",
  "Myopic (pubs)",
  "Myopic (pubs + grant)",
  "Myopic (pubs + seed)",
  "Forward (pubs)",
  "Forward (pubs + grant)",
  "Forward (pubs + seed)"
)

# ============================================================
# Core model code
# ============================================================

# ----- Distributions -----

# Power law (Pareto Type I) sampler
rpareto <- function(n, x_min, shape) {
  u <- runif(n)
  x_min / (u)^(1 / shape)
}

# ----- Gaussian copula for K-R correlation -----

draw_initial_population <- function(n, k_min, k_shape, r_min, r_shape, rho_kr) {
  if (abs(rho_kr) < 1e-10) {
    u_k <- runif(n)
    u_r <- runif(n)
  } else {
    z1 <- rnorm(n)
    z2 <- rho_kr * z1 + sqrt(1 - rho_kr^2) * rnorm(n)
    u_k <- pnorm(z1)
    u_r <- pnorm(z2)
  }
  K0 <- k_min / (1 - u_k)^(1 / k_shape)
  R0 <- r_min / (1 - u_r)^(1 / r_shape)
  list(K0 = K0, R0 = R0)
}

# ----- Model primitives -----

lambda_rate <- function(K, R, gamma) {
  gamma * (K * R) / (K + R)
}

update_knowledge <- function(K, R, epsilon) {
  # Harmonic mean (Michaelis-Menten) growth: K grows in proportion to current
  # output rate λ = K·R/(K+R). No ceiling. Always non-negative.
  K + epsilon * K * R / (K + R)
}

# ----- Bottleneck measures -----

compute_bottleneck <- function(K, R) {
  # D: direction of bottleneck (−1 fully resource-bottlenecked, +1 fully K-bottlenecked)
  # S: severity = D²
  D <- (K - R) / (K + R)
  S <- D^2
  list(D = D, S = S)
}

# ----- Signals -----

draw_signals <- function(K0, R0, tau_r, tau_k) {
  n <- length(K0)
  sigma_r <- R0 + rnorm(n, 0, tau_r)
  sigma_k <- K0 + rnorm(n, 0, tau_k)
  list(sigma_r = sigma_r, sigma_k = sigma_k)
}

# ----- Observation model (log-likelihood) -----

loglik_pubs <- function(p_obs, K, R, gamma) {
  lam <- lambda_rate(K, R, gamma)
  dpois(p_obs, lambda = pmax(lam, 1e-12), log = TRUE)
}

# Continuous extension of the Poisson log-pmf used for CE reweighting in
# the Forward planner. The CE planner conditions on an imagined
# p̄_1 = E^(0)[λ_1], which is non-integer in general. The natural extension
# log f(x | λ) = x·log(λ) − λ − log Γ(x+1) coincides with dpois at integer
# x and varies smoothly off the integers.
loglik_pubs_continuous <- function(p_value, K, R, gamma) {
  lam <- pmax(lambda_rate(K, R, gamma), 1e-12)
  p_value * log(lam) - lam - lgamma(p_value + 1)
}

loglik_resource_signal <- function(sigma_r, R, tau_r) {
  dnorm(sigma_r, mean = R, sd = tau_r, log = TRUE)
}

loglik_grant_signal <- function(sigma_k, K, tau_k) {
  dnorm(sigma_k, mean = K, sd = tau_k, log = TRUE)
}

# ----- Posterior inference (importance sampling) -----

posterior_samples_single <- function(
    p_obs, sigma_r, sigma_k,
    M, k_min, k_shape, r_min, r_shape, gamma,
    tau_r, tau_k,
    use_resource_signal, use_grant_signal,
    p1_obs = NULL, g1_i = 0
) {
  K_s <- rpareto(M, k_min, k_shape)
  R_s <- rpareto(M, r_min, r_shape)
  
  ll <- loglik_pubs(p_obs, K_s, R_s, gamma)
  if (!is.null(p1_obs)) {
    ll <- ll + loglik_pubs(p1_obs, K_s, R_s + g1_i, gamma)
  }
  if (use_resource_signal) {
    ll <- ll + loglik_resource_signal(sigma_r, R_s, tau_r)
  }
  if (use_grant_signal) {
    ll <- ll + loglik_grant_signal(sigma_k, K_s, tau_k)
  }
  
  ll <- ll - max(ll)
  w <- exp(ll)
  w <- w / sum(w)
  
  list(K0 = K_s, R0 = R_s, w = w)
}

build_posteriors <- function(
    p_obs, sigma_r, sigma_k,
    M, k_min, k_shape, r_min, r_shape, gamma,
    tau_r, tau_k,
    use_resource_signal, use_grant_signal,
    p1_obs = NULL, g1 = NULL
) {
  n <- length(p_obs)
  posts <- vector("list", n)
  for (i in seq_len(n)) {
    posts[[i]] <- posterior_samples_single(
      p_obs = p_obs[i],
      sigma_r = sigma_r[i],
      sigma_k = sigma_k[i],
      M = M,
      k_min = k_min, k_shape = k_shape,
      r_min = r_min, r_shape = r_shape,
      gamma = gamma,
      tau_r = tau_r, tau_k = tau_k,
      use_resource_signal = use_resource_signal,
      use_grant_signal    = use_grant_signal,
      p1_obs = if (!is.null(p1_obs)) p1_obs[i] else NULL,
      g1_i   = if (!is.null(g1)) g1[i] else 0
    )
  }
  posts
}

# ----- Expected output for one round given (K, R, gamma) -----

posterior_expected_lambda <- function(post, g, gamma) {
  R <- post$R0 + g
  lam <- lambda_rate(post$K0, R, gamma)
  sum(post$w * lam)
}

posterior_marginal_lambda <- function(post, g, dg, gamma) {
  (posterior_expected_lambda(post, g + dg, gamma) -
     posterior_expected_lambda(post, g, gamma)) / dg
}

# Round-2 expected lambda with knowledge update
posterior_expected_lambda_r2 <- function(post, g1, g2, gamma, epsilon) {
  K2 <- update_knowledge(post$K0, post$R0 + g1, epsilon)
  R2 <- post$R0 + g2
  lam <- lambda_rate(K2, R2, gamma)
  sum(post$w * lam)
}

posterior_marginal_lambda_r2 <- function(post, g1, g2, dg, gamma, epsilon) {
  (posterior_expected_lambda_r2(post, g1, g2 + dg, gamma, epsilon) -
     posterior_expected_lambda_r2(post, g1, g2,      gamma, epsilon)) / dg
}

# Two-round expected output (for forward-looking strategies)
posterior_expected_two_round <- function(post, g1, g2, gamma, epsilon) {
  K1 <- post$K0
  R1 <- post$R0 + g1
  lam1 <- lambda_rate(K1, R1, gamma)
  
  K2 <- update_knowledge(K1, R1, epsilon)
  R2 <- post$R0 + g2
  lam2 <- lambda_rate(K2, R2, gamma)
  
  sum(post$w * (lam1 + lam2))
}

# ----- CE machinery for the Forward planner -----
#
# A fully optimal Bayesian planner would compute
#   V_2*(π0, g1) = E^(0)[ V_2*(π1) | g1 ],
# where π1 is the (random) post-round-1 posterior. We use a
# certainty-equivalent (CE) approximation: replace the expectation
# over realizations of p1 with a single deterministic conditioning
# at the prior-predictive mean p̄_1,i(g1) = E^(0)[λ(K_i, R_{0,i}+g1)].
#
# CE bias is upward on α (CE underestimates the option value of round 2,
# so puts too much in round 1). Acceptable for a heuristic planner;
# diagnostics below include ESS to flag where importance weights collapse.

effective_sample_size <- function(w) {
  s <- sum(w)
  if (s <= 0) return(0)
  w_norm <- w / s
  1 / sum(w_norm^2)
}

# CE reweight: imagine observing p̄_1,i at current g1[i] and update π0 → π̄1.
# Keeps the same (K, R) sample atoms; only the importance weights change.
#
# DETERMINISTIC by construction — no resampling. Earlier versions did a
# stochastic SIR resample when the reweighted ESS fell below M/4. That broke
# the Forward planner: its marginals are finite differences that call this
# twice (at g1 and g1+dg); when SIR fired, the two calls drew DIFFERENT random
# atoms, so the difference was dominated by O(1) resampling noise rather than
# the O(dg) signal. The effect was mild at 2 rounds / coarse granularity but
# diverged with horizon and finer granularity, and was grant-signal-specific
# (a sharp K posterior lowers ESS, firing SIR). Keeping the reweight fully
# deterministic makes the finite-difference marginals exact; sample degeneracy
# at heavy tails is controlled by raising M instead (see ESS assertion).
ce_reweight_posterior <- function(post0, g1_i, gamma, ess_floor = NULL) {
  K_s <- post0$K0
  R_s <- post0$R0
  w_s <- post0$w

  # Imagined observation: prior-predictive mean of round-1 publications
  lam_at_g <- lambda_rate(K_s, R_s + g1_i, gamma)
  p_bar <- sum(w_s * lam_at_g)

  # Reweight by Poisson(p_bar | λ) — weights only, atoms fixed
  log_lik <- loglik_pubs_continuous(p_bar, K_s, R_s + g1_i, gamma)
  log_w <- log(pmax(w_s, 1e-300)) + log_lik
  log_w <- log_w - max(log_w)
  w_new <- exp(log_w); w_new <- w_new / sum(w_new)

  list(K0 = K_s, R0 = R_s, w = w_new, p_bar = p_bar)
}

# Type-A marginal under CE: effect of adding δ to g1[i].
# Captures (a) direct round-1 expected return at π0,
#          (b) knowledge dynamics into round 2 (K2 depends on g1),
#          (c) information value of g1 via the imagined posterior π̄1.
# Total derivative is computed by finite difference: re-reweight at g1+δ.
posterior_marginal_forward_A <- function(post0, g1, g2, dg, gamma, epsilon) {
  # Round-1 expectation always evaluated at π0 (the actual round-1 posterior
  # before any imagined observation); this isolates Force A from Forces D/E.
  v1_base <- posterior_expected_lambda(post0, g1,        gamma)
  v1_plus <- posterior_expected_lambda(post0, g1 + dg,   gamma)
  
  # Round-2 expectation evaluated at the CE-reweighted posterior π̄1
  pb_base <- ce_reweight_posterior(post0, g1,      gamma)
  pb_plus <- ce_reweight_posterior(post0, g1 + dg, gamma)
  v2_base <- posterior_expected_lambda_r2(pb_base, g1,      g2, gamma, epsilon)
  v2_plus <- posterior_expected_lambda_r2(pb_plus, g1 + dg, g2, gamma, epsilon)
  
  ((v1_plus + v2_plus) - (v1_base + v2_base)) / dg
}

# Type-B marginal under CE: effect of adding δ to g2[i].
# Round-2 only; π̄1 is unchanged by g2.
posterior_marginal_forward_B <- function(post_bar1, g1, g2, dg, gamma, epsilon) {
  v_base <- posterior_expected_lambda_r2(post_bar1, g1, g2,      gamma, epsilon)
  v_plus <- posterior_expected_lambda_r2(post_bar1, g1, g2 + dg, gamma, epsilon)
  (v_plus - v_base) / dg
}

# ============================================================
# Allocation strategies
# ============================================================

# ----- Strategy 1: No funding -----
allocate_no_funding <- function(n, B) {
  list(g1 = rep(0, n), g2 = rep(0, n))
}

# ----- Strategy 2: Uniform seed -----
allocate_uniform <- function(n, B) {
  g <- rep(B / n, n)
  list(g1 = g, g2 = g)
}

# ----- Strategy 3: Naive (proportional to pubs) -----
allocate_naive <- function(p_round, B) {
  w <- p_round + 1e-6
  B * (w / sum(w))
}

allocate_naive_two_round <- function(n, B, p_init, K_current, R0, g1_alloc, gamma) {
  g1 <- allocate_naive(p_init, B)
  R1 <- R0 + g1
  lam1 <- lambda_rate(K_current, R1, gamma)
  g2 <- allocate_naive(lam1, B)
  list(g1 = g1, g2 = g2, lam1 = lam1)
}

# ----- Greedy one-round (with optional baseline) -----
allocate_greedy_one_round <- function(posts, B, delta, gamma, init_g = NULL) {
  n <- length(posts)
  g <- if (is.null(init_g)) rep(0, n) else init_g
  remaining <- B
  dg <- delta
  
  mv <- vapply(seq_len(n), function(i) {
    posterior_marginal_lambda(posts[[i]], g[i], dg, gamma)
  }, numeric(1))
  
  while (remaining >= delta) {
    i_star <- which.max(mv)
    g[i_star] <- g[i_star] + delta
    remaining <- remaining - delta
    mv[i_star] <- posterior_marginal_lambda(posts[[i_star]], g[i_star], dg, gamma)
  }
  
  g
}

# ----- Greedy round-2 using knowledge-updated marginals -----
allocate_greedy_round2 <- function(posts, B, delta, gamma, g1_vec, epsilon) {
  n <- length(posts)
  g <- rep(0, n)
  remaining <- B
  dg <- delta
  
  mv <- vapply(seq_len(n), function(i) {
    posterior_marginal_lambda_r2(posts[[i]], g1_vec[i], g[i], dg, gamma, epsilon)
  }, numeric(1))
  
  while (remaining >= delta) {
    i_star <- which.max(mv)
    g[i_star] <- g[i_star] + delta
    remaining <- remaining - delta
    mv[i_star] <- posterior_marginal_lambda_r2(
      posts[[i_star]], g1_vec[i_star], g[i_star], dg, gamma, epsilon
    )
  }
  
  g
}

# ----- Myopic Bayes: Strategies 4, 5 -----
run_myopic_bayes <- function(
    K_true, R0, p_init, sigma_r, sigma_k,
    B, delta, gamma, epsilon,
    M, k_min, k_shape, r_min, r_shape,
    tau_r, tau_k,
    use_resource_signal, use_grant_signal
) {
  n <- length(K_true)
  
  posts1 <- build_posteriors(
    p_obs = p_init, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal    = use_grant_signal
  )
  g1 <- allocate_greedy_one_round(posts1, B, delta, gamma)
  
  R1        <- R0 + g1
  lam1_true <- lambda_rate(K_true, R1, gamma)
  p1        <- rpois(n, pmax(lam1_true, 1e-12))
  K2_true   <- update_knowledge(K_true, R1, epsilon)
  
  posts2 <- build_posteriors(
    p_obs = p_init, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal    = use_grant_signal,
    p1_obs = p1, g1 = g1
  )
  g2 <- allocate_greedy_round2(posts2, B, delta, gamma, g1, epsilon)
  
  R2        <- R0 + g2
  lam2_true <- lambda_rate(K2_true, R2, gamma)
  p2        <- rpois(n, pmax(lam2_true, 1e-12))
  
  list(
    g1 = g1, g2 = g2,
    p1 = p1, p2 = p2,
    lam1 = lam1_true, lam2 = lam2_true,
    K1 = K_true, K2 = K2_true,
    R1 = R1, R2 = R2,
    total_output   = sum(p1 + p2),
    total_expected = sum(lam1_true + lam2_true)
  )
}

# ----- Myopic Bayes with seed: Strategy 6 -----
run_myopic_seed <- function(
    K_true, R0, p_init, sigma_r, sigma_k,
    B, delta, gamma, epsilon, x_seed,
    M, k_min, k_shape, r_min, r_shape,
    tau_r, tau_k, use_resource_signal
) {
  n <- length(K_true)
  
  g1_uniform <- rep(x_seed * B / n, n)
  remaining_B <- (1 - x_seed) * B
  
  posts1 <- build_posteriors(
    p_obs = p_init, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal    = FALSE
  )
  g1 <- allocate_greedy_one_round(posts1, remaining_B, delta, gamma,
                                  init_g = g1_uniform)
  
  R1        <- R0 + g1
  lam1_true <- lambda_rate(K_true, R1, gamma)
  p1        <- rpois(n, pmax(lam1_true, 1e-12))
  K2_true   <- update_knowledge(K_true, R1, epsilon)
  
  posts2 <- build_posteriors(
    p_obs = p_init, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal    = FALSE,
    p1_obs = p1, g1 = g1
  )
  g2 <- allocate_greedy_round2(posts2, B, delta, gamma, g1, epsilon)
  
  R2        <- R0 + g2
  lam2_true <- lambda_rate(K2_true, R2, gamma)
  p2        <- rpois(n, pmax(lam2_true, 1e-12))
  
  list(
    g1 = g1, g2 = g2,
    p1 = p1, p2 = p2,
    lam1 = lam1_true, lam2 = lam2_true,
    K1 = K_true, K2 = K2_true,
    R1 = R1, R2 = R2,
    total_output   = sum(p1 + p2),
    total_expected = sum(lam1_true + lam2_true)
  )
}

# ----- Forward-looking Bayes (CE PLANNER, CLOSED-LOOP) -----
#                                            Strategies 7, 8, 9
#
# Planning: for each researcher i, an imagined posterior π̄1[i] is
# maintained alongside the pre-round posterior π0[i]. π̄1[i] is obtained
# by reweighting π0[i] with the Poisson likelihood at the prior-predictive
# mean p̄_1,i(g1[i]) = E^(0)[λ(K_i, R_{0,i} + g1[i])]. This lets the
# planner internalize Forces D and E (information value of round-1
# grants) in addition to A, B, C — making it forward-looking in the
# full Bayesian sense.
#
# Greedy operates over 2n candidate moves under remaining budget
# (2B − Σg1 − Σg2). Type-A marginals evaluate round-1 at π0 and round-2
# at π̄1; Type-B marginals evaluate round-2 only at π̄1. After a Type-A
# move on i*, π̄1[i*] is refreshed; Type-B moves leave π̄1 untouched.
#
# Closed-loop execution: after observing the real p1, the posterior is
# rebuilt incorporating the actual Poisson likelihood (π1) and g2 is
# re-greedied under remaining budget. The planning g2 is discarded;
# only g1 carries over from planning to execution.
#
# CE is an approximation. The exact Bayesian planner would integrate
# over realizations of p1 (Monte Carlo or stochastic programming),
# which is too expensive for interactive Shiny. CE biases α upward
# (Jensen on round-2 option value). The MC reference is a separate
# validation script.
#
# Strategy variants:
#   use_grant_signal = TRUE   → strategy 8 (pubs + grant)
#   use_seed = TRUE           → strategy 9 (pubs + seed, no grant)
#   neither                   → strategy 7 (pubs only)
run_forward_bayes <- function(
    K_true, R0, p_init, sigma_r, sigma_k,
    B, delta, gamma, epsilon,
    M, k_min, k_shape, r_min, r_shape,
    tau_r, tau_k,
    use_resource_signal, use_grant_signal,
    x_seed = 0, use_seed = FALSE
) {
  n <- length(K_true)
  
  # ---- Planning stage ----
  posts0 <- build_posteriors(
    p_obs = p_init, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal    = use_grant_signal
  )
  
  total_budget <- 2 * B
  dg <- delta
  
  if (use_seed) {
    g1 <- rep(x_seed * B / n, n)        # uniform seed counts toward round 1
  } else {
    g1 <- rep(0, n)
  }
  g2_plan <- rep(0, n)
  remaining <- total_budget - sum(g1) - sum(g2_plan)
  
  # CE-imagined posteriors at current g1 (one per researcher)
  posts_bar1 <- lapply(seq_len(n), function(i) {
    ce_reweight_posterior(posts0[[i]], g1[i], gamma)
  })
  
  mv_A <- vapply(seq_len(n), function(i) {
    posterior_marginal_forward_A(posts0[[i]], g1[i], g2_plan[i], dg, gamma, epsilon)
  }, numeric(1))
  mv_B <- vapply(seq_len(n), function(i) {
    posterior_marginal_forward_B(posts_bar1[[i]], g1[i], g2_plan[i], dg, gamma, epsilon)
  }, numeric(1))
  
  while (remaining >= delta) {
    iA <- which.max(mv_A); vA <- mv_A[iA]
    iB <- which.max(mv_B); vB <- mv_B[iB]
    
    if (vA >= vB) {
      i_star <- iA
      g1[i_star] <- g1[i_star] + delta
      # Type-A move shifts g1[i_star], so refresh π̄1[i_star]
      posts_bar1[[i_star]] <- ce_reweight_posterior(
        posts0[[i_star]], g1[i_star], gamma
      )
    } else {
      i_star <- iB
      g2_plan[i_star] <- g2_plan[i_star] + delta
      # Type-B move leaves π̄1[i_star] unchanged
    }
    remaining <- remaining - delta
    
    mv_A[i_star] <- posterior_marginal_forward_A(
      posts0[[i_star]], g1[i_star], g2_plan[i_star], dg, gamma, epsilon
    )
    mv_B[i_star] <- posterior_marginal_forward_B(
      posts_bar1[[i_star]], g1[i_star], g2_plan[i_star], dg, gamma, epsilon
    )
  }
  
  # ---- Round-1 realization under TRUE state ----
  R1        <- R0 + g1
  lam1_true <- lambda_rate(K_true, R1, gamma)
  p1        <- rpois(n, pmax(lam1_true, 1e-12))
  K2_true   <- update_knowledge(K_true, R1, epsilon)
  
  # ---- Closed-loop reallocation ----
  # Rebuild posterior with the actual p1 observation, then re-greedy g2.
  posts1 <- build_posteriors(
    p_obs = p_init, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal    = use_grant_signal,
    p1_obs = p1, g1 = g1
  )
  remaining_budget <- max(total_budget - sum(g1), 0)
  g2 <- allocate_greedy_round2(posts1, remaining_budget, delta, gamma, g1, epsilon)
  
  # ---- Round-2 realization under TRUE state ----
  R2        <- R0 + g2
  lam2_true <- lambda_rate(K2_true, R2, gamma)
  p2        <- rpois(n, pmax(lam2_true, 1e-12))
  
  total_spend <- sum(g1) + sum(g2)
  alpha <- if (total_spend > 0) sum(g1) / total_spend else NA_real_
  
  list(
    g1 = g1, g2 = g2,
    p1 = p1, p2 = p2,
    lam1 = lam1_true, lam2 = lam2_true,
    K1 = K_true, K2 = K2_true,
    R1 = R1, R2 = R2,
    total_output   = sum(p1 + p2),
    total_expected = sum(lam1_true + lam2_true),
    alpha          = alpha
  )
}

# ============================================================
# Main simulation pipeline
# ============================================================

run_simulation <- function(
    seed = 1,
    n = 100,
    k_min = 1.0, k_shape = 2.0,
    r_min = 1.0, r_shape = 2.0,
    rho_kr = 0,
    gamma = 1.0,
    epsilon = 0.3,
    b = 0.5,
    n_steps = 50,
    tau_r = 1.0,
    tau_k = 1.0,
    use_resource_signal = TRUE,
    n_pre_rounds = 0,
    x_seed = 0.5,
    M = 1500,
    strategies = 1:9,
    verbose = FALSE
) {
  set.seed(seed)
  
  # Compute absolute budget and greedy step size
  E_R   <- r_min * r_shape / (r_shape - 1)
  B     <- b * n * E_R
  delta <- B / n_steps
  
  # Step 1: Draw initial conditions
  pop <- draw_initial_population(n, k_min, k_shape, r_min, r_shape, rho_kr)
  K0_init <- pop$K0
  R0_init <- pop$R0
  
  initial_state <- list(K = K0_init, R = R0_init)
  
  # Step 2: Generate initial publication counts
  lam_init <- lambda_rate(K0_init, R0_init, gamma)
  p_init <- rpois(n, pmax(lam_init, 1e-12))
  
  # Step 3: Run pre-rounds (Naive allocation)
  K_current <- K0_init
  R0_current <- R0_init
  p_cumul <- p_init
  
  if (n_pre_rounds > 0) {
    for (r in seq_len(n_pre_rounds)) {
      g_pre <- allocate_naive(p_cumul, B)
      R_round <- R0_current + g_pre
      lam_pre <- lambda_rate(K_current, R_round, gamma)
      p_round <- rpois(n, pmax(lam_pre, 1e-12))
      p_cumul <- p_cumul + p_round
      K_current <- update_knowledge(K_current, R_round, epsilon)
    }
  }
  
  post_preround_state <- list(K = K_current, R = R0_current)
  
  # Step 4: Draw signals from post-pre-round state
  sigs <- draw_signals(K_current, R0_current, tau_r, tau_k)
  sigma_r <- sigs$sigma_r
  sigma_k <- sigs$sigma_k
  
  # Step 5: Run main 2 rounds for each selected strategy
  results <- list()
  
  for (s in strategies) {
    set.seed(seed * 1000)
    
    if (s == 1) {
      alloc <- allocate_no_funding(n, B)
      R1 <- R0_current + alloc$g1
      lam1 <- lambda_rate(K_current, R1, gamma)
      p1 <- rpois(n, pmax(lam1, 1e-12))
      K2 <- update_knowledge(K_current, R1, epsilon)
      R2 <- R0_current + alloc$g2
      lam2 <- lambda_rate(K2, R2, gamma)
      p2 <- rpois(n, pmax(lam2, 1e-12))
      results[[s]] <- list(
        name = STRATEGY_NAMES[s],
        g1 = alloc$g1, g2 = alloc$g2,
        p1 = p1, p2 = p2,
        lam1 = lam1, lam2 = lam2,
        K1 = K_current, K2 = K2,
        R1 = R1, R2 = R2,
        total_output   = sum(p1 + p2),
        total_expected = sum(lam1 + lam2)
      )
      
    } else if (s == 2) {
      alloc <- allocate_uniform(n, B)
      R1 <- R0_current + alloc$g1
      lam1 <- lambda_rate(K_current, R1, gamma)
      p1 <- rpois(n, pmax(lam1, 1e-12))
      K2 <- update_knowledge(K_current, R1, epsilon)
      R2 <- R0_current + alloc$g2
      lam2 <- lambda_rate(K2, R2, gamma)
      p2 <- rpois(n, pmax(lam2, 1e-12))
      results[[s]] <- list(
        name = STRATEGY_NAMES[s],
        g1 = alloc$g1, g2 = alloc$g2,
        p1 = p1, p2 = p2,
        lam1 = lam1, lam2 = lam2,
        K1 = K_current, K2 = K2,
        R1 = R1, R2 = R2,
        total_output   = sum(p1 + p2),
        total_expected = sum(lam1 + lam2)
      )
      
    } else if (s == 3) {
      nv <- allocate_naive_two_round(n, B, p_cumul, K_current, R0_current, NULL, gamma)
      g1 <- nv$g1
      R1 <- R0_current + g1
      lam1 <- lambda_rate(K_current, R1, gamma)
      p1 <- rpois(n, pmax(lam1, 1e-12))
      K2 <- update_knowledge(K_current, R1, epsilon)
      g2 <- allocate_naive(p1, B)
      R2 <- R0_current + g2
      lam2 <- lambda_rate(K2, R2, gamma)
      p2 <- rpois(n, pmax(lam2, 1e-12))
      results[[s]] <- list(
        name = STRATEGY_NAMES[s],
        g1 = g1, g2 = g2,
        p1 = p1, p2 = p2,
        lam1 = lam1, lam2 = lam2,
        K1 = K_current, K2 = K2,
        R1 = R1, R2 = R2,
        total_output   = sum(p1 + p2),
        total_expected = sum(lam1 + lam2)
      )
      
    } else if (s == 4) {
      res_s <- run_myopic_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal    = FALSE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s
      
    } else if (s == 5) {
      res_s <- run_myopic_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal    = TRUE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s
      
    } else if (s == 6) {
      res_s <- run_myopic_seed(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon, x_seed,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k, use_resource_signal
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s
      
    } else if (s == 7) {
      res_s <- run_forward_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal    = FALSE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s
      
    } else if (s == 8) {
      res_s <- run_forward_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal    = TRUE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s
      
    } else if (s == 9) {
      res_s <- run_forward_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal    = FALSE,
        x_seed = x_seed, use_seed = TRUE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s
    }
  }
  
  # Attach descriptive alpha (round-1 share) to every strategy.
  # Forward-looking strategies already compute alpha; other strategies
  # get a descriptive value here (0.5 when sum(g1)=sum(g2)=B; NA when
  # total spend is zero, i.e. strategy 1).
  for (s in seq_along(results)) {
    r <- results[[s]]
    if (is.null(r) || !is.null(r$alpha)) next
    total_spend <- sum(r$g1) + sum(r$g2)
    r$alpha <- if (total_spend > 0) sum(r$g1) / total_spend else NA_real_
    results[[s]] <- r
  }
  
  # Compute bottleneck measures
  bn_initial  <- compute_bottleneck(initial_state$K, initial_state$R)
  bn_post_pre <- compute_bottleneck(post_preround_state$K, post_preround_state$R)
  
  list(
    params = list(
      seed = seed, n = n,
      k_min = k_min, k_shape = k_shape,
      r_min = r_min, r_shape = r_shape,
      rho_kr = rho_kr,
      gamma = gamma, epsilon = epsilon,
      b = b, B = B, E_R = E_R, n_steps = n_steps, delta = delta,
      tau_r = tau_r, tau_k = tau_k,
      use_resource_signal = use_resource_signal,
      n_pre_rounds = n_pre_rounds,
      x_seed = x_seed, M = M
    ),
    initial_state = initial_state,
    post_preround_state = post_preround_state,
    signals = list(sigma_r = sigma_r, sigma_k = sigma_k),
    p_init = p_init,
    p_cumul = p_cumul,
    K_at_start = K_current,
    R0_at_start = R0_current,
    strategies = results,
    bottleneck = list(
      initial = bn_initial,
      post_preround = bn_post_pre
    )
  )
}


# ============================================================
# ============================================================
# T-ROUND GENERALIZATION  (the default model)
# ------------------------------------------------------------
# Generalizes the 2-round planners above to T rounds. At T=2 the
# forward / myopic / seed planners reduce to run_simulation() bit-
# identically (see tests/test_T2_reduction.R).
# ============================================================
# ============================================================

# ---- Named constants ----
LAMBDA_FLOOR <- 1e-12   # floor on Poisson rate (matches app.R)

# ============================================================
# Posterior with full observation history
# ============================================================
# Importance-sampling posterior over (K0, R0) for each researcher,
# conditioning on baseline pubs p_init, every past round's pubs
# p_hist[[s]] under that round's grant g_hist[[s]], and the one-time
# signals. rpareto draw order matches v5 posterior_samples_single so
# the T=2 path consumes RNG identically.
#
# INPUTS  p_init  baseline pub counts (n)
#         sigma_r, sigma_k  one-time signals (n)
#         p_hist, g_hist    lists of past-round pub / grant vectors
#                           (NULL or empty for the first planning round)
# OUTPUT  list of n posteriors, each list(K0, R0, w)
build_posteriors_hist <- function(
    p_init, sigma_r, sigma_k,
    M, k_min, k_shape, r_min, r_shape, gamma,
    tau_r, tau_k, use_resource_signal, use_grant_signal,
    p_hist = NULL, g_hist = NULL
) {
  n <- length(p_init)
  posts <- vector("list", n)
  for (i in seq_len(n)) {
    K_s <- rpareto(M, k_min, k_shape)
    R_s <- rpareto(M, r_min, r_shape)
    ll  <- loglik_pubs(p_init[i], K_s, R_s, gamma)
    if (!is.null(p_hist) && length(p_hist) > 0) {
      for (s in seq_along(p_hist)) {
        ll <- ll + loglik_pubs(p_hist[[s]][i], K_s, R_s + g_hist[[s]][i], gamma)
      }
    }
    if (use_resource_signal) ll <- ll + loglik_resource_signal(sigma_r[i], R_s, tau_r)
    if (use_grant_signal)    ll <- ll + loglik_grant_signal(sigma_k[i], K_s, tau_k)
    ll <- ll - max(ll)
    w  <- exp(ll); w <- w / sum(w)
    posts[[i]] <- list(K0 = K_s, R0 = R_s, w = w)
  }
  posts
}

# ============================================================
# Myopic machinery (round-t greedy, K compounded through history)
# ============================================================
# Expected round-t λ under `post`, given this researcher's executed
# past grants g_hist_i (rounds 1..t-1) and current-round grant g_cur.
# t=1 (empty history) == posterior_expected_lambda;
# t=2 (history [g1])  == posterior_expected_lambda_r2.
post_lambda_round_t <- function(post, g_hist_i, g_cur, gamma, epsilon) {
  K <- post$K0
  if (length(g_hist_i) > 0) {
    for (gh in g_hist_i) {
      K <- update_knowledge(K, post$R0 + gh, epsilon)
    }
  }
  sum(post$w * lambda_rate(K, post$R0 + g_cur, gamma))
}

post_marginal_round_t <- function(post, g_hist_i, g_cur, dg, gamma, epsilon) {
  (post_lambda_round_t(post, g_hist_i, g_cur + dg, gamma, epsilon) -
     post_lambda_round_t(post, g_hist_i, g_cur,      gamma, epsilon)) / dg
}

# Greedy fill of `budget` for a single round, maximizing round-t expected
# λ. g_hist is a list-of-length-(t-1) of full grant vectors (n) already
# executed; init_g optionally seeds the round (seed floor). delta = step.
# t=1                 == allocate_greedy_one_round;
# t=2 with g_hist=[g1] == allocate_greedy_round2.
greedy_round_t <- function(posts, budget, delta, gamma, epsilon,
                           g_hist = list(), init_g = NULL) {
  n <- length(posts)
  g <- if (is.null(init_g)) rep(0, n) else init_g
  remaining <- budget
  ghi <- function(i) vapply(g_hist, function(gh) gh[i], numeric(1))
  mv <- vapply(seq_len(n), function(i)
    post_marginal_round_t(posts[[i]], ghi(i), g[i], delta, gamma, epsilon),
    numeric(1))
  while (remaining >= delta) {
    i_star <- which.max(mv)
    g[i_star] <- g[i_star] + delta
    remaining <- remaining - delta
    mv[i_star] <- post_marginal_round_t(posts[[i_star]], ghi(i_star),
                                        g[i_star], delta, gamma, epsilon)
  }
  g
}

# ============================================================
# Forward (receding-horizon CE) planner
# ============================================================
# Value of researcher i's H-round plan g_path under the CE rule:
#   round 1 expectation at post0 (π0);
#   rounds 2..H expectation at post_bar1 (π̄1 = CE-reweight of π0 at g_path[1]),
#   with K compounding along g_path through post_bar1's atoms.
# g_hist_i are this researcher's EXECUTED past grants (rounds before the plan);
# K is compounded through them first, so the plan's "round 1" starts from the
# correct current knowledge. At t=1 (empty history) H=2 == v1@π0 + v2@π̄1 (v5).
fwd_researcher_value <- function(post0, post_bar1, g_hist_i, g_path, gamma, epsilon) {
  H <- length(g_path)
  K0 <- post0$K0
  for (gh in g_hist_i) K0 <- update_knowledge(K0, post0$R0 + gh, epsilon)
  val <- sum(post0$w * lambda_rate(K0, post0$R0 + g_path[1], gamma))
  if (H >= 2) {
    Kb <- post_bar1$K0
    for (gh in g_hist_i) Kb <- update_knowledge(Kb, post_bar1$R0 + gh, epsilon)
    Kb <- update_knowledge(Kb, post_bar1$R0 + g_path[1], epsilon)
    for (s in 2:H) {
      val <- val + sum(post_bar1$w * lambda_rate(Kb, post_bar1$R0 + g_path[s], gamma))
      if (s < H) Kb <- update_knowledge(Kb, post_bar1$R0 + g_path[s], epsilon)
    }
  }
  val
}

# Marginal of adding dg to slot s of researcher i's plan.
#   s == 1: round-1 move — recomputes π̄1 at g1 and g1+dg (info value),
#           re-evaluates the whole path. At H=2 == forward_A.
#   s >= 2: future-round move — π̄1 fixed, only rounds ≥ s change via
#           compounding. At H=2 (s=2) == forward_B.
fwd_marginal <- function(post0, post_bar1, g_hist_i, g_path, s, dg, gamma, epsilon) {
  if (s == 1L) {
    g_plus <- g_path; g_plus[1] <- g_plus[1] + dg
    pb_base <- ce_reweight_posterior(post0, g_path[1], gamma)
    pb_plus <- ce_reweight_posterior(post0, g_plus[1], gamma)
    v_base <- fwd_researcher_value(post0, pb_base, g_hist_i, g_path, gamma, epsilon)
    v_plus <- fwd_researcher_value(post0, pb_plus, g_hist_i, g_plus, gamma, epsilon)
  } else {
    g_plus <- g_path; g_plus[s] <- g_plus[s] + dg
    v_base <- fwd_researcher_value(post0, post_bar1, g_hist_i, g_path, gamma, epsilon)
    v_plus <- fwd_researcher_value(post0, post_bar1, g_hist_i, g_plus, gamma, epsilon)
  }
  (v_plus - v_base) / dg
}

# Receding-horizon CE plan over H rounds for the current planning round.
# Returns an n x H matrix; the caller executes column 1 only.
#   H == 1 (last round): NO forward machinery — plain info-updated greedy on
#          the current posterior. Matches v5 closed-loop g2 and avoids the
#          ce_reweight RNG that would otherwise diverge from v5.
#   H >= 2: greedy over H*n (researcher, round) slots with the forward
#           marginals; global argmax, ties broken (round asc, researcher asc)
#           to match v5's "A over B, lowest index" rule.
plan_forward_ce <- function(post0_list, B_rem, delta, gamma, epsilon, H,
                            g_hist = list(), seed_round1 = NULL) {
  n <- length(post0_list)
  ghi <- function(i) vapply(g_hist, function(gh) gh[i], numeric(1))
  if (H == 1L) {
    # Single remaining round: plain info-updated greedy. A seed floor already
    # counts toward the budget, so greedy adds only B_rem − Σseed on top of it.
    seed_sum <- if (is.null(seed_round1)) 0 else sum(seed_round1)
    g <- greedy_round_t(post0_list, B_rem - seed_sum, delta, gamma, epsilon,
                        g_hist = g_hist, init_g = seed_round1)
    return(matrix(g, nrow = n, ncol = 1))
  }

  plan <- matrix(0, n, H)
  if (!is.null(seed_round1)) plan[, 1] <- seed_round1
  remaining <- B_rem - sum(plan)

  bar1 <- lapply(seq_len(n), function(i)
    ce_reweight_posterior(post0_list[[i]], plan[i, 1], gamma))

  mv <- matrix(0, n, H)
  for (i in seq_len(n))
    for (s in seq_len(H))
      mv[i, s] <- fwd_marginal(post0_list[[i]], bar1[[i]], ghi(i), plan[i, ], s,
                               delta, gamma, epsilon)

  while (remaining >= delta) {
    # global argmax with tie-break (smallest round s, then smallest i)
    best_val <- -Inf; best_i <- 1L; best_s <- 1L
    for (s in seq_len(H)) {
      i <- which.max(mv[, s])
      if (mv[i, s] > best_val) { best_val <- mv[i, s]; best_i <- i; best_s <- s }
    }
    plan[best_i, best_s] <- plan[best_i, best_s] + delta
    remaining <- remaining - delta
    if (best_s == 1L)
      bar1[[best_i]] <- ce_reweight_posterior(post0_list[[best_i]], plan[best_i, 1], gamma)
    for (s in seq_len(H))
      mv[best_i, s] <- fwd_marginal(post0_list[[best_i]], bar1[[best_i]], ghi(best_i),
                                    plan[best_i, ], s, delta, gamma, epsilon)
  }
  plan
}

# ============================================================
# Per-strategy allocation for the CURRENT round t
# ============================================================
# Returns the executed grant vector g_t (n) for round t under strategy S.
# tranche = B_total / T (per-round budget for the non-forward strategies).
# Forward strategies ignore `tranche` and plan the full remaining budget.
allocate_round <- function(
    S, t, T_rounds, posts, tranche, B_total, delta, gamma, epsilon,
    x_seed, p_prev, B_rem, g_hist
) {
  n <- length(p_prev)                              # posts is NULL for S1/S2/S3
  if (S == 1L) {                                   # No funding
    return(rep(0, n))
  } else if (S == 2L) {                            # Uniform every round
    return(rep(B_total / (n * T_rounds), n))
  } else if (S == 3L) {                            # Naive ∝ observed pubs
    return(allocate_naive(p_prev, tranche))
  } else if (S %in% c(4L, 5L)) {                   # Myopic (pubs [+ grant])
    return(greedy_round_t(posts, tranche, delta, gamma, epsilon, g_hist = g_hist))
  } else if (S == 6L) {                            # Myopic + round-1 seed floor
    if (t == 1L) {
      seed <- rep(x_seed * tranche / n, n)
      return(greedy_round_t(posts, (1 - x_seed) * tranche, delta, gamma, epsilon,
                            g_hist = g_hist, init_g = seed))
    }
    return(greedy_round_t(posts, tranche, delta, gamma, epsilon, g_hist = g_hist))
  } else if (S %in% c(7L, 8L)) {                   # Forward (pubs [+ grant])
    H <- T_rounds - t + 1L
    plan <- plan_forward_ce(posts, B_rem, delta, gamma, epsilon, H, g_hist = g_hist)
    return(plan[, 1])
  } else if (S == 9L) {                            # Forward + round-1 seed floor
    H <- T_rounds - t + 1L
    seed <- if (t == 1L) rep(x_seed * tranche / n, n) else NULL
    plan <- plan_forward_ce(posts, B_rem, delta, gamma, epsilon, H,
                            g_hist = g_hist, seed_round1 = seed)
    return(plan[, 1])
  }
  stop("unknown strategy ", S)
}

# whether a strategy uses each signal / needs a posterior
strategy_uses_grant   <- function(S) S %in% c(5L, 8L)
strategy_uses_posterior <- function(S) S %in% c(4L, 5L, 6L, 7L, 8L, 9L)
strategy_seed_frac    <- function(S, x_seed) if (S %in% c(6L, 9L)) x_seed else 0

# ============================================================
# Simulate one strategy over T rounds
# ============================================================
# Realizes rounds 1..T under the TRUE state, rebuilding the posterior from
# accumulated observations each round. Returns per-round grants, expected
# output (sum of true λ, the reported metric), and schedule diagnostics.
simulate_trial_T <- function(
    S, T_rounds, K_true, R0, p_init, sigma_r, sigma_k,
    B_total, delta, gamma, epsilon, x_seed,
    M, k_min, k_shape, r_min, r_shape, tau_r, tau_k,
    use_resource_signal, detail = FALSE
) {
  n       <- length(K_true)
  tranche <- B_total / T_rounds
  use_grant <- strategy_uses_grant(S)

  K_cur   <- K_true
  B_rem   <- B_total
  p_hist  <- list()          # realized pubs, rounds 1..t-1
  g_hist  <- list()          # executed grants, rounds 1..t-1
  g_rounds  <- vector("list", T_rounds)
  lam_rounds <- numeric(T_rounds)
  # detail-mode trajectories (populated from already-computed values only —
  # no extra RNG, so the T=2 anchor is preserved)
  if (detail) {
    K_rounds <- vector("list", T_rounds)   # knowledge ENTERING round t
    R_rounds <- vector("list", T_rounds)   # resources this round (R0 + g_t)
    lam_full <- vector("list", T_rounds)   # per-researcher expected output
  }

  for (t in seq_len(T_rounds)) {
    posts <- NULL
    if (strategy_uses_posterior(S)) {
      posts <- build_posteriors_hist(
        p_init, sigma_r, sigma_k, M, k_min, k_shape, r_min, r_shape, gamma,
        tau_r, tau_k, use_resource_signal, use_grant,
        p_hist = p_hist, g_hist = g_hist
      )
    }
    p_prev <- if (t == 1L) p_init else p_hist[[t - 1L]]

    g_t <- allocate_round(S, t, T_rounds, posts, tranche, B_total, delta,
                          gamma, epsilon, x_seed, p_prev, B_rem, g_hist)

    R_t   <- R0 + g_t
    lam_t <- lambda_rate(K_cur, R_t, gamma)          # TRUE expected output
    p_t   <- rpois(n, pmax(lam_t, LAMBDA_FLOOR))     # realization → next posterior

    if (detail) { K_rounds[[t]] <- K_cur; R_rounds[[t]] <- R_t; lam_full[[t]] <- lam_t }

    g_rounds[[t]]  <- g_t
    lam_rounds[t]  <- sum(lam_t)
    p_hist[[t]]    <- p_t
    g_hist[[t]]    <- g_t
    B_rem          <- B_rem - sum(g_t)
    K_cur          <- update_knowledge(K_cur, R_t, epsilon)   # compound for t+1
  }

  # ---- schedule diagnostics ----
  spend_t     <- vapply(g_rounds, sum, numeric(1))
  total_spend <- sum(spend_t)
  # Budget-share schedule α_t = (Σ g_t)/B_total and its center-of-mass b_idx
  # (spec §4). b_idx = 0.5 under an equal-tranche schedule.
  alpha_t <- if (B_total > 0) spend_t / B_total else rep(NA_real_, T_rounds)
  b_idx   <- if (B_total > 0) sum(seq_len(T_rounds) * alpha_t) / (T_rounds + 1) else NA_real_
  # Back-compat scalar α: round-1 share of ACTUAL SPEND (v5 definition; NA when
  # nothing is spent). Matches run_simulation()'s alpha for the T=2 anchor.
  alpha_v5 <- if (total_spend > 0) spend_t[1] / total_spend else NA_real_

  out <- list(
    strategy       = S,
    name           = STRATEGY_NAMES[S],
    g_rounds       = g_rounds,
    total_expected = sum(lam_rounds),
    total_output   = sum(unlist(p_hist)),   # realized pubs (always returned; app aggregation reads it)
    alpha          = alpha_v5,
    alpha_t        = alpha_t,
    b_idx          = b_idx,
    gini_g1        = gini(g_rounds[[1]])
  )
  if (detail) {
    out$K_rounds   <- K_rounds              # K entering each round (K_rounds[[1]] == K1, [[2]] == K2)
    out$R_rounds   <- R_rounds              # R0 + g_t each round
    out$lam_rounds <- lam_full              # per-researcher expected output each round
    out$p_rounds   <- p_hist                # realized pub counts each round
  }
  out
}

# Gini coefficient of a non-negative vector (0 if all equal / all zero).
gini <- function(x) {
  x <- sort(x); n <- length(x); s <- sum(x)
  if (s <= 0 || n < 2) return(0)
  sum((2 * seq_len(n) - n - 1) * x) / (n * s)
}

# ============================================================
# Top-level T-round simulation (mirrors run_simulation, all 9 strategies)
# ============================================================
run_simulation_T <- function(
    seed = 1, T_rounds = 2, n = 50,
    k_min = 1.0, k_shape = 2.0, r_min = 1.0, r_shape = 2.0,
    rho_kr = 0, gamma = 1.0, epsilon = 0.1, b = 0.5,
    n_steps = 50, tau_r = 1.0, tau_k = 1.0,
    use_resource_signal = TRUE, n_pre_rounds = 0,
    x_seed = 0.25, M = 200, strategies = 1:9, verbose = FALSE,
    detail = FALSE
) {
  set.seed(seed)

  E_R     <- r_min * r_shape / (r_shape - 1)     # Pareto mean (mean-E[R] norm.)
  B       <- b * n * E_R
  B_total <- 2 * B                               # total budget across T rounds
  delta   <- B / n_steps

  pop <- draw_initial_population(n, k_min, k_shape, r_min, r_shape, rho_kr)
  K0_init <- pop$K0; R0_init <- pop$R0
  rho_s   <- suppressWarnings(cor(K0_init, R0_init, method = "spearman"))

  lam_init <- lambda_rate(K0_init, R0_init, gamma)
  p_init   <- rpois(n, pmax(lam_init, LAMBDA_FLOOR))

  # Pre-rounds (naive), updating K and cumulative pubs — matches v5.
  K_current <- K0_init; R0_current <- R0_init; p_cumul <- p_init
  if (n_pre_rounds > 0) {
    for (r in seq_len(n_pre_rounds)) {
      g_pre   <- allocate_naive(p_cumul, B)
      R_round <- R0_current + g_pre
      lam_pre <- lambda_rate(K_current, R_round, gamma)
      p_round <- rpois(n, pmax(lam_pre, LAMBDA_FLOOR))
      p_cumul <- p_cumul + p_round
      K_current <- update_knowledge(K_current, R_round, epsilon)
    }
  }

  sigs    <- draw_signals(K_current, R0_current, tau_r, tau_k)
  sigma_r <- sigs$sigma_r; sigma_k <- sigs$sigma_k

  results <- vector("list", 9)
  for (S in strategies) {
    set.seed(seed * 1000)                        # common random numbers / strategy
    results[[S]] <- simulate_trial_T(
      S, T_rounds, K_current, R0_current, p_cumul, sigma_r, sigma_k,
      B_total, delta, gamma, epsilon, x_seed,
      M, k_min, k_shape, r_min, r_shape, tau_r, tau_k, use_resource_signal,
      detail = detail
    )
  }

  res <- list(
    params = list(seed = seed, T_rounds = T_rounds, n = n, b = b, B = B,
                  B_total = B_total, epsilon = epsilon, tau_k = tau_k,
                  tau_r = tau_r, k_shape = k_shape, r_shape = r_shape,
                  k_min = k_min, r_min = r_min, gamma = gamma,
                  rho_kr = rho_kr, x_seed = x_seed, M = M, n_steps = n_steps,
                  n_pre_rounds = n_pre_rounds, use_resource_signal = use_resource_signal),
    rho_s = rho_s,
    strategies = results
  )
  # detail mode: attach the population/state trajectories the Shiny app reads
  if (detail) {
    res$initial_state       <- list(K = K0_init, R = R0_init)
    res$post_preround_state <- list(K = K_current, R = R0_current)
    res$K_at_start          <- K_current
    res$R0_at_start         <- R0_current
    res$signals             <- list(sigma_r = sigma_r, sigma_k = sigma_k)
    res$p_init              <- p_init
    res$p_cumul             <- p_cumul
  }
  res
}
