# app.R
# ============================================================
# Grant-funding model explorer (Shiny) — v5 (CE-forward consolidated)
# ------------------------------------------------------------
# Changes from v4:
#
#   The CE machinery introduced in v4 as a separate Strategy 10
#   ("Optimal Bayes") is now folded directly into the Forward
#   planner. All three forward strategies (7, 8, 9) inherit it.
#
#   Rationale: "Forward" is intended to denote the maximally
#   forward-looking Bayesian planner. The v3/v4 forward planner
#   omitted Forces D and E (information value of round-1 grants),
#   which made it not truly forward-looking. Keeping a separate
#   "Optimal" strategy implicitly conceded that Forward was
#   suboptimal — an awkward taxonomy.
#
#   After v5, the three Forward variants form a clean orthogonal
#   ablation along (signal, intervention):
#     S7: pubs only,         no seed
#     S8: pubs + grant,      no seed
#     S9: pubs only,         seed  (no grant signal)
#
#   The planner is CE-approximate: imagined round-1 observations
#   are replaced by their prior-predictive mean rather than
#   Monte-Carlo-integrated. CE introduces an upward bias on α
#   (Jensen). An exact MC reference is a separate validation
#   script, not a strategy.
#
#   Carryovers from v4 (still in force):
#     - Closed-loop execution: round 2 is re-greedied under the
#       actual posterior π1 built from the realized p1.
#     - Harmonic-mean knowledge dynamics. κ fully absent.
#     - Free inter-round budget split. α emerges from optimization.
# ============================================================

options(shiny.sanitize.errors = FALSE)

suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(scales)
})

# ============================================================
# CONSTANTS
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

STRATEGY_COLORS <- c(
  "No funding"            = "#bdbdbd",
  "Uniform seed"          = "#a5d6a7",
  "Naive (prop. to pubs)" = "#90caf9",
  "Myopic (pubs)"         = "#42a5f5",
  "Myopic (pubs + grant)" = "#1e88e5",
  "Myopic (pubs + seed)"  = "#7e57c2",
  "Forward (pubs)"        = "#ef5350",
  "Forward (pubs + grant)"= "#c62828",
  "Forward (pubs + seed)" = "#e65100"
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
# Keeps the same (K, R) sample atoms; only weights change.
# If ESS falls below ess_floor (default M/4), perform an SIR resampling
# step: draw M atoms with replacement using the current weights, reset
# weights to 1/M, then re-apply the CE Poisson likelihood. This refreshes
# the sample pool without leaving the existing posterior support.
ce_reweight_posterior <- function(post0, g1_i, gamma, ess_floor = NULL) {
  M <- length(post0$w)
  if (is.null(ess_floor)) ess_floor <- M / 4
  
  K_s <- post0$K0
  R_s <- post0$R0
  w_s <- post0$w
  
  # Imagined observation: prior-predictive mean of round-1 publications
  lam_at_g <- lambda_rate(K_s, R_s + g1_i, gamma)
  p_bar <- sum(w_s * lam_at_g)
  
  # Reweight by Poisson(p_bar | λ)
  log_lik <- loglik_pubs_continuous(p_bar, K_s, R_s + g1_i, gamma)
  log_w <- log(pmax(w_s, 1e-300)) + log_lik
  log_w <- log_w - max(log_w)
  w_new <- exp(log_w); w_new <- w_new / sum(w_new)
  
  if (effective_sample_size(w_new) < ess_floor) {
    # SIR resample from π0 using current π0 weights, then re-apply CE
    idx <- sample.int(M, size = M, replace = TRUE, prob = w_s)
    K_s <- K_s[idx]; R_s <- R_s[idx]
    lam_at_g <- lambda_rate(K_s, R_s + g1_i, gamma)
    p_bar <- mean(lam_at_g)                        # uniform weights post-SIR
    log_lik <- loglik_pubs_continuous(p_bar, K_s, R_s + g1_i, gamma)
    log_w <- log_lik - max(log_lik)
    w_new <- exp(log_w); w_new <- w_new / sum(w_new)
  }
  
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
# Plotting utilities
# ============================================================

theme_sim <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = base_size + 1, margin = margin(b = 6)),
      axis.text = element_text(color = "grey20"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "plain"),
      plot.margin = margin(12, 12, 12, 12)
    )
}

# Helper for "Value of X" comparison bar plots (signals, forward looking, seed).
# Expects a tibble with columns: Comparison (factor), Without, With, Gain.
# Pull se_total_expected from a strategy result; return 0 if absent / NA.
se_or_zero <- function(strategy_result) {
  v <- strategy_result$se_total_expected
  if (is.null(v) || is.na(v)) 0 else v
}

make_value_plot <- function(comparisons, title, subtitle,
                            baseline_label, treatment_label,
                            baseline_color, treatment_color) {
  # SE columns optional; default to 0 (no error bars)
  if (is.null(comparisons$SE_Without)) comparisons$SE_Without <- 0
  if (is.null(comparisons$SE_With)) comparisons$SE_With <- 0
  
  df_long <- comparisons %>%
    pivot_longer(cols = c(Without, With), names_to = "Info", values_to = "Output") %>%
    mutate(SE = ifelse(Info == "Without", SE_Without, SE_With))
  df_long$Info <- factor(df_long$Info, levels = c("Without", "With"))
  
  y_max <- max(c(comparisons$With, comparisons$Without))
  has_errors <- any(df_long$SE > 0)
  
  p <- ggplot(df_long, aes(x = Comparison, y = Output, fill = Info)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9)
  
  if (has_errors) {
    p <- p + geom_errorbar(aes(ymin = Output - SE, ymax = Output + SE),
                           position = position_dodge(width = 0.7),
                           width = 0.2, alpha = 0.7, linewidth = 0.5)
  }
  
  p +
    geom_text(data = comparisons,
              aes(x = Comparison,
                  y = pmax(With, Without) + y_max * 0.04,
                  label = sprintf("%+.1f", Gain)),
              inherit.aes = FALSE, size = 4, fontface = "bold", color = "#c62828") +
    scale_fill_manual(values = c("Without" = baseline_color, "With" = treatment_color),
                      labels = c("Without" = baseline_label, "With" = treatment_label)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = title, subtitle = subtitle, x = NULL,
         y = "Total expected output", fill = NULL) +
    theme_sim(base_size = 13) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(color = "grey40", size = 11))
}

# ============================================================
# Shiny UI
# ============================================================

ui <- fluidPage(
  withMathJax(),
  
  tags$head(
    tags$style(HTML("
      .control-label { font-weight: 600; }
      .shiny-output-error-validation { color: #b00020; font-weight: 600; }
      .description-box {
        background-color: #f8f9fa;
        border-left: 4px solid #495057;
        padding: 16px 20px;
        margin-bottom: 20px;
        font-size: 14px;
        line-height: 1.6;
      }
      .description-box p { margin: 0 0 10px 0; }
      .description-box p:last-child { margin-bottom: 0; }
      .panel-group { margin-bottom: 12px; }
      .panel { border: 1px solid #ddd; border-radius: 4px; margin-bottom: 8px; }
      .panel-heading {
        background-color: #f5f5f5;
        padding: 10px 14px;
        cursor: pointer;
        border-radius: 3px;
      }
      .panel-heading:hover { background-color: #e8e8e8; }
      .panel-title {
        font-size: 13px;
        font-weight: 600;
        margin: 0;
        color: #333;
      }
      .panel-title .caret-icon {
        float: right;
        color: #666;
        transition: transform 0.2s;
      }
      .panel-body {
        padding: 12px 14px;
        border-top: 1px solid #ddd;
        background-color: #fafafa;
      }
      .param-help {
        font-size: 11px;
        color: #666;
        margin-top: -8px;
        margin-bottom: 12px;
        line-height: 1.4;
      }
      .plot-explanation {
        background-color: #f0f4f8;
        border-left: 3px solid #4a90d9;
        padding: 10px 14px;
        margin-bottom: 12px;
        font-size: 13px;
        line-height: 1.5;
      }
      .btn-run {
        width: 100%;
        margin-top: 8px;
        margin-bottom: 8px;
      }
      .summary-box {
        font-size: 13px;
        background-color: #f8f9fa;
        border: 1px solid #e0e0e0;
        border-radius: 4px;
        padding: 12px;
        margin-bottom: 12px;
      }
      .summary-title {
        margin: 0 0 6px 0;
        font-size: 15px;
        font-weight: 600;
        color: #2c3e50;
      }
      .summary-box pre {
        margin: 0;
        padding: 6px 0 0 0;
        background: transparent;
        border: none;
      }
    "))
  ),
  
  titlePanel("A Model of Grant Funding and Scientific Output"),
  
  fluidRow(
    column(12,
           div(class = "description-box",
               p("This model investigates how funders should allocate grants across
           two discrete funding rounds when researchers differ in knowledge
           (\\(K\\)) and resources (\\(R\\)). Output is bottlenecked by whichever
           input is scarcer. We compare 9 allocation strategies under varying
           information regimes."),
               p("The funder's challenge: a researcher with high \\(K\\) and low \\(R\\)
           can look identical in publication count to one with low \\(K\\) and
           high \\(R\\). We compare non-Bayesian baselines, myopic Bayesian
           strategies, and forward-looking strategies that plan the full
           inter-round allocation to account for knowledge growth.")
           )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      actionButton("run", "Run Simulation", class = "btn-primary btn-run"),
      actionButton("reset", "Reset to Defaults", class = "btn-default btn-run"),
      
      numericInput("seed", "Random seed", value = 1, min = 1, max = 1000000, step = 1),
      div(class = "param-help", "Auto-advances to a new random value (1–10^6) after each run; set manually to reproduce a specific draw."),
      
      numericInput("n_trials", "Trials per run", value = 1, min = 1, max = 30, step = 1),
      div(class = "param-help", "Number of trials to average. Higher values reduce noise but increase compute time. Per-strategy plots use the first seed's state; comparison plots show means with error bars."),
      
      hr(),
      
      tags$div(class = "panel",
               tags$div(class = "panel-heading",
                        `data-toggle` = "collapse", `data-target` = "#panel-strategies",
                        tags$h4(class = "panel-title",
                                "Strategies to compare",
                                tags$span(class = "caret-icon", HTML("&#9662;"))
                        )
               ),
               tags$div(id = "panel-strategies", class = "collapse panel-body",
                        checkboxGroupInput("strategies", NULL,
                                           choices = setNames(1:9, STRATEGY_NAMES),
                                           selected = 1:9
                        )
               )
      ),
      
      tags$div(class = "panel",
               tags$div(class = "panel-heading",
                        `data-toggle` = "collapse", `data-target` = "#panel-popbudget",
                        tags$h4(class = "panel-title",
                                "Population & Budget",
                                tags$span(class = "caret-icon", HTML("&#9662;"))
                        )
               ),
               tags$div(id = "panel-popbudget", class = "collapse panel-body",
                        sliderInput("n", "Number of researchers", min = 10, max = 500, value = 100, step = 10),
                        div(class = "param-help", "Population size."),
                        
                        sliderInput("b", "Budget fraction (b = B / n·E[R])", min = 0.0, max = 2.0, value = 0.5, step = 0.05),
                        div(class = "param-help", "Per-round budget as a multiple of aggregate expected resource: b = 1 means total budget equals the population's total expected resource endowment."),
                        
                        sliderInput("gamma", "Output scaling (\\(\\gamma\\))", min = 0.1, max = 5, value = 1.0, step = 0.1),
                        div(class = "param-help", "Scales the paper production rate \\(\\lambda = \\gamma K R / (K+R)\\).")
               )
      ),
      
      tags$div(class = "panel",
               tags$div(class = "panel-heading",
                        `data-toggle` = "collapse", `data-target` = "#panel-knowledge",
                        tags$h4(class = "panel-title",
                                "Knowledge Distribution",
                                tags$span(class = "caret-icon", HTML("&#9662;"))
                        )
               ),
               tags$div(id = "panel-knowledge", class = "collapse panel-body",
                        sliderInput("k_shape", "Pareto shape (\\(\\alpha_K\\))", min = 1.1, max = 5, value = 2.0, step = 0.1),
                        div(class = "param-help", "Shape of the power law for knowledge; lower = heavier tail."),
                        
                        sliderInput("k_min", "Pareto scale (\\(k_{min}\\))", min = 0.1, max = 5, value = 1.0, step = 0.1),
                        div(class = "param-help", "Minimum knowledge value."),
                        
                        sliderInput("epsilon", "Growth rate (\\(\\varepsilon\\))", min = 0.01, max = 1, value = 0.3, step = 0.01),
                        div(class = "param-help", "Rate of knowledge growth between rounds. K grows in proportion to current output rate: \\(K_2 = K_1 + \\varepsilon \\cdot K_1 R_1 / (K_1 + R_1)\\).")
               )
      ),
      
      tags$div(class = "panel",
               tags$div(class = "panel-heading",
                        `data-toggle` = "collapse", `data-target` = "#panel-resources",
                        tags$h4(class = "panel-title",
                                "Resource Distribution",
                                tags$span(class = "caret-icon", HTML("&#9662;"))
                        )
               ),
               tags$div(id = "panel-resources", class = "collapse panel-body",
                        sliderInput("r_shape", "Pareto shape (\\(\\alpha_R\\))", min = 1.1, max = 5, value = 2.0, step = 0.1),
                        div(class = "param-help", "Shape of the power law for resources; lower = heavier tail."),
                        
                        sliderInput("r_min", "Pareto scale (\\(r_{min}\\))", min = 0.1, max = 5, value = 1.0, step = 0.1),
                        div(class = "param-help", "Minimum resource value.")
               )
      ),
      
      tags$div(class = "panel",
               tags$div(class = "panel-heading",
                        `data-toggle` = "collapse", `data-target` = "#panel-initial",
                        tags$h4(class = "panel-title",
                                "Initial Conditions",
                                tags$span(class = "caret-icon", HTML("&#9662;"))
                        )
               ),
               tags$div(id = "panel-initial", class = "collapse panel-body",
                        sliderInput("rho_kr", "K-R correlation (\\(\\rho_{KR}\\))", min = -0.9, max = 0.9, value = 0, step = 0.1),
                        div(class = "param-help", "Correlation between initial K and R. 0 = independent."),
                        
                        sliderInput("n_pre_rounds", "Pre-rounds", min = 0, max = 20, value = 0, step = 1),
                        div(class = "param-help", "Rounds of Naive funding before main strategies begin.")
               )
      ),
      
      tags$div(class = "panel",
               tags$div(class = "panel-heading",
                        `data-toggle` = "collapse", `data-target` = "#panel-signals",
                        tags$h4(class = "panel-title",
                                "Signal Parameters",
                                tags$span(class = "caret-icon", HTML("&#9662;"))
                        )
               ),
               tags$div(id = "panel-signals", class = "collapse panel-body",
                        checkboxInput("use_resource_signal", "Resource signal enabled", value = TRUE),
                        div(class = "param-help", "When enabled, Bayesian strategies observe a noisy signal of R."),
                        
                        sliderInput("tau_r", "Resource signal noise (\\(\\tau_R\\))", min = 0.01, max = 10, value = 1.0, step = 0.1),
                        div(class = "param-help", "SD of resource signal noise. Lower = more informative."),
                        
                        sliderInput("tau_k", "Grant signal noise (\\(\\tau_K\\))", min = 0.01, max = 10, value = 1.0, step = 0.1),
                        div(class = "param-help", "SD of grant proposal signal noise. Lower = more informative.")
               )
      ),
      
      tags$div(class = "panel",
               tags$div(class = "panel-heading",
                        `data-toggle` = "collapse", `data-target` = "#panel-seed-strat",
                        tags$h4(class = "panel-title",
                                "Seed Strategy Parameters",
                                tags$span(class = "caret-icon", HTML("&#9662;"))
                        )
               ),
               tags$div(id = "panel-seed-strat", class = "collapse panel-body",
                        sliderInput("x_seed", "Seed fraction (x)", min = 0, max = 1, value = 0.5, step = 0.05),
                        div(class = "param-help", "Fraction of round-1 budget allocated uniformly in seed strategies.")
               )
      ),
      
      hr(),
      helpText("Tip: reduce n or M if computation is slow.")
    ),
    
    mainPanel(
      width = 9,
      
      div(class = "summary-box",
          tags$h4("Simulation results", class = "summary-title"),
          verbatimTextOutput("summary", placeholder = TRUE)
      ),
      
      tabsetPanel(id = "main_tabs",
                  
                  tabPanel("Strategy Comparison",
                           div(class = "plot-explanation",
                               HTML("<strong>What this shows:</strong> Total expected research output
                  (sum of \\(\\lambda\\) over both rounds and all researchers) for
                  each selected strategy. Higher bars = more total output.")
                           ),
                           plotOutput("fig_strategy_comparison", height = 520)
                  ),
                  
                  tabPanel("Funding Effects",
                           div(class = "plot-explanation",
                               HTML("<strong>What this shows:</strong> How the two main funding rounds
                  reshape the \\((K, R)\\) distribution for a selected strategy.
                  <strong>Left:</strong> state entering round 1 (post-pre-rounds).
                  <strong>Center:</strong> effective \\((K, R)\\) during round 1
                  (\\(R\\) includes round-1 grant). <strong>Right:</strong> effective
                  \\((K, R)\\) during round 2 (knowledge has grown; \\(R\\) includes
                  round-2 grant). Arrows on the right panel show each researcher's
                  trajectory from start to round 2.")
                           ),
                           fluidRow(
                             column(4,
                                    selectInput("fe_strategy", "Strategy",
                                                choices = setNames(1:9, STRATEGY_NAMES),
                                                selected = 4
                                    )
                             )
                           ),
                           plotOutput("fig_funding_effects", height = 500)
                  ),
                  
                  tabPanel("Value of Grant Signals",
                           div(class = "plot-explanation",
                               HTML("<strong>What this shows:</strong> Compares strategies that differ
                  only in their information set to quantify the value of the grant
                  signal.")
                           ),
                           plotOutput("fig_signal_value", height = 480)
                  ),
                  
                  tabPanel("Value of Seed Grants",
                           div(class = "plot-explanation",
                               HTML("<strong>What this shows:</strong> Compares no-seed and uniform-seed
                  variants of the pubs-only strategies to isolate the value of the
                  seeding intervention. Two matched pairs: 4↔6 (myopic) and 7↔9
                  (forward). Both sides of each pair use only the publication
                  signal, so the comparison is clean.")
                           ),
                           plotOutput("fig_seed_value", height = 480)
                  ),
                  
                  tabPanel("Value of Forward Looking",
                           div(class = "plot-explanation",
                               HTML("<strong>What this shows:</strong> Compares myopic and forward
                  (CE) planners at matched information / intervention settings to
                  isolate the value of optimizing over two rounds rather than one.
                  Three matched pairs: 4↔7 (pubs only), 5↔8 (pubs + grant),
                  6↔9 (pubs + seed).")
                           ),
                           plotOutput("fig_forward_value", height = 480)
                  ),
                  
                  tabPanel("Bottleneck Analysis",
                           div(class = "plot-explanation",
                               HTML("<strong>What this shows:</strong> Distribution of bottleneck
                  measures across researchers at initial conditions, after round 1,
                  and after round 2. Select a strategy to view.")
                           ),
                           fluidRow(
                             column(4,
                                    selectInput("bn_strategy", "Strategy",
                                                choices = setNames(1:9, STRATEGY_NAMES),
                                                selected = 4
                                    )
                             )
                           ),
                           plotOutput("fig_bottleneck", height = 520)
                  ),
                  
                  tabPanel("Pre-Round Effects",
                           div(class = "plot-explanation",
                               HTML("<strong>What this shows:</strong> How pre-rounds of naive
                  (publication-proportional) funding reshape the population
                  <em>before</em> the main strategies begin. <strong>Left:</strong>
                  initial \\((K, R)\\) as drawn from the power law distributions.
                  <strong>Right:</strong> \\((K, R)\\) after all pre-rounds complete.
                  Points colored by bottleneck direction \\(D_i = (K-R)/(K+R)\\):
                  blue = knowledge-bottlenecked, red = resource-bottlenecked.
                  Pearson correlation annotated.")
                           ),
                           plotOutput("fig_preround_effects", height = 480)
                  )
      )
    )
  ),
  
  tags$script(HTML("
    $(document).on('click', '.panel-heading', function() {
      $(this).next('.collapse').slideToggle(200);
      $(this).find('.caret-icon').toggleClass('collapsed');
    });
  "))
)

# ============================================================
# Shiny Server
# ============================================================

server <- function(input, output, session) {
  
  observeEvent(input$reset, {
    updateSliderInput(session, "n", value = 100)
    updateSliderInput(session, "b", value = 0.5)
    updateSliderInput(session, "gamma", value = 1.0)
    updateCheckboxGroupInput(session, "strategies", selected = 1:9)
    updateNumericInput(session, "seed", value = 1)
    updateNumericInput(session, "n_trials", value = 1)
    updateSliderInput(session, "k_shape", value = 2.0)
    updateSliderInput(session, "k_min", value = 1.0)
    updateSliderInput(session, "epsilon", value = 0.3)
    updateSliderInput(session, "r_shape", value = 2.0)
    updateSliderInput(session, "r_min", value = 1.0)
    updateSliderInput(session, "rho_kr", value = 0)
    updateSliderInput(session, "n_pre_rounds", value = 0)
    updateCheckboxInput(session, "use_resource_signal", value = TRUE)
    updateSliderInput(session, "tau_r", value = 1.0)
    updateSliderInput(session, "tau_k", value = 1.0)
    updateSliderInput(session, "x_seed", value = 0.5)
  })
  
  sim_result <- reactiveVal(NULL)
  
  observe({
    input$run
    isolate({
      strats <- as.integer(input$strategies)
      if (length(strats) == 0) strats <- 1
      
      n_runs <- max(1L, as.integer(input$n_trials))
      base_seed <- as.integer(input$seed)
      
      withProgress(message = "Running simulation...", value = 0.0, {
        all_runs <- vector("list", n_runs)
        for (i in seq_len(n_runs)) {
          all_runs[[i]] <- run_simulation(
            seed = base_seed + (i - 1L),
            n = input$n,
            k_min = input$k_min, k_shape = input$k_shape,
            r_min = input$r_min, r_shape = input$r_shape,
            rho_kr = input$rho_kr,
            gamma = input$gamma,
            epsilon = input$epsilon,
            b = input$b,
            n_steps = 50,
            tau_r = input$tau_r,
            tau_k = input$tau_k,
            use_resource_signal = input$use_resource_signal,
            n_pre_rounds = input$n_pre_rounds,
            x_seed = input$x_seed,
            M = 1500,
            strategies = strats
          )
          setProgress(value = i / n_runs,
                      detail = sprintf("seed %d of %d", i, n_runs))
        }
        
        # Use first run as the representative for state-dependent plots
        # (Funding Effects, Bottleneck, Pre-Round) and as the base for aggregation
        res <- all_runs[[1]]
        res$n_trials <- n_runs
        res$seed_range <- c(base_seed, base_seed + n_runs - 1L)
        
        # Aggregate per-strategy metrics across runs
        if (n_runs > 1L) {
          for (s in seq_along(res$strategies)) {
            if (is.null(res$strategies[[s]])) next
            outs <- vapply(all_runs, function(r) {
              v <- r$strategies[[s]]$total_output;     if (is.null(v)) NA_real_ else v
            }, numeric(1))
            exps <- vapply(all_runs, function(r) {
              v <- r$strategies[[s]]$total_expected;   if (is.null(v)) NA_real_ else v
            }, numeric(1))
            alphas <- vapply(all_runs, function(r) {
              v <- r$strategies[[s]]$alpha;            if (is.null(v)) NA_real_ else v
            }, numeric(1))
            
            res$strategies[[s]]$total_output     <- mean(outs, na.rm = TRUE)
            res$strategies[[s]]$se_total_output  <- sd(outs,   na.rm = TRUE) / sqrt(sum(!is.na(outs)))
            res$strategies[[s]]$total_expected   <- mean(exps, na.rm = TRUE)
            res$strategies[[s]]$se_total_expected<- sd(exps,   na.rm = TRUE) / sqrt(sum(!is.na(exps)))
            res$strategies[[s]]$alpha            <- mean(alphas, na.rm = TRUE)
            res$strategies[[s]]$se_alpha         <- sd(alphas, na.rm = TRUE) / sqrt(sum(!is.na(alphas)))
          }
        }
        
        sim_result(res)
      })
      
      # Auto-advance seed to a fresh random value for the next run.
      # Advance past the range we just used to avoid overlap.
      updateNumericInput(session, "seed", value = sample.int(1000000, 1))
    })
  })
  
  output$summary <- renderPrint({
    res <- sim_result()
    req(res)
    
    multi <- !is.null(res$n_trials) && res$n_trials > 1
    
    if (multi) {
      cat(sprintf("Trials: %d (seeds %d–%d) | n = %d | b = %.2f (B = %.1f/round) | Pre-rounds: %d\n",
                  re