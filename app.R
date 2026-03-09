# app.R
# ============================================================
# Grant-funding model explorer (Shiny) — v2
# ------------------------------------------------------------
# A self-contained Shiny app for a 2-round discrete grant
# funding model with 9 allocation strategies, two noisy
# signals (resource + grant), bottleneck measures, and
# pre-round staging.
#
# Run:
#   shiny::runApp("app.R")
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
  # Draw correlated uniforms via Gaussian copula
  if (abs(rho_kr) < 1e-10) {
    u_k <- runif(n)
    u_r <- runif(n)
  } else {
    z1 <- rnorm(n)
    z2 <- rho_kr * z1 + sqrt(1 - rho_kr^2) * rnorm(n)
    u_k <- pnorm(z1)
    u_r <- pnorm(z2)
  }
  # Inverse CDF for Pareto: F^{-1}(u) = x_min / (1-u)^(1/shape)
  K0 <- k_min / (1 - u_k)^(1 / k_shape)
  R0 <- r_min / (1 - u_r)^(1 / r_shape)
  list(K0 = K0, R0 = R0)
}

# ----- Model primitives -----

# Output rate (harmonic-mean bottleneck)
lambda_rate <- function(K, R, gamma) {
  gamma * (K * R) / (K + R)
}

# Knowledge update (discrete, between rounds)
update_knowledge <- function(K, R, epsilon, kappa) {
  K + epsilon * K * (1 - K / (kappa * R))
}

# ----- Bottleneck measures -----

compute_bottleneck <- function(K, R, kappa) {
  D <- (K - R) / (K + R)
  S <- D^2
  Tght <- K / (kappa * R)
  list(D = D, S = S, Tght = Tght)
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

loglik_resource_signal <- function(sigma_r, R, tau_r) {
  dnorm(sigma_r, mean = R, sd = tau_r, log = TRUE)
}

loglik_grant_signal <- function(sigma_k, K, tau_k) {
  dnorm(sigma_k, mean = K, sd = tau_k, log = TRUE)
}

# ----- Posterior inference (importance sampling) -----

# Build posterior for a single researcher
posterior_samples_single <- function(
    p_obs, sigma_r, sigma_k,
    M, k_min, k_shape, r_min, r_shape, gamma,
    tau_r, tau_k,
    use_resource_signal, use_grant_signal
) {
  # Prior samples
  K_s <- rpareto(M, k_min, k_shape)
  R_s <- rpareto(M, r_min, r_shape)

  # Log-likelihood: publications

  ll <- loglik_pubs(p_obs, K_s, R_s, gamma)

  # Optionally add resource signal

  if (use_resource_signal) {
    ll <- ll + loglik_resource_signal(sigma_r, R_s, tau_r)
  }

  # Optionally add grant signal
  if (use_grant_signal) {
    ll <- ll + loglik_grant_signal(sigma_k, K_s, tau_k)
  }

  # Normalize weights
  ll <- ll - max(ll)
  w <- exp(ll)
  w <- w / sum(w)

  list(K0 = K_s, R0 = R_s, w = w)
}

# Build posteriors for all researchers
build_posteriors <- function(
    p_obs, sigma_r, sigma_k,
    M, k_min, k_shape, r_min, r_shape, gamma,
    tau_r, tau_k,
    use_resource_signal, use_grant_signal
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
      use_grant_signal = use_grant_signal
    )
  }
  posts
}

# ----- Expected output for one round given (K, R, gamma) -----

# Posterior expected lambda for researcher i given grant g
posterior_expected_lambda <- function(post, g, gamma) {
  R <- post$R0 + g
  lam <- lambda_rate(post$K0, R, gamma)
  sum(post$w * lam)
}

# Posterior marginal value of grant for researcher i
posterior_marginal_lambda <- function(post, g, dg, gamma) {
  (posterior_expected_lambda(post, g + dg, gamma) -
     posterior_expected_lambda(post, g, gamma)) / dg
}

# ----- Two-round expected output (for forward-looking strategies) -----

# Expected total output over both rounds for a researcher given grants g1, g2
# Uses posterior samples; returns expected (lambda1 + lambda2)
posterior_expected_two_round <- function(post, g1, g2, gamma, epsilon, kappa) {
  n_s <- length(post$K0)
  K1 <- post$K0
  R1 <- post$R0 + g1
  lam1 <- lambda_rate(K1, R1, gamma)

  # Knowledge update between rounds
  K2 <- update_knowledge(K1, R1, epsilon, kappa)
  R2 <- post$R0 + g2
  lam2 <- lambda_rate(K2, R2, gamma)

  sum(post$w * (lam1 + lam2))
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
  # Round 1: proportional to initial pubs
  g1 <- allocate_naive(p_init, B)
  # Simulate round 1 output (expected, not stochastic for strategy comparison)
  R1 <- R0 + g1
  lam1 <- lambda_rate(K_current, R1, gamma)
  # Round 2: proportional to round-1 output (use lam1 as proxy for pubs)
  g2 <- allocate_naive(lam1, B)
  list(g1 = g1, g2 = g2, lam1 = lam1)
}

# ----- Greedy allocation for one round (Bayesian strategies) -----
allocate_greedy_one_round <- function(posts, B, delta, gamma) {
  n <- length(posts)
  g <- rep(0, n)
  remaining <- B
  dg <- delta

  # Initial marginal values

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

# ----- Myopic Bayes (pubs only): Strategy 4 -----
# ----- Myopic Bayes (pubs + grant): Strategy 5 -----
# Both use allocate_greedy_one_round; differ only in posterior construction

run_myopic_bayes <- function(
    K_true, R0, p_init, sigma_r, sigma_k,
    B, delta, gamma, epsilon, kappa,
    M, k_min, k_shape, r_min, r_shape,
    tau_r, tau_k,
    use_resource_signal, use_grant_signal
) {
  n <- length(K_true)

  # -- Round 1 --
  posts1 <- build_posteriors(
    p_obs = p_init, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal = use_grant_signal
  )
  g1 <- allocate_greedy_one_round(posts1, B, delta, gamma)

  # Realize round 1 under TRUE state
  R1 <- R0 + g1
  lam1_true <- lambda_rate(K_true, R1, gamma)
  p1 <- rpois(n, pmax(lam1_true, 1e-12))

  # Knowledge update
  K2_true <- update_knowledge(K_true, R1, epsilon, kappa)

  # -- Round 2: update posteriors with round-1 pubs --
  # Cumulative pubs = p_init + p1
  p_cumul <- p_init + p1
  posts2 <- build_posteriors(
    p_obs = p_cumul, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal = use_grant_signal
  )
  g2 <- allocate_greedy_one_round(posts2, B, delta, gamma)

  R2 <- R0 + g2
  lam2_true <- lambda_rate(K2_true, R2, gamma)
  p2 <- rpois(n, pmax(lam2_true, 1e-12))

  list(
    g1 = g1, g2 = g2,
    p1 = p1, p2 = p2,
    lam1 = lam1_true, lam2 = lam2_true,
    K1 = K_true, K2 = K2_true,
    R1 = R1, R2 = R2,
    total_output = sum(p1 + p2),
    total_expected = sum(lam1_true + lam2_true)
  )
}

# ----- Myopic Bayes (pubs + seed): Strategy 6 -----
run_myopic_seed <- function(
    K_true, R0, p_init, sigma_r, sigma_k,
    B, delta, gamma, epsilon, kappa, x_seed,
    M, k_min, k_shape, r_min, r_shape,
    tau_r, tau_k, use_resource_signal
) {
  n <- length(K_true)

  # -- Round 1: x_seed fraction uniform, rest myopic --
  g1_uniform <- rep(x_seed * B / n, n)
  remaining_B <- (1 - x_seed) * B

  posts1 <- build_posteriors(
    p_obs = p_init, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal = FALSE
  )
  g1_opt <- allocate_greedy_one_round(posts1, remaining_B, delta, gamma)
  g1 <- g1_uniform + g1_opt

  # Realize round 1
  R1 <- R0 + g1
  lam1_true <- lambda_rate(K_true, R1, gamma)
  p1 <- rpois(n, pmax(lam1_true, 1e-12))
  K2_true <- update_knowledge(K_true, R1, epsilon, kappa)

  # -- Round 2: full myopic, no seed --
  p_cumul <- p_init + p1
  posts2 <- build_posteriors(
    p_obs = p_cumul, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal = FALSE
  )
  g2 <- allocate_greedy_one_round(posts2, B, delta, gamma)

  R2 <- R0 + g2
  lam2_true <- lambda_rate(K2_true, R2, gamma)
  p2 <- rpois(n, pmax(lam2_true, 1e-12))

  list(
    g1 = g1, g2 = g2,
    p1 = p1, p2 = p2,
    lam1 = lam1_true, lam2 = lam2_true,
    K1 = K_true, K2 = K2_true,
    R1 = R1, R2 = R2,
    total_output = sum(p1 + p2),
    total_expected = sum(lam1_true + lam2_true)
  )
}

# ----- Forward-looking Bayes: Strategies 7, 8, 9 -----
# Approximate via MC: for each candidate g1 allocation,
# simulate round-1 outcomes, compute optimal round-2 response, average.
# For tractability, use grid search over per-researcher g1 candidates.

run_forward_bayes <- function(
    K_true, R0, p_init, sigma_r, sigma_k,
    B, delta, gamma, epsilon, kappa,
    M, k_min, k_shape, r_min, r_shape,
    tau_r, tau_k,
    use_resource_signal, use_grant_signal,
    x_seed = 0, use_seed = FALSE,
    n_mc_forward = 20
) {
  n <- length(K_true)

  # Build round-1 posteriors
  posts1 <- build_posteriors(
    p_obs = p_init, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal = use_grant_signal
  )

  # For forward-looking: evaluate two-round expected value
  # We approximate by computing the greedy allocation that maximizes
  # E[lam1 + lam2] using the posterior

  # Marginal two-round value for researcher i given current g1
  mv_two_round <- function(posts, g1_vec, dg, gamma, epsilon, kappa, B_r2) {
    n <- length(posts)
    vapply(seq_len(n), function(i) {
      # Current two-round value
      # For round 2, assume greedy optimal given round-1 allocation
      # Approximate: give researcher i dg more in round 1, measure gain
      v_base <- posterior_expected_two_round(
        posts[[i]], g1_vec[i], 0, gamma, epsilon, kappa
      )
      v_plus <- posterior_expected_two_round(
        posts[[i]], g1_vec[i] + dg, 0, gamma, epsilon, kappa
      )
      (v_plus - v_base) / dg
    }, numeric(1))
  }

  # Greedy allocation for round 1 using two-round marginal values
  if (use_seed) {
    g1_uniform <- rep(x_seed * B / n, n)
    remaining_B1 <- (1 - x_seed) * B
  } else {
    g1_uniform <- rep(0, n)
    remaining_B1 <- B
  }

  g1_opt <- rep(0, n)
  remaining <- remaining_B1
  dg <- delta

  mv <- mv_two_round(posts1, g1_uniform + g1_opt, dg, gamma, epsilon, kappa, B)

  while (remaining >= delta) {
    i_star <- which.max(mv)
    g1_opt[i_star] <- g1_opt[i_star] + delta
    remaining <- remaining - delta
    # Recompute only for funded researcher
    v_base <- posterior_expected_two_round(
      posts1[[i_star]], g1_uniform[i_star] + g1_opt[i_star], 0, gamma, epsilon, kappa
    )
    v_plus <- posterior_expected_two_round(
      posts1[[i_star]], g1_uniform[i_star] + g1_opt[i_star] + dg, 0, gamma, epsilon, kappa
    )
    mv[i_star] <- (v_plus - v_base) / dg
  }

  g1 <- g1_uniform + g1_opt

  # Realize round 1 under TRUE state
  R1 <- R0 + g1
  lam1_true <- lambda_rate(K_true, R1, gamma)
  p1 <- rpois(n, pmax(lam1_true, 1e-12))
  K2_true <- update_knowledge(K_true, R1, epsilon, kappa)

  # -- Round 2: myopic optimal given updated posteriors --
  p_cumul <- p_init + p1
  posts2 <- build_posteriors(
    p_obs = p_cumul, sigma_r = sigma_r, sigma_k = sigma_k,
    M = M, k_min = k_min, k_shape = k_shape,
    r_min = r_min, r_shape = r_shape, gamma = gamma,
    tau_r = tau_r, tau_k = tau_k,
    use_resource_signal = use_resource_signal,
    use_grant_signal = use_grant_signal
  )

  if (use_seed) {
    g2 <- allocate_greedy_one_round(posts2, B, delta, gamma)
  } else {
    g2 <- allocate_greedy_one_round(posts2, B, delta, gamma)
  }

  R2 <- R0 + g2
  lam2_true <- lambda_rate(K2_true, R2, gamma)
  p2 <- rpois(n, pmax(lam2_true, 1e-12))

  list(
    g1 = g1, g2 = g2,
    p1 = p1, p2 = p2,
    lam1 = lam1_true, lam2 = lam2_true,
    K1 = K_true, K2 = K2_true,
    R1 = R1, R2 = R2,
    total_output = sum(p1 + p2),
    total_expected = sum(lam1_true + lam2_true)
  )
}

# ============================================================
# Main simulation pipeline
# ============================================================

run_simulation <- function(
    seed = 1,
    n = 100,
    # Knowledge distribution (power law)
    k_min = 1.0, k_shape = 2.0,
    # Resource distribution (power law)
    r_min = 1.0, r_shape = 2.0,
    # K-R correlation
    rho_kr = 0,
    # Dynamics
    gamma = 1.0,
    epsilon = 0.1,
    kappa = 1.0,
    # Funding
    B = 50,
    delta = 1,
    # Signal parameters
    tau_r = 1.0,
    tau_k = 1.0,
    use_resource_signal = TRUE,
    # Pre-rounds
    n_pre_rounds = 0,
    # Seed strategy
    x_seed = 0.5,
    # Posterior
    M = 1500,
    # Which strategies to run
    strategies = 1:9,
    verbose = FALSE
) {
  set.seed(seed)

  # ---- Step 1: Draw initial conditions ----
  pop <- draw_initial_population(n, k_min, k_shape, r_min, r_shape, rho_kr)
  K0_init <- pop$K0
  R0_init <- pop$R0

  # Record initial distribution (before pre-rounds)
  initial_state <- list(K = K0_init, R = R0_init)

  # ---- Step 2: Draw signals (once, before any rounds) ----
  sigs <- draw_signals(K0_init, R0_init, tau_r, tau_k)
  sigma_r <- sigs$sigma_r
  sigma_k <- sigs$sigma_k

  # ---- Step 3: Generate initial publication counts ----
  lam_init <- lambda_rate(K0_init, R0_init, gamma)
  p_init <- rpois(n, pmax(lam_init, 1e-12))

  # ---- Step 4: Run pre-rounds (Naive allocation) ----
  K_current <- K0_init
  R0_current <- R0_init
  p_cumul <- p_init

  if (n_pre_rounds > 0) {
    for (r in seq_len(n_pre_rounds)) {
      # Allocate via Naive
      g_pre <- allocate_naive(p_cumul, B)
      # Update resources for this round
      R_round <- R0_current + g_pre
      # Draw output
      lam_pre <- lambda_rate(K_current, R_round, gamma)
      p_round <- rpois(n, pmax(lam_pre, 1e-12))
      p_cumul <- p_cumul + p_round
      # Update knowledge
      K_current <- update_knowledge(K_current, R_round, epsilon, kappa)
      # Baseline resources persist (grants are per-round)
    }
  }

  # Record post-pre-round state
  post_preround_state <- list(K = K_current, R = R0_current)

  # ---- Step 5: Run main 2 rounds for each selected strategy ----
  results <- list()

  for (s in strategies) {
    # Use same RNG state offset for fairness
    set.seed(seed * 1000 + s)

    if (s == 1) {
      # No funding
      alloc <- allocate_no_funding(n, B)
      R1 <- R0_current + alloc$g1
      lam1 <- lambda_rate(K_current, R1, gamma)
      p1 <- rpois(n, pmax(lam1, 1e-12))
      K2 <- update_knowledge(K_current, R1, epsilon, kappa)
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
        total_output = sum(p1 + p2),
        total_expected = sum(lam1 + lam2)
      )

    } else if (s == 2) {
      # Uniform seed
      alloc <- allocate_uniform(n, B)
      R1 <- R0_current + alloc$g1
      lam1 <- lambda_rate(K_current, R1, gamma)
      p1 <- rpois(n, pmax(lam1, 1e-12))
      K2 <- update_knowledge(K_current, R1, epsilon, kappa)
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
        total_output = sum(p1 + p2),
        total_expected = sum(lam1 + lam2)
      )

    } else if (s == 3) {
      # Naive
      nv <- allocate_naive_two_round(n, B, p_cumul, K_current, R0_current, NULL, gamma)
      g1 <- nv$g1
      R1 <- R0_current + g1
      lam1 <- lambda_rate(K_current, R1, gamma)
      p1 <- rpois(n, pmax(lam1, 1e-12))
      K2 <- update_knowledge(K_current, R1, epsilon, kappa)
      # Round 2 naive uses round-1 pubs
      g2 <- allocate_naive(p_cumul + p1, B)
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
        total_output = sum(p1 + p2),
        total_expected = sum(lam1 + lam2)
      )

    } else if (s == 4) {
      # Myopic Bayes (pubs only)
      res_s <- run_myopic_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon, kappa,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal = FALSE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s

    } else if (s == 5) {
      # Myopic Bayes (pubs + grant)
      res_s <- run_myopic_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon, kappa,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal = TRUE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s

    } else if (s == 6) {
      # Myopic Bayes (pubs + seed)
      res_s <- run_myopic_seed(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon, kappa, x_seed,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k, use_resource_signal
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s

    } else if (s == 7) {
      # Forward Bayes (pubs only)
      res_s <- run_forward_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon, kappa,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal = FALSE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s

    } else if (s == 8) {
      # Forward Bayes (pubs + grant)
      res_s <- run_forward_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon, kappa,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal = TRUE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s

    } else if (s == 9) {
      # Forward Bayes (pubs + seed)
      res_s <- run_forward_bayes(
        K_current, R0_current, p_cumul, sigma_r, sigma_k,
        B, delta, gamma, epsilon, kappa,
        M, k_min, k_shape, r_min, r_shape,
        tau_r, tau_k,
        use_resource_signal = use_resource_signal,
        use_grant_signal = FALSE,
        x_seed = x_seed, use_seed = TRUE
      )
      res_s$name <- STRATEGY_NAMES[s]
      results[[s]] <- res_s
    }
  }

  # ---- Compute bottleneck measures ----
  bn_initial <- compute_bottleneck(initial_state$K, initial_state$R, kappa)
  bn_post_pre <- compute_bottleneck(post_preround_state$K, post_preround_state$R, kappa)

  list(
    params = list(
      seed = seed, n = n,
      k_min = k_min, k_shape = k_shape,
      r_min = r_min, r_shape = r_shape,
      rho_kr = rho_kr,
      gamma = gamma, epsilon = epsilon, kappa = kappa,
      B = B, delta = delta,
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
      .primary-controls {
        background-color: #fff;
        padding: 12px 0;
        margin-bottom: 8px;
      }
      .primary-controls .form-group { margin-bottom: 14px; }
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
           strategies, and forward-looking strategies that account for
           knowledge growth between rounds.")
      )
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      # ==== Primary Parameters (always visible) ====
      div(class = "primary-controls",
        sliderInput("n", "Number of researchers", min = 10, max = 500, value = 100, step = 10),
        div(class = "param-help", "Population size."),

        sliderInput("B", "Per-round budget (B)", min = 1, max = 200, value = 50, step = 5),
        div(class = "param-help", "Total funding available each round."),

        sliderInput("gamma", "Output scaling (\\(\\gamma\\))", min = 0.1, max = 5, value = 1.0, step = 0.1),
        div(class = "param-help", "Scales the paper production rate \\(\\lambda = \\gamma K R / (K+R)\\)."),

        checkboxGroupInput("strategies", "Strategies to compare",
          choices = setNames(1:9, STRATEGY_NAMES),
          selected = 1:9
        ),

        numericInput("seed", "Random seed", value = 1, min = 1, step = 1),
        div(class = "param-help", "Controls random draws; change to see different populations.")
      ),

      actionButton("run", "Run Simulation", class = "btn-primary btn-run"),

      actionButton("reset", "Reset to Defaults", class = "btn-default btn-run"),

      hr(),

      # ==== Knowledge Distribution ====
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

          sliderInput("epsilon", "Growth rate (\\(\\varepsilon\\))", min = 0.01, max = 1, value = 0.1, step = 0.01),
          div(class = "param-help", "Rate of knowledge growth between rounds."),

          sliderInput("kappa", "Knowledge ceiling / resource (\\(\\kappa\\))", min = 0.1, max = 5, value = 1.0, step = 0.1),
          div(class = "param-help", "Knowledge ceiling = \\(\\kappa \\times R\\).")
        )
      ),

      # ==== Resource Distribution ====
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

      # ==== Initial Conditions ====
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

      # ==== Signal Parameters ====
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

      # ==== Seed Strategy Parameters ====
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
        verbatimTextOutput("summary", placeholder = TRUE)
      ),

      tabsetPanel(id = "main_tabs",

        # ==== Tab 1: Strategy Comparison ====
        tabPanel("Strategy Comparison",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Total expected research output
                  (sum of \\(\\lambda\\) over both rounds and all researchers) for
                  each selected strategy. Higher bars = more total output.")
          ),
          plotOutput("fig_strategy_comparison", height = 520)
        ),

        # ==== Tab 2: Pre-Round Effects ====
        tabPanel("Pre-Round Effects",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> How pre-rounds of naive
                  (publication-proportional) funding reshape the population
                  <em>before</em> the main strategies begin. <strong>Left:</strong>
                  initial \\((K, R)\\) as drawn from the power law distributions.
                  <strong>Right:</strong> \\((K, R)\\) after all pre-rounds complete.
                  Pre-rounds let correlation emerge organically through the
                  funding&ndash;knowledge feedback loop. Points colored by bottleneck
                  direction \\(D_i = (K-R)/(K+R)\\): blue = knowledge-bottlenecked,
                  red = resource-bottlenecked. Pearson correlation annotated.")
          ),
          plotOutput("fig_preround_effects", height = 480)
        ),

        # ==== Tab 3: Funding Effects ====
        tabPanel("Funding Effects",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> How the two main funding rounds
                  reshape the \\((K, R)\\) distribution for a selected strategy.
                  <strong>Left:</strong> state entering round 1 (post-pre-rounds).
                  <strong>Center:</strong> effective \\((K, R)\\) during round 1
                  (\\(R\\) includes round-1 grant). <strong>Right:</strong> effective
                  \\((K, R)\\) during round 2 (knowledge has grown; \\(R\\) includes
                  round-2 grant). Arrows on the right panel show each researcher's
                  trajectory from start to round 2. Compare how different strategies
                  move researchers in this space.")
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

        # ==== Tab 4: Bottleneck Analysis ====
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

        # ==== Tab 5: Value of Signals ====
        tabPanel("Value of Signals",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Compares strategies that differ
                  only in their information set to quantify the value of the grant
                  signal and resource signal.")
          ),
          plotOutput("fig_signal_value", height = 480)
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

  # ---- Reset to defaults ----
  observeEvent(input$reset, {
    updateSliderInput(session, "n", value = 100)
    updateSliderInput(session, "B", value = 50)
    updateSliderInput(session, "gamma", value = 1.0)
    updateCheckboxGroupInput(session, "strategies", selected = 1:9)
    updateNumericInput(session, "seed", value = 1)
    updateSliderInput(session, "k_shape", value = 2.0)
    updateSliderInput(session, "k_min", value = 1.0)
    updateSliderInput(session, "epsilon", value = 0.1)
    updateSliderInput(session, "kappa", value = 1.0)
    updateSliderInput(session, "r_shape", value = 2.0)
    updateSliderInput(session, "r_min", value = 1.0)
    updateSliderInput(session, "rho_kr", value = 0)
    updateSliderInput(session, "n_pre_rounds", value = 0)
    updateCheckboxInput(session, "use_resource_signal", value = TRUE)
    updateSliderInput(session, "tau_r", value = 1.0)
    updateSliderInput(session, "tau_k", value = 1.0)
    updateSliderInput(session, "x_seed", value = 0.5)
  })

  # ---- Main simulation reactive ----
  sim_result <- reactiveVal(NULL)

  observe({
    input$run
    isolate({
      strats <- as.integer(input$strategies)
      if (length(strats) == 0) strats <- 1

      withProgress(message = "Running simulation...", value = 0.1, {
        res <- run_simulation(
          seed = input$seed,
          n = input$n,
          k_min = input$k_min, k_shape = input$k_shape,
          r_min = input$r_min, r_shape = input$r_shape,
          rho_kr = input$rho_kr,
          gamma = input$gamma,
          epsilon = input$epsilon,
          kappa = input$kappa,
          B = input$B,
          delta = 1,
          tau_r = input$tau_r,
          tau_k = input$tau_k,
          use_resource_signal = input$use_resource_signal,
          n_pre_rounds = input$n_pre_rounds,
          x_seed = input$x_seed,
          M = 1500,
          strategies = strats
        )
        incProgress(0.9)
        sim_result(res)
      })
    })
  })

  # ---- Summary output ----
  output$summary <- renderPrint({
    res <- sim_result()
    req(res)

    cat("==================================================================\n")
    cat("                         KEY RESULTS\n")
    cat("==================================================================\n\n")
    cat(sprintf("Population: n = %d | Budget/round: B = %.0f | Pre-rounds: %d\n",
                res$params$n, res$params$B, res$params$n_pre_rounds))
    cat(sprintf("K ~ Pareto(min=%.1f, shape=%.1f) | R ~ Pareto(min=%.1f, shape=%.1f) | rho = %.1f\n\n",
                res$params$k_min, res$params$k_shape,
                res$params$r_min, res$params$r_shape,
                res$params$rho_kr))

    cat("TOTAL EXPECTED OUTPUT (lambda_1 + lambda_2 summed over researchers)\n")
    cat("------------------------------------------------------------------\n")

    for (s in seq_along(res$strategies)) {
      r <- res$strategies[[s]]
      if (is.null(r)) next
      cat(sprintf("  %-30s %7.1f\n", r$name, r$total_expected))
    }
    cat("\n")

    # Correlation info
    cat(sprintf("K-R Pearson correlation: initial = %.3f, post-pre-rounds = %.3f\n",
                cor(res$initial_state$K, res$initial_state$R),
                cor(res$post_preround_state$K, res$post_preround_state$R)))
  })

  # ---- Tab 1: Strategy Comparison ----
  output$fig_strategy_comparison <- renderPlot({
    res <- sim_result()
    req(res)

    df <- tibble(
      Strategy = character(),
      Total = numeric()
    )
    for (s in seq_along(res$strategies)) {
      r <- res$strategies[[s]]
      if (is.null(r)) next
      df <- bind_rows(df, tibble(Strategy = r$name, Total = r$total_expected))
    }

    df$Strategy <- factor(df$Strategy, levels = df$Strategy)

    baseline_val <- df$Total[df$Strategy == "No funding"]
    if (length(baseline_val) == 0) baseline_val <- 0

    ggplot(df, aes(x = Strategy, y = Total, fill = Strategy)) +
      geom_col(width = 0.7, alpha = 0.9) +
      {if (baseline_val > 0) geom_hline(yintercept = baseline_val, linetype = "dashed", alpha = 0.5)} +
      geom_text(aes(label = sprintf("%.1f", Total)),
                vjust = -0.5, size = 3.5, fontface = "bold") +
      scale_fill_manual(values = STRATEGY_COLORS, drop = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        title = "Total expected output by strategy (2 rounds)",
        subtitle = sprintf("n = %d, B = %.0f/round, gamma = %.1f",
                           res$params$n, res$params$B, res$params$gamma),
        x = NULL,
        y = "Total expected output (sum of lambdas)"
      ) +
      theme_sim(base_size = 13) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 9, angle = 30, hjust = 1),
        plot.subtitle = element_text(color = "grey40", size = 11)
      )
  }, res = 110)

  # ---- Tab 2: Pre-Round Effects ----
  output$fig_preround_effects <- renderPlot({
    res <- sim_result()
    req(res)

    # Initial state (raw draws)
    df_init <- tibble(
      K = res$initial_state$K,
      R = res$initial_state$R,
      D = (res$initial_state$K - res$initial_state$R) /
          (res$initial_state$K + res$initial_state$R),
      Phase = "Initial (raw draws)"
    )
    cor_init <- cor(res$initial_state$K, res$initial_state$R)

    # Post-pre-round state
    n_pre <- res$params$n_pre_rounds
    phase_label <- if (n_pre > 0) {
      sprintf("After %d pre-round%s", n_pre, if (n_pre > 1) "s" else "")
    } else {
      "After 0 pre-rounds (unchanged)"
    }

    df_post <- tibble(
      K = res$post_preround_state$K,
      R = res$post_preround_state$R,
      D = (res$post_preround_state$K - res$post_preround_state$R) /
          (res$post_preround_state$K + res$post_preround_state$R),
      Phase = phase_label
    )
    cor_post <- cor(res$post_preround_state$K, res$post_preround_state$R)

    df_all <- bind_rows(df_init, df_post)
    df_all$Phase <- factor(df_all$Phase, levels = c("Initial (raw draws)", phase_label))

    # Common axis limits across both panels
    max_val <- max(c(df_all$K, df_all$R), na.rm = TRUE) * 1.05

    # Correlation annotations
    ann <- tibble(
      Phase = factor(c("Initial (raw draws)", phase_label),
                     levels = levels(df_all$Phase)),
      label = c(sprintf("r = %.3f", cor_init), sprintf("r = %.3f", cor_post))
    )

    ggplot(df_all, aes(x = R, y = K, color = D)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) +
      geom_point(size = 1.8, alpha = 0.8) +
      scale_color_gradient2(
        low = "#2166ac", mid = "grey80", high = "#b2182b",
        midpoint = 0, name = "Bottleneck\ndirection (D)",
        limits = c(-1, 1)
      ) +
      geom_text(data = ann, aes(x = max_val * 0.05, y = max_val * 0.95, label = label),
                inherit.aes = FALSE, hjust = 0, size = 4, fontface = "bold") +
      facet_wrap(~ Phase) +
      coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
      labs(
        title = "Pre-round effects: how naive funding history reshapes the population",
        subtitle = sprintf("K ~ Pareto(%.1f, %.1f), R ~ Pareto(%.1f, %.1f), rho = %.1f",
                           res$params$k_min, res$params$k_shape,
                           res$params$r_min, res$params$r_shape,
                           res$params$rho_kr),
        x = "Resources (R)",
        y = "Knowledge (K)"
      ) +
      theme_sim(base_size = 12) +
      theme(plot.subtitle = element_text(color = "grey40", size = 10))
  }, res = 110)

  # ---- Tab 3: Funding Effects ----
  output$fig_funding_effects <- renderPlot({
    res <- sim_result()
    req(res)

    s_idx <- as.integer(input$fe_strategy)
    strat <- res$strategies[[s_idx]]
    req(strat)

    n <- res$params$n

    # Three snapshots:
    # 1. Start of main rounds (post-pre-round K, baseline R)
    # 2. During round 1 (K1, R1 = R0 + g1)
    # 3. During round 2 (K2, R2 = R0 + g2) — K has grown from round-1 funding

    df_start <- tibble(
      K = res$K_at_start,
      R = res$R0_at_start,
      D = (res$K_at_start - res$R0_at_start) /
          (res$K_at_start + res$R0_at_start),
      id = seq_len(n),
      Stage = "Before funding"
    )

    df_r1 <- tibble(
      K = strat$K1,
      R = strat$R1,
      D = (strat$K1 - strat$R1) / (strat$K1 + strat$R1),
      id = seq_len(n),
      Stage = "Round 1"
    )

    df_r2 <- tibble(
      K = strat$K2,
      R = strat$R2,
      D = (strat$K2 - strat$R2) / (strat$K2 + strat$R2),
      id = seq_len(n),
      Stage = "Round 2"
    )

    df_all <- bind_rows(df_start, df_r1, df_r2)
    df_all$Stage <- factor(df_all$Stage,
                           levels = c("Before funding", "Round 1", "Round 2"))

    # Common axis limits
    max_k <- max(df_all$K, na.rm = TRUE) * 1.05
    max_r <- max(df_all$R, na.rm = TRUE) * 1.05
    max_val <- max(max_k, max_r)

    # Correlation annotations per panel
    cor_start <- cor(res$K_at_start, res$R0_at_start)
    cor_r1 <- cor(strat$K1, strat$R1)
    cor_r2 <- cor(strat$K2, strat$R2)
    ann <- tibble(
      Stage = factor(c("Before funding", "Round 1", "Round 2"),
                     levels = levels(df_all$Stage)),
      label = c(sprintf("r = %.3f", cor_start),
                sprintf("r = %.3f", cor_r1),
                sprintf("r = %.3f", cor_r2))
    )

    # Build arrow segments: start → round 2 (the net movement)
    df_arrows <- tibble(
      x_start = res$R0_at_start,
      y_start = res$K_at_start,
      x_end = strat$R2,
      y_end = strat$K2,
      # Only show arrows for researchers who actually moved
      moved = abs(strat$R2 - res$R0_at_start) > 0.01 | abs(strat$K2 - res$K_at_start) > 0.01,
      Stage = factor("Round 2", levels = levels(df_all$Stage))
    ) %>% filter(moved)

    ggplot(df_all, aes(x = R, y = K, color = D)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) +
      # Arrows only on the Round 2 panel
      geom_segment(data = df_arrows,
                   aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                   inherit.aes = FALSE,
                   arrow = arrow(length = unit(0.08, "inches"), type = "closed"),
                   color = "grey50", alpha = 0.25, linewidth = 0.3) +
      geom_point(size = 1.8, alpha = 0.8) +
      scale_color_gradient2(
        low = "#2166ac", mid = "grey80", high = "#b2182b",
        midpoint = 0, name = "Bottleneck\ndirection (D)",
        limits = c(-1, 1)
      ) +
      geom_text(data = ann,
                aes(x = max_val * 0.05, y = max_val * 0.95, label = label),
                inherit.aes = FALSE, hjust = 0, size = 3.8, fontface = "bold") +
      facet_wrap(~ Stage) +
      coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
      labs(
        title = sprintf("Funding effects on (K, R) distribution: %s", strat$name),
        subtitle = sprintf("B = %.0f/round | Arrows (right panel) show start-to-round-2 movement",
                           res$params$B),
        x = "Resources (R)",
        y = "Knowledge (K)"
      ) +
      theme_sim(base_size = 12) +
      theme(plot.subtitle = element_text(color = "grey40", size = 10))
  }, res = 110)

  # ---- Tab 4: Bottleneck Analysis ----
  output$fig_bottleneck <- renderPlot({
    res <- sim_result()
    req(res)

    s_idx <- as.integer(input$bn_strategy)
    strat <- res$strategies[[s_idx]]
    req(strat)

    kappa_val <- res$params$kappa

    # Initial (post-pre-round)
    bn0 <- compute_bottleneck(res$K_at_start, res$R0_at_start, kappa_val)
    # After round 1
    bn1 <- compute_bottleneck(strat$K1, strat$R1, kappa_val)
    # After round 2
    bn2 <- compute_bottleneck(strat$K2, strat$R2, kappa_val)

    df_bn <- bind_rows(
      tibble(D = bn0$D, S = bn0$S, Tght = bn0$Tght, Stage = "Initial"),
      tibble(D = bn1$D, S = bn1$S, Tght = bn1$Tght, Stage = "After Round 1"),
      tibble(D = bn2$D, S = bn2$S, Tght = bn2$Tght, Stage = "After Round 2")
    )
    df_bn$Stage <- factor(df_bn$Stage, levels = c("Initial", "After Round 1", "After Round 2"))

    # Three-panel plot: D, S, T
    df_d <- df_bn %>% select(Stage, value = D) %>% mutate(Measure = "Direction (D)")
    df_s <- df_bn %>% select(Stage, value = S) %>% mutate(Measure = "Severity (S)")
    df_t <- df_bn %>% select(Stage, value = Tght) %>% mutate(Measure = "Ceiling tightness (T)")

    df_long <- bind_rows(df_d, df_s, df_t)
    df_long$Measure <- factor(df_long$Measure,
                              levels = c("Direction (D)", "Severity (S)", "Ceiling tightness (T)"))

    ggplot(df_long, aes(x = value, fill = Stage)) +
      geom_density(alpha = 0.4, linewidth = 0.4) +
      facet_wrap(~ Measure, scales = "free") +
      scale_fill_manual(values = c("Initial" = "#bdbdbd", "After Round 1" = "#42a5f5", "After Round 2" = "#1565c0")) +
      labs(
        title = sprintf("Bottleneck measures: %s", strat$name),
        x = "Value",
        y = "Density",
        fill = "Stage"
      ) +
      theme_sim(base_size = 12) +
      theme(legend.position = "bottom")
  }, res = 110)

  # ---- Tab 5: Value of Signals ----
  output$fig_signal_value <- renderPlot({
    res <- sim_result()
    req(res)

    # Build comparison pairs
    comparisons <- tibble(
      Comparison = character(),
      Without = numeric(),
      With = numeric(),
      Gain = numeric()
    )

    # Myopic: pubs vs pubs+grant
    s4 <- res$strategies[[4]]
    s5 <- res$strategies[[5]]
    if (!is.null(s4) && !is.null(s5)) {
      comparisons <- bind_rows(comparisons, tibble(
        Comparison = "Myopic:\ngrant signal",
        Without = s4$total_expected,
        With = s5$total_expected,
        Gain = s5$total_expected - s4$total_expected
      ))
    }

    # Forward: pubs vs pubs+grant
    s7 <- res$strategies[[7]]
    s8 <- res$strategies[[8]]
    if (!is.null(s7) && !is.null(s8)) {
      comparisons <- bind_rows(comparisons, tibble(
        Comparison = "Forward:\ngrant signal",
        Without = s7$total_expected,
        With = s8$total_expected,
        Gain = s8$total_expected - s7$total_expected
      ))
    }

    if (nrow(comparisons) == 0) {
      plot.new()
      text(0.5, 0.5, "Run strategies 4, 5, 7, 8 to see signal value comparisons",
           cex = 1.2)
      return(invisible(NULL))
    }

    comparisons$Comparison <- factor(comparisons$Comparison,
                                     levels = comparisons$Comparison)

    df_long <- comparisons %>%
      pivot_longer(cols = c(Without, With), names_to = "Info", values_to = "Output")
    df_long$Info <- factor(df_long$Info, levels = c("Without", "With"))

    ggplot(df_long, aes(x = Comparison, y = Output, fill = Info)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
      geom_text(data = comparisons,
                aes(x = Comparison, y = With + max(comparisons$With) * 0.02,
                    label = sprintf("+%.1f", Gain)),
                inherit.aes = FALSE, size = 4, fontface = "bold", color = "#c62828") +
      scale_fill_manual(values = c("Without" = "#90caf9", "With" = "#1565c0"),
                        labels = c("Without" = "Without signal", "With" = "With signal")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
      labs(
        title = "Value of additional signals",
        subtitle = "Gain from adding grant signal to pubs-only strategies",
        x = NULL,
        y = "Total expected output",
        fill = NULL
      ) +
      theme_sim(base_size = 13) +
      theme(
        legend.position = "bottom",
        plot.subtitle = element_text(color = "grey40", size = 11)
      )
  }, res = 110)
}

shinyApp(ui, server)
