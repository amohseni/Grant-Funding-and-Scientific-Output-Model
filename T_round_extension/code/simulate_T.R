# ============================================================
# simulate_T.R  —  T-round generalization of the v5 model
# ------------------------------------------------------------
# Extends app.R (v5, 2-round) to T ∈ {1,...,5} rounds. Prime
# directive: the T=2 path reproduces v5 exactly.
#
# DEPENDS ON app.R primitives (source app.R first, stripped of
# the trailing shinyApp() call): rpareto, draw_initial_population,
# lambda_rate, update_knowledge, draw_signals, loglik_*,
# ce_reweight_posterior, effective_sample_size, allocate_naive,
# STRATEGY_NAMES.
#
# Design of the generalization (each reduces to v5 at T=2):
#   * Posterior accumulates ALL past-round pub observations
#     (build_posteriors_hist). At t=2 == v5 build_posteriors(p1,g1).
#   * Myopic round-t marginal compounds K through executed history
#     (post_lambda_round_t). At t=1 == posterior_expected_lambda;
#     at t=2 == posterior_expected_lambda_r2.
#   * Forward = receding-horizon CE plan. Round 1 valued at π0,
#     rounds ≥2 at the CE-reweighted π̄1 (one-step info anticipation),
#     CE-compounded beyond. At H=2 the round-1 / round-2 marginals
#     are exactly v5's forward_A / forward_B; the last round (H=1)
#     is the plain info-updated greedy == v5 closed-loop g2.
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
    use_resource_signal
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

  list(
    strategy       = S,
    name           = STRATEGY_NAMES[S],
    g_rounds       = g_rounds,
    total_expected = sum(lam_rounds),
    alpha          = alpha_v5,
    alpha_t        = alpha_t,
    b_idx          = b_idx,
    gini_g1        = gini(g_rounds[[1]])
  )
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
    x_seed = 0.25, M = 200, strategies = 1:9, verbose = FALSE
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
      M, k_min, k_shape, r_min, r_shape, tau_r, tau_k, use_resource_signal
    )
  }

  list(
    params = list(seed = seed, T_rounds = T_rounds, n = n, b = b, B = B,
                  B_total = B_total, epsilon = epsilon, tau_k = tau_k,
                  tau_r = tau_r, k_shape = k_shape, r_shape = r_shape,
                  rho_kr = rho_kr, x_seed = x_seed, M = M, n_steps = n_steps),
    rho_s = rho_s,
    strategies = results
  )
}
