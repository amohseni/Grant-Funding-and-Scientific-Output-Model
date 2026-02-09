# ============================================================
# Minimal Grant-Funding Model + Standard Plots
# ------------------------------------------------------------
#
# Contents:
#   - A fast, deterministic expected-value simulator.
#   - Bayesian posterior approximation by importance sampling.
#   - Greedy funding allocator using posterior marginal values.
#   - Key plots of simulation outputs using ggplot2.
#
# Notes:
#   - No packages are required for simulation/inference; plots require:
#       ggplot2, dplyr, tidyr, purrr, scales
#     (Optional) viridisLite for nicer continuous scales.
#   - Knowledge dynamics use closed-form logistic growth (no decay path).
# ============================================================

# -----------------------------
# Packages (plots only)
# -----------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(scales)
})

set.seed(1)

# ============================================================
# 1) Distributions (no packages)
# ============================================================

rpareto <- function(n, k_min, shape) {
  # Pareto Type I on [k_min, ∞):  P(K > k) = (k_min/k)^shape
  u <- runif(n)
  k_min / (u)^(1 / shape)
}

rtruncnorm_pos <- function(n, mu, sd) {
  # Truncated Normal on (0, ∞) via rejection sampling.
  out <- numeric(0)
  while (length(out) < n) {
    x <- rnorm(n, mu, sd)
    out <- c(out, x[x > 0])
  }
  out[1:n]
}

# ============================================================
# 2) Model primitives + observation model
# ============================================================

lambda_rate <- function(K, R, alpha_out) {
  # Smooth bottleneck: approximately alpha_out * min(K, R) for K far from R
  alpha_out * (K * R) / (K + R)
}

loglik_obs <- function(p_obs, s_obs, K0, R0, tau, alpha_out, xi) {
  # Vectorized log-likelihood for observed (p_obs, s_obs) given (K0, R0).
  lam0 <- lambda_rate(K0, R0, alpha_out)
  ll_p <- dpois(p_obs, lambda = tau * lam0, log = TRUE)
  ll_s <- if (is.finite(xi)) dnorm(s_obs, mean = K0, sd = xi, log = TRUE) else 0
  ll_p + ll_s
}

posterior_samples <- function(p_obs, s_obs, M,
                              k_min, k_shape, mu_R, sd_R,
                              tau, alpha_out, xi) {
  # Importance sampling approximation to posterior over (K0, R0).
  K0 <- rpareto(M, k_min, k_shape)
  R0 <- rtruncnorm_pos(M, mu_R, sd_R)
  
  ll <- loglik_obs(p_obs, s_obs, K0, R0, tau, alpha_out, xi)
  
  # Stable normalization of weights
  ll <- ll - max(ll)
  w <- exp(ll)
  w <- w / sum(w)
  
  list(K0 = K0, R0 = R0, w = w)
}

# ============================================================
# 3) Fast expected value engine (no MC)
#    Knowledge growth: logistic with cap proportional to resources
# ============================================================

K_logistic_closed_form <- function(t, K0, cap, eps) {
  # Closed-form logistic solution for:
  #   dK/dt = eps*K*(1 - K/cap), with cap>0
  #   K(t) = cap / (1 + (cap/K0 - 1)*exp(-eps t))
  out <- numeric(length(K0))
  idx <- K0 > 0
  out[idx] <- cap[idx] / (1 + (cap[idx] / K0[idx] - 1) * exp(-eps * t))
  out
}

expected_papers_EV_no_decay <- function(K0, R, T, dt, eps, kappa, alpha_out) {
  # Expected papers over horizon T:
  #   E[P] = ∫_0^T lambda(K(t), R) dt
  # Knowledge path uses closed-form logistic; integral via trapezoid rule.
  stopifnot(length(K0) == length(R))
  stopifnot(T > 0, dt > 0, eps >= 0, kappa > 0, alpha_out >= 0)
  
  times <- seq(0, T, by = dt)
  if (tail(times, 1) < T) times <- c(times, T)
  
  cap <- kappa * R
  
  lam_prev <- lambda_rate(K_logistic_closed_form(times[1], K0, cap, eps), R, alpha_out)
  integral <- rep(0, length(K0))
  
  for (j in 2:length(times)) {
    t <- times[j]
    dt_j <- times[j] - times[j - 1]
    K_t <- K_logistic_closed_form(t, K0, cap, eps)
    lam_next <- lambda_rate(K_t, R, alpha_out)
    integral <- integral + 0.5 * (lam_prev + lam_next) * dt_j
    lam_prev <- lam_next
  }
  
  integral
}

posterior_expected_value <- function(post, g, T, dt, eps, kappa, alpha_out) {
  # E[∫ lambda dt | posterior] given grant g (resources increase linearly: R = R0 + g)
  R <- post$R0 + g
  V <- expected_papers_EV_no_decay(post$K0, R, T, dt, eps, kappa, alpha_out)
  sum(post$w * V)
}

posterior_marginal_value <- function(post, g, dg, T, dt, eps, kappa, alpha_out) {
  # Finite difference marginal value: (V(g+dg)-V(g))/dg
  (posterior_expected_value(post, g + dg, T, dt, eps, kappa, alpha_out) -
     posterior_expected_value(post, g,      T, dt, eps, kappa, alpha_out)) / dg
}

# ============================================================
# 4) Greedy funding allocation
# ============================================================

allocate_greedy <- function(posts, B, delta, T, dt, eps, kappa, alpha_out,
                            dg = delta, verbose = TRUE) {
  n <- length(posts)
  stopifnot(B >= 0, delta > 0, dg > 0)
  
  g <- rep(0, n)
  remaining <- B
  step <- 0
  
  mv <- vapply(seq_len(n), function(i) {
    posterior_marginal_value(posts[[i]], g[i], dg, T, dt, eps, kappa, alpha_out)
  }, numeric(1))
  
  while (remaining >= delta) {
    step <- step + 1
    i_star <- which.max(mv)
    
    g[i_star] <- g[i_star] + delta
    remaining <- remaining - delta
    
    # Separable objective: only recompute for the funded researcher
    mv[i_star] <- posterior_marginal_value(posts[[i_star]], g[i_star], dg, T, dt, eps, kappa, alpha_out)
    
    if (verbose && step %% 10 == 0) {
      cat(sprintf("step=%d  remaining=%.2f  funded=%d  g[%d]=%.2f  mv=%.6f\n",
                  step, remaining, i_star, i_star, g[i_star], mv[i_star]))
    }
  }
  
  list(g = g, remaining = remaining, mv = mv)
}

# ============================================================
# 5) One-run simulation (synthetic data) returning EVERYTHING
# ============================================================

run_demo_fast <- function(
    n = 40,
    # priors
    k_min = 1, k_shape = 2.2,
    mu_R = 10, sd_R = 3,
    # dynamics + output
    alpha_out = 0.12,
    eps = 0.25,
    kappa = 2.5,
    # horizon and integration step
    T = 24,
    dt = 0.25,
    # observation model
    tau = 6,
    xi = 3.0,          # set to Inf to ignore proposal signals
    # funding
    B = 60,
    delta = 1,
    # posterior approximation
    M = 1500,
    verbose = TRUE
) {
  # --- True population
  K0_true <- rpareto(n, k_min, k_shape)
  R0_true <- rtruncnorm_pos(n, mu_R, sd_R)
  
  # --- Observed signals
  lam0_true <- lambda_rate(K0_true, R0_true, alpha_out)
  p_obs <- rpois(n, tau * lam0_true)
  s_obs <- if (is.finite(xi)) rnorm(n, mean = K0_true, sd = xi) else rep(0, n)
  
  # --- Posteriors
  posts <- vector("list", n)
  for (i in seq_len(n)) {
    posts[[i]] <- posterior_samples(
      p_obs = p_obs[i], s_obs = s_obs[i], M = M,
      k_min = k_min, k_shape = k_shape, mu_R = mu_R, sd_R = sd_R,
      tau = tau, alpha_out = alpha_out, xi = xi
    )
  }
  
  # --- Allocate funds
  alloc <- allocate_greedy(posts, B, delta, T, dt, eps, kappa, alpha_out, dg = delta, verbose = verbose)
  
  # --- Evaluate expected outcomes under TRUE states (baseline vs funded)
  baseline <- sum(expected_papers_EV_no_decay(K0_true, R0_true,           T, dt, eps, kappa, alpha_out))
  funded   <- sum(expected_papers_EV_no_decay(K0_true, R0_true + alloc$g, T, dt, eps, kappa, alpha_out))
  
  list(
    params = list(
      n = n, B = B, delta = delta,
      k_min = k_min, k_shape = k_shape, mu_R = mu_R, sd_R = sd_R,
      alpha_out = alpha_out, eps = eps, kappa = kappa,
      T = T, dt = dt,
      tau = tau, xi = xi,
      M = M
    ),
    truth = list(K0 = K0_true, R0 = R0_true),
    obs = list(p = p_obs, s = s_obs),
    posts = posts,
    allocation = alloc,
    eval = list(
      expected_total_baseline = baseline,
      expected_total_funded = funded,
      expected_gain = funded - baseline
    )
  )
}

# ============================================================
# 6) Plotting: standardized aesthetics + Figures 1–7
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

scale_color_continuous_safe <- function(name = waiver(), ...) {
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    scale_color_viridis_c(name = name, option = "C", ...)
  } else {
    scale_color_gradient(name = name, low = "grey70", high = "grey10", ...)
  }
}

scale_fill_continuous_safe <- function(name = waiver(), ...) {
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    scale_fill_viridis_c(name = name, option = "C", ...)
  } else {
    scale_fill_gradient(name = name, low = "grey70", high = "grey10", ...)
  }
}

build_researcher_df <- function(res) {
  tibble(
    id = seq_along(res$truth$K0),
    K0 = res$truth$K0,
    R0 = res$truth$R0,
    p_obs = res$obs$p,
    g = res$allocation$g
  ) %>%
    mutate(
      funded = g > 0,
      R_funded = R0 + g
    )
}

wquantile <- function(x, w, probs = c(0.05, 0.95)) {
  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  approx(cw, x, xout = probs, ties = "ordered")$y
}

make_plots <- function(res, base_size = 12, contour_bins = 12, exemplars_n = 3) {
  df <- build_researcher_df(res)
  par <- res$params
  
  alpha_out <- par$alpha_out
  eps <- par$eps
  kappa <- par$kappa
  T <- par$T
  dt <- par$dt
  
  # ---------- Figure 1: (K,R) contours ----------
  R_grid <- seq(min(df$R0), max(df$R_funded), length.out = 160)
  K_grid <- seq(min(df$K0), max(df$K0), length.out = 160)
  grid <- expand_grid(K = K_grid, R = R_grid) %>% mutate(lambda = lambda_rate(K, R, alpha_out))
  
  fig1 <- ggplot() +
    geom_contour(data = grid, aes(x = R, y = K, z = lambda), bins = contour_bins, linewidth = 0.3, alpha = 0.85) +
    geom_point(data = df, aes(x = R0, y = K0, shape = funded), size = 2.1, alpha = 0.9) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19)) +
    labs(
      title = "Figure 1. Bottleneck geometry in (K,R)-space",
      x = "Baseline resources R0",
      y = "Initial knowledge K0",
      shape = "Funded"
    ) +
    theme_sim(base_size)
  
  cap1 <- paste(
    "Contours show the paper-production rate λ(K,R). Points show researchers at (K0,R0).",
    "Funding is valuable when it increases resources in the resource-limited region (high K relative to R)."
  )
  
  # ---------- Figure 2: p vs K (confounding) ----------
  fig2 <- ggplot(df, aes(x = K0, y = p_obs, color = R0, shape = funded)) +
    geom_point(size = 2.1, alpha = 0.9) +
    scale_color_continuous_safe("Baseline resources R0") +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19)) +
    labs(
      title = "Figure 2. Publications confound knowledge and resources",
      x = "True initial knowledge K0",
      y = "Observed past publications p",
      shape = "Funded"
    ) +
    theme_sim(base_size)
  
  cap2 <- paste(
    "Publication counts depend on both knowledge K and resources R. For fixed p, researchers can differ substantially in K",
    "because low resources suppress output. This is why ‘fund the most published’ is generally suboptimal."
  )
  
  # ---------- Figure 3: posterior vs truth (diagnostic) ----------
  post_summ <- map_dfr(seq_len(nrow(df)), function(i) {
    post <- res$posts[[i]]
    K_s <- post$K0
    w <- post$w
    m <- sum(w * K_s)
    qs <- wquantile(K_s, w, probs = c(0.05, 0.95))
    tibble(id = i, K_post_mean = m, K_post_lo = qs[1], K_post_hi = qs[2])
  })
  df3 <- df %>% left_join(post_summ, by = "id")
  
  fig3 <- ggplot(df3, aes(x = K0, y = K_post_mean, color = funded)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.35, alpha = 0.6) +
    geom_errorbar(aes(ymin = K_post_lo, ymax = K_post_hi), width = 0, alpha = 0.55) +
    geom_point(size = 2.0, alpha = 0.9) +
    scale_color_manual(values = c(`FALSE` = "grey35", `TRUE` = "black")) +
    labs(
      title = "Figure 3. Posterior recovery of knowledge (calibration check)",
      x = "True initial knowledge K0",
      y = "Posterior mean E[K0 | signals]",
      color = "Funded"
    ) +
    theme_sim(base_size)
  
  cap3 <- paste(
    "Posterior mean E[K0|signals] vs true K0, with 90% intervals. This separates inference quality from allocation logic:",
    "weak performance can reflect uninformative signals rather than a flawed allocator."
  )
  
  # ---------- Figure 4: marginal value curves for exemplars ----------
  ids_topfund <- df %>% arrange(desc(g)) %>% slice_head(n = exemplars_n) %>% pull(id)
  ids_highK_unfund <- df %>% filter(!funded) %>% arrange(desc(K0)) %>% slice_head(n = exemplars_n) %>% pull(id)
  ids_ex <- unique(c(ids_topfund, ids_highK_unfund))
  
  g_grid <- seq(0, max(df$g) + 10, by = 1)
  dg <- 1
  
  mv_df <- expand_grid(id = ids_ex, g = g_grid) %>%
    mutate(
      mv = map2_dbl(id, g, ~posterior_marginal_value(res$posts[[.x]], .y, dg, T, dt, eps, kappa, alpha_out)),
      label = paste0("i=", id)
    )
  
  fig4 <- ggplot(mv_df, aes(x = g, y = mv, color = label)) +
    geom_line(linewidth = 0.8, alpha = 0.95) +
    labs(
      title = "Figure 4. Marginal returns to funding (diminishing returns)",
      x = "Grant size g",
      y = "Posterior marginal value (expected papers per $)",
      color = "Researcher"
    ) +
    theme_sim(base_size) +
    theme(legend.position = "bottom")
  
  cap4 <- paste(
    "Curves show posterior expected marginal gains dE[papers]/dg versus grant size. Optimal funding spreads funds until",
    "marginal returns approximately equalize across funded researchers, rather than concentrating funding past diminishing returns."
  )
  
  # ---------- Figure 5: allocation vs baseline EV ----------
  df5 <- df %>%
    mutate(
      baseline_EV = expected_papers_EV_no_decay(K0, R0, T, dt, eps, kappa, alpha_out),
      funded_EV = expected_papers_EV_no_decay(K0, R0 + g, T, dt, eps, kappa, alpha_out)
    )
  
  fig5 <- ggplot(df5, aes(x = baseline_EV, y = g, color = K0, shape = funded)) +
    geom_point(size = 2.1, alpha = 0.9) +
    scale_color_continuous_safe("True knowledge K0") +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19)) +
    labs(
      title = "Figure 5. Funding is not simply proportional to baseline productivity",
      x = "Baseline expected papers over horizon",
      y = "Allocated grant g",
      shape = "Funded"
    ) +
    theme_sim(base_size)
  
  cap5 <- paste(
    "Baseline output reflects both knowledge and resources. Model-based funding targets high marginal gains—often researchers",
    "who are productive but still resource-bottlenecked—rather than scaling grants mechanically with baseline output."
  )
  
  # ---------- Figure 6: time paths for exemplars ----------
  times <- seq(0, T, by = dt)
  if (tail(times, 1) < T) times <- c(times, T)
  
  cumtrap <- function(t, y) {
    out <- numeric(length(t))
    for (j in 2:length(t)) out[j] <- out[j-1] + 0.5 * (y[j-1] + y[j]) * (t[j] - t[j-1])
    out
  }
  
  time_paths_one <- function(K0, R0, g, eps, kappa, alpha_out, times) {
    cap0 <- kappa * R0
    cap1 <- kappa * (R0 + g)
    
    K0_t <- K_logistic_closed_form(times, rep(K0, length(times)), rep(cap0, length(times)), eps)
    K1_t <- K_logistic_closed_form(times, rep(K0, length(times)), rep(cap1, length(times)), eps)
    
    lam0 <- lambda_rate(K0_t, R0, alpha_out)
    lam1 <- lambda_rate(K1_t, R0 + g, alpha_out)
    
    tibble(
      t = rep(times, 2),
      K = c(K0_t, K1_t),
      cumP = c(cumtrap(times, lam0), cumtrap(times, lam1)),
      scenario = rep(c("No funding", "With funding"), each = length(times))
    )
  }
  
  ids_t <- df %>% arrange(desc(g)) %>% slice_head(n = exemplars_n) %>% pull(id)
  ts_df <- map_dfr(ids_t, function(i) {
    dat <- time_paths_one(df$K0[i], df$R0[i], df$g[i], eps, kappa, alpha_out, times)
    dat$id <- i
    dat
  })
  
  fig6a <- ggplot(ts_df, aes(x = t, y = K, color = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ id, scales = "free_y") +
    labs(
      title = "Figure 6a. Knowledge trajectories for exemplars",
      x = "Time",
      y = "Knowledge K(t)",
      color = ""
    ) +
    theme_sim(base_size) +
    theme(legend.position = "bottom")
  
  fig6b <- ggplot(ts_df, aes(x = t, y = cumP, color = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ id, scales = "free_y") +
    labs(
      title = "Figure 6b. Cumulative expected output for exemplars",
      x = "Time",
      y = "Cumulative expected papers",
      color = ""
    ) +
    theme_sim(base_size) +
    theme(legend.position = "bottom")
  
  cap6 <- paste(
    "Exemplar trajectories compare K(t) and cumulative expected output with vs without funding. Funding raises resources directly",
    "and increases the resource-dependent knowledge ceiling, producing higher trajectories when researchers are resource-bottlenecked."
  )
  
  # ---------- Figure 7 helpers (aggregate across runs; potentially expensive) ----------
  allocate_naive_publications <- function(p, B) {
    w <- p + 1e-6
    B * (w / sum(w))
  }
  
  make_aggregate_df <- function(R_reps = 60, seed = 1, verbose = FALSE) {
    set.seed(seed)
    map_dfr(seq_len(R_reps), function(r) {
      rr <- run_demo_fast(
        n = par$n,
        k_min = par$k_min, k_shape = par$k_shape,
        mu_R = par$mu_R, sd_R = par$sd_R,
        alpha_out = par$alpha_out, eps = par$eps, kappa = par$kappa,
        T = par$T, dt = par$dt,
        tau = par$tau, xi = par$xi,
        B = par$B, delta = par$delta,
        M = par$M,
        verbose = verbose
      )
      d <- tibble(K0 = rr$truth$K0, R0 = rr$truth$R0, p = rr$obs$p, g_opt = rr$allocation$g)
      base <- sum(expected_papers_EV_no_decay(d$K0, d$R0, rr$params$T, rr$params$dt, rr$params$eps, rr$params$kappa, rr$params$alpha_out))
      opt  <- sum(expected_papers_EV_no_decay(d$K0, d$R0 + d$g_opt, rr$params$T, rr$params$dt, rr$params$eps, rr$params$kappa, rr$params$alpha_out))
      g_nv <- allocate_naive_publications(d$p, rr$params$B)
      nv   <- sum(expected_papers_EV_no_decay(d$K0, d$R0 + g_nv, rr$params$T, rr$params$dt, rr$params$eps, rr$params$kappa, rr$params$alpha_out))
      tibble(baseline = base, optimal = opt, naive = nv)
    }) %>%
      pivot_longer(everything(), names_to = "policy", values_to = "total_EV")
  }
  
  plot_aggregate <- function(agg_long) {
    ggplot(agg_long, aes(x = total_EV, fill = policy)) +
      geom_density(alpha = 0.35, linewidth = 0.3) +
      labs(
        title = "Figure 7. Aggregate welfare comparison across runs",
        x = "Total expected papers over horizon",
        y = "Density",
        fill = "Policy"
      ) +
      theme_sim(base_size) +
      theme(legend.position = "bottom")
  }
  
  cap7 <- paste(
    "Distributions compare total expected output under no funding, model-based allocation, and a naive publication-proportional rule.",
    "This summarizes typical welfare gains and variability rather than a single illustrative run."
  )
  
  list(
    figures = list(
      fig1_KR_contours = fig1,
      fig2_p_vs_K_confounds = fig2,
      fig3_posterior_vs_truth = fig3,
      fig4_marginal_value_curves = fig4,
      fig5_allocation_vs_baseline = fig5,
      fig6a_K_timepaths = fig6a,
      fig6b_output_timepaths = fig6b
    ),
    captions = list(
      fig1 = cap1,
      fig2 = cap2,
      fig3 = cap3,
      fig4 = cap4,
      fig5 = cap5,
      fig6 = cap6,
      fig7 = cap7
    ),
    aggregate_helpers = list(
      make_aggregate_df = make_aggregate_df,
      plot_aggregate = plot_aggregate
    )
  )
}

# ============================================================
# 7) Convenience wrapper: run + attach plots + captions
# ============================================================

run_demo_fast_with_plots <- function(...) {
  res <- run_demo_fast(...)
  out <- make_plots(res)
  res$plots <- out$figures
  res$captions <- out$captions
  res$aggregate_helpers <- out$aggregate_helpers
  res
}

# ============================================================
# Example usage (uncomment)
# ============================================================

# res <- run_demo_fast_with_plots(verbose = TRUE)
# res$eval
# res$plots$fig1_KR_contours
# res$captions$fig1
#
# # Aggregate Figure 7 (expensive; start with small reps):
# agg <- res$aggregate_helpers$make_aggregate_df(R_reps = 30, seed = 1, verbose = FALSE)
# fig7 <- res$aggregate_helpers$plot_aggregate(agg)
# fig7
# res$captions$fig7
