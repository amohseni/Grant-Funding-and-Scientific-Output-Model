# app.R
# ============================================================
# Grant-funding model explorer (Shiny)
# ------------------------------------------------------------
# A fast, self-contained Shiny app for:
#   - sampling a researcher population (K0, R0),
#   - generating noisy signals (publications p, optional proposal signal s),
#   - posterior approximation by importance sampling,
#   - greedy grant allocation by posterior marginal value,
#   - standardized plots (Figures 1–7) that update with parameters.
#
# Run:
#   shiny::runApp("app.R")   # or just Run App in RStudio
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
# Core model code (self-contained)
# ============================================================

# ----- Distributions -----

rpareto <- function(n, k_min, shape) {
  u <- runif(n)
  k_min / (u)^(1 / shape)
}

rtruncnorm_pos <- function(n, mu, sd) {
  out <- numeric(0)
  while (length(out) < n) {
    x <- rnorm(n, mu, sd)
    out <- c(out, x[x > 0])
  }
  out[1:n]
}

# ----- Model primitives -----

lambda_rate <- function(K, R, alpha_out) {
  alpha_out * (K * R) / (K + R)
}

loglik_obs <- function(p_obs, s_obs, K0, R0, tau, alpha_out, xi) {
  lam0 <- lambda_rate(K0, R0, alpha_out)
  ll_p <- dpois(p_obs, lambda = tau * lam0, log = TRUE)
  ll_s <- if (is.finite(xi)) dnorm(s_obs, mean = K0, sd = xi, log = TRUE) else 0
  ll_p + ll_s
}

posterior_samples <- function(p_obs, s_obs, M,
                              k_min, k_shape, mu_R, sd_R,
                              tau, alpha_out, xi) {
  K0 <- rpareto(M, k_min, k_shape)
  R0 <- rtruncnorm_pos(M, mu_R, sd_R)
  ll <- loglik_obs(p_obs, s_obs, K0, R0, tau, alpha_out, xi)
  
  ll <- ll - max(ll)
  w <- exp(ll)
  w <- w / sum(w)
  
  list(K0 = K0, R0 = R0, w = w)
}

# ----- Fast expected value engine (no MC) -----

K_logistic_closed_form <- function(t, K0, cap, eps) {
  out <- numeric(length(K0))
  idx <- K0 > 0
  out[idx] <- cap[idx] / (1 + (cap[idx] / K0[idx] - 1) * exp(-eps * t))
  out
}

expected_papers_EV_no_decay <- function(K0, R, T, dt, eps, kappa, alpha_out) {
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
  R <- post$R0 + g
  V <- expected_papers_EV_no_decay(post$K0, R, T, dt, eps, kappa, alpha_out)
  sum(post$w * V)
}

posterior_marginal_value <- function(post, g, dg, T, dt, eps, kappa, alpha_out) {
  (posterior_expected_value(post, g + dg, T, dt, eps, kappa, alpha_out) -
     posterior_expected_value(post, g,      T, dt, eps, kappa, alpha_out)) / dg
}

# ----- Allocation -----

allocate_greedy <- function(posts, B, delta, T, dt, eps, kappa, alpha_out,
                            dg = delta, verbose = FALSE) {
  n <- length(posts)
  stopifnot(B >= 0, delta > 0, dg > 0)
  
  g <- rep(0, n)
  remaining <- B
  
  mv <- vapply(seq_len(n), function(i) {
    posterior_marginal_value(posts[[i]], g[i], dg, T, dt, eps, kappa, alpha_out)
  }, numeric(1))
  
  step <- 0
  while (remaining >= delta) {
    step <- step + 1
    i_star <- which.max(mv)
    g[i_star] <- g[i_star] + delta
    remaining <- remaining - delta
    mv[i_star] <- posterior_marginal_value(posts[[i_star]], g[i_star], dg, T, dt, eps, kappa, alpha_out)
    
    if (verbose && step %% 10 == 0) {
      message(sprintf("step=%d remaining=%.2f funded=%d g[%d]=%.2f mv=%.6f",
                      step, remaining, i_star, i_star, g[i_star], mv[i_star]))
    }
  }
  
  list(g = g, remaining = remaining, mv = mv)
}

# ----- One-run simulation -----

run_demo_fast <- function(
    seed = 1,
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
    xi = 3.0,
    # funding
    B = 60,
    delta = 1,
    # posterior approximation
    M = 1500,
    verbose = FALSE
) {
  set.seed(seed)
  
  K0_true <- rpareto(n, k_min, k_shape)
  R0_true <- rtruncnorm_pos(n, mu_R, sd_R)
  
  lam0_true <- lambda_rate(K0_true, R0_true, alpha_out)
  p_obs <- rpois(n, tau * lam0_true)
  s_obs <- if (is.finite(xi)) rnorm(n, mean = K0_true, sd = xi) else rep(0, n)
  
  posts <- vector("list", n)
  for (i in seq_len(n)) {
    posts[[i]] <- posterior_samples(
      p_obs = p_obs[i], s_obs = s_obs[i], M = M,
      k_min = k_min, k_shape = k_shape, mu_R = mu_R, sd_R = sd_R,
      tau = tau, alpha_out = alpha_out, xi = xi
    )
  }
  
  alloc <- allocate_greedy(posts, B, delta, T, dt, eps, kappa, alpha_out, dg = delta, verbose = verbose)
  
  baseline <- sum(expected_papers_EV_no_decay(K0_true, R0_true,           T, dt, eps, kappa, alpha_out))
  funded   <- sum(expected_papers_EV_no_decay(K0_true, R0_true + alloc$g, T, dt, eps, kappa, alpha_out))
  
  list(
    params = list(
      seed = seed,
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
# Plotting utilities (standardized aesthetics)
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

source("sweep.R", local = FALSE)

make_plots <- function(res, base_size = 12, contour_bins = 12, exemplars_n = 3) {
  df <- build_researcher_df(res)
  par <- res$params
  
  alpha_out <- par$alpha_out
  eps <- par$eps
  kappa <- par$kappa
  T <- par$T
  dt <- par$dt
  
  # --- Figure 1: (K,R) contours ---
  R_grid <- seq(min(df$R0), max(df$R_funded), length.out = 160)
  K_grid <- seq(min(df$K0), max(df$K0), length.out = 160)
  grid <- expand_grid(K = K_grid, R = R_grid) %>% mutate(lambda = lambda_rate(K, R, alpha_out))
  
  fig1 <- ggplot() +
    geom_contour(data = grid, aes(x = R, y = K, z = lambda),
                 bins = contour_bins, linewidth = 0.3, alpha = 0.85) +
    geom_point(data = df, aes(x = R0, y = K0, shape = funded), size = 2.1, alpha = 0.9) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19)) +
    labs(
      title = "Figure 1. Bottleneck geometry in (K,R)-space",
      x = "Baseline resources R0",
      y = "Initial knowledge K0",
      shape = "Funded"
    ) +
    theme_sim(base_size)
  
  # --- Figure 2: p vs K (confounding) ---
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
  
  # --- Figure 3: posterior vs truth (K) ---
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
  
  # --- Figure 4: marginal value curves for exemplars ---
  ids_topfund <- df %>% arrange(desc(g)) %>% slice_head(n = exemplars_n) %>% pull(id)
  ids_highK_unfund <- df %>% filter(!funded) %>% arrange(desc(K0)) %>% slice_head(n = exemplars_n) %>% pull(id)
  ids_ex <- unique(c(ids_topfund, ids_highK_unfund))
  
  g_grid <- seq(0, max(df$g) + 10, by = 1)
  dg <- 1
  
  mv_df <- expand_grid(id = ids_ex, g = g_grid) %>%
    mutate(
      mv = map2_dbl(id, g, ~posterior_marginal_value(res$posts[[.x]], .y, dg,
                                                     T, dt, eps, kappa, alpha_out)),
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
  
  # --- Figure 5: allocation vs baseline EV ---
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
  
  # --- Figure 6: time paths for exemplars (K(t) and cumulative papers) ---
  times <- seq(0, T, by = dt)
  if (tail(times, 1) < T) times <- c(times, T)
  
  cumtrap <- function(t, y) {
    out <- numeric(length(t))
    for (j in 2:length(t)) out[j] <- out[j - 1] + 0.5 * (y[j - 1] + y[j]) * (t[j] - t[j - 1])
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
  
  # --- Figure 7 helpers: aggregate across runs (optional) ---
  allocate_naive_publications <- function(p, B) {
    w <- p + 1e-6
    B * (w / sum(w))
  }
  
  make_aggregate_df <- function(R_reps = 40, seed = 1, verbose = FALSE) {
    set.seed(seed)
    map_dfr(seq_len(R_reps), function(r) {
      rr <- run_demo_fast(
        seed = seed + r,
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
  
  list(
    df = df,
    figures = list(
      fig1 = fig1, fig2 = fig2, fig3 = fig3, fig4 = fig4,
      fig5 = fig5, fig6a = fig6a, fig6b = fig6b
    ),
    aggregate_helpers = list(make_aggregate_df = make_aggregate_df, plot_aggregate = plot_aggregate)
  )
}

# ============================================================
# Shiny app
# ============================================================

ui <- fluidPage(
  # MathJax support
  withMathJax(),
  
  tags$head(
    tags$style(HTML("
      .control-label { font-weight: 600; }
      .shiny-output-error-validation { color: #b00020; font-weight: 600; }
      
      /* Description box */
      .description-box {
        background-color: #f8f9fa;
        border-left: 4px solid #495057;
        padding: 16px 20px;
        margin-bottom: 20px;
        font-size: 14px;
        line-height: 1.6;
      }
      .description-box p {
        margin: 0 0 10px 0;
      }
      .description-box p:last-child {
        margin-bottom: 0;
      }
      
      /* Collapsible panels */
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
      
      /* Parameter help text */
      .param-help {
        font-size: 11px;
        color: #666;
        margin-top: -8px;
        margin-bottom: 12px;
        line-height: 1.4;
      }
      
      /* Primary controls section */
      .primary-controls {
        background-color: #fff;
        padding: 12px 0;
        margin-bottom: 8px;
      }
      .primary-controls .form-group { margin-bottom: 14px; }
      
      /* Plot explanation box */
      .plot-explanation {
        background-color: #f0f4f8;
        border-left: 3px solid #4a90d9;
        padding: 10px 14px;
        margin-bottom: 12px;
        font-size: 13px;
        line-height: 1.5;
      }
      
      /* Run button */
      .btn-run {
        width: 100%;
        margin-top: 8px;
        margin-bottom: 8px;
      }
      
      /* Summary box styling */
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
  
  # Description block
  fluidRow(
    column(12,
      div(class = "description-box",
        p("This model investigates how funders should allocate grants when researchers 
           differ in knowledge (\\(K\\)) and resources (\\(R\\)), but only noisy signals 
           are observable. We assume researchers produce knowledge as a function of their 
           knowledge and resources—output is bottlenecked by whichever input is scarcer.
           Funding relaxes the resource constraint."),
        p("The funder's challenge: a researcher with high \\(K\\) and low \\(R\\) 
           (resource-constrained, high marginal return) can look identical in publication 
           count to one with low \\(K\\) and high \\(R\\) (knowledge-constrained, low 
           marginal return). We compare various allocation rules under differing information 
           regimes: (1) observing only past publications, (2) observing an additional noisy 
           signal of knowledge (e.g., grant proposal quality).")
      )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # ---- Primary controls (always visible) ----
      div(class = "primary-controls",
        sliderInput("n", "Researchers", min = 10, max = 200, value = 40, step = 5),
        div(class = "param-help", "Number of researchers in the population."),
        
        sliderInput("B", "Total budget", min = 0, max = 300, value = 60, step = 5),
        div(class = "param-help", "Total funding to allocate across all researchers."),
        
        checkboxInput("use_proposal", "Use proposal signal", value = TRUE),
        div(class = "param-help", 
            "When checked, funder observes a noisy signal of each researcher's true knowledge 
             (e.g., grant proposal quality). When unchecked, funder sees only publication counts."),
        
        numericInput("seed", "Random seed", value = 1, min = 1, step = 1),
        div(class = "param-help", "Controls random draws; change to see different populations.")
      ),
      
      actionButton("run", "Run simulation", class = "btn-primary btn-run"),
      
      hr(),
      
      # ---- Collapsible: Prior distributions ----
      tags$div(class = "panel",
        tags$div(class = "panel-heading", 
          `data-toggle` = "collapse", `data-target` = "#panel-priors",
          tags$h4(class = "panel-title", 
            "Prior distributions",
            tags$span(class = "caret-icon", HTML("&#9662;"))
          )
        ),
        tags$div(id = "panel-priors", class = "collapse panel-body",
          sliderInput("k_min", "Knowledge minimum", min = 0.1, max = 5, value = 1, step = 0.1),
          div(class = "param-help", "Lower bound of Pareto distribution for initial knowledge \\(K_0\\)."),
          
          sliderInput("k_shape", "Knowledge shape", min = 1.1, max = 5, value = 2.2, step = 0.1),
          div(class = "param-help", "Pareto tail parameter; lower = heavier tail (more high-\\(K\\) researchers)."),
          
          sliderInput("mu_R", "Resources mean", min = 1, max = 30, value = 10, step = 0.5),
          div(class = "param-help", "Mean of truncated normal distribution for baseline resources \\(R_0\\)."),
          
          sliderInput("sd_R", "Resources SD", min = 0.1, max = 15, value = 3, step = 0.5),
          div(class = "param-help", "Standard deviation of baseline resources.")
        )
      ),
      
      # ---- Collapsible: Knowledge dynamics ----
      tags$div(class = "panel",
        tags$div(class = "panel-heading",
          `data-toggle` = "collapse", `data-target` = "#panel-dynamics",
          tags$h4(class = "panel-title",
            "Knowledge dynamics",
            tags$span(class = "caret-icon", HTML("&#9662;"))
          )
        ),
        tags$div(id = "panel-dynamics", class = "collapse panel-body",
          sliderInput("alpha_out", "Output scale (α)", min = 0.01, max = 1, value = 0.12, step = 0.01),
          div(class = "param-help", "Scales paper production rate; higher = more papers per unit time."),
          
          sliderInput("eps", "Growth rate (ε)", min = 0.0, max = 1.0, value = 0.25, step = 0.01),
          div(class = "param-help", "Rate at which knowledge grows over time (logistic)."),
          
          sliderInput("kappa", "Resource cap (κ)", min = 0.5, max = 10, value = 2.5, step = 0.1),
          div(class = "param-help", "Knowledge ceiling is \\(\\kappa \\times R\\); more resources → higher potential knowledge."),
          
          sliderInput("T", "Time horizon", min = 1, max = 60, value = 24, step = 1),
          div(class = "param-help", "How far into the future to project output (e.g., months or years).")
        )
      ),
      
      # ---- Collapsible: Observation model ----
      tags$div(class = "panel",
        tags$div(class = "panel-heading",
          `data-toggle` = "collapse", `data-target` = "#panel-observation",
          tags$h4(class = "panel-title",
            "Observation model",
            tags$span(class = "caret-icon", HTML("&#9662;"))
          )
        ),
        tags$div(id = "panel-observation", class = "collapse panel-body",
          sliderInput("tau", "Observation window (τ)", min = 1, max = 24, value = 6, step = 1),
          div(class = "param-help", "Duration over which past publications are observed (pre-grant)."),
          
          sliderInput("xi", "Proposal noise (ξ)", min = 0.1, max = 20, value = 3, step = 0.1),
          div(class = "param-help", "SD of proposal signal around true \\(K\\); smaller = more informative proposals.")
        )
      ),
      
      # ---- Collapsible: Simulation settings ----
      tags$div(class = "panel",
        tags$div(class = "panel-heading",
          `data-toggle` = "collapse", `data-target` = "#panel-simulation",
          tags$h4(class = "panel-title",
            "Simulation settings",
            tags$span(class = "caret-icon", HTML("&#9662;"))
          )
        ),
        tags$div(id = "panel-simulation", class = "collapse panel-body",
          sliderInput("dt", "Integration step", min = 0.05, max = 2, value = 0.25, step = 0.05),
          div(class = "param-help", "Numerical precision; smaller = more accurate but slower."),
          
          sliderInput("delta", "Allocation step", min = 0.5, max = 10, value = 1, step = 0.5),
          div(class = "param-help", "Granularity of greedy allocation; smaller = finer but slower."),
          
          sliderInput("M", "Posterior samples", min = 200, max = 6000, value = 1500, step = 100),
          div(class = "param-help", "Monte Carlo samples for Bayesian inference; more = smoother posteriors."),
          
          sliderInput("exemplars_n", "Exemplar count", min = 1, max = 8, value = 3, step = 1),
          div(class = "param-help", "Number of individual researchers shown in trajectory plots.")
        )
      ),
      
      hr(),
      helpText("Tip: if slow, increase dt, decrease n, or decrease M.")
    ),
    
    mainPanel(
      width = 9,
      
      div(class = "summary-box",
        verbatimTextOutput("summary", placeholder = TRUE)
      ),
      
      tabsetPanel(id = "main_tabs",
        # ---- NEW DEFAULT TAB: Signal value comparison ----
        tabPanel("Value of Signals",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Compares total expected research output under three 
                  allocation policies: (1) <em>Baseline</em>—no additional funding; 
                  (2) <em>Publications only</em>—funder allocates based solely on observed papers; 
                  (3) <em>With proposal signal</em>—funder also observes a noisy signal of true knowledge.
                  The gap between policies quantifies the value of better information.")
          ),
          plotOutput("fig_signal_value", height = 480)
        ),
        
        tabPanel("(K,R) Geometry",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Each point is a researcher in \\((R_0, K_0)\\) space.
                  Contour lines show iso-output curves (\\(\\lambda = \\alpha K R / (K+R)\\)).
                  Funded researchers (filled) tend to lie above the diagonal—high \\(K\\) relative to \\(R\\)—
                  because they have the highest marginal return to additional resources.")
          ),
          plotOutput("fig1", height = 480)
        ),
        
        tabPanel("Confounding",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Publication count (\\(p\\)) vs. true knowledge (\\(K_0\\)), 
                  colored by resources (\\(R_0\\)). The key insight: researchers with different \\((K,R)\\) 
                  combinations can have similar publication records. This confounding makes it hard to 
                  identify who is resource-constrained vs. knowledge-constrained from publications alone.")
          ),
          plotOutput("fig2", height = 480)
        ),
        
        tabPanel("Posterior Check",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> How well the funder's posterior beliefs about \\(K_0\\) 
                  match truth. Points near the diagonal indicate accurate inference. Error bars show 
                  90% credible intervals. With the proposal signal, posteriors should be tighter and 
                  better calibrated.")
          ),
          plotOutput("fig3", height = 480)
        ),
        
        tabPanel("Marginal Returns",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Expected marginal value of funding (additional papers 
                  per dollar) as a function of grant size, for selected researchers. Curves slope downward 
                  due to diminishing returns: as \\(R\\) increases, the bottleneck shifts to \\(K\\).")
          ),
          plotOutput("fig4", height = 480)
        ),
        
        tabPanel("Allocation Pattern",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Grant allocation (\\(g\\)) vs. baseline expected output.
                  Optimal funding is <em>not</em> proportional to baseline productivity—high-\\(K\\), 
                  low-\\(R\\) researchers may have low baseline output but high marginal value.")
          ),
          plotOutput("fig5", height = 480)
        ),
        
        tabPanel("Trajectories",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Time paths of knowledge \\(K(t)\\) and cumulative 
                  expected papers for top-funded researchers, with vs. without their grant.
                  Funding raises the knowledge ceiling (\\(\\kappa R\\)), enabling sustained growth.")
          ),
          fluidRow(
            column(6, plotOutput("fig6a", height = 400)),
            column(6, plotOutput("fig6b", height = 400))
          )
        ),
        
        tabPanel("Aggregates",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Distribution of total expected output across many
                  simulation runs, comparing policies. This reveals robustness: does optimal allocation
                  consistently outperform naive allocation across different random populations?")
          ),
          fluidRow(
            column(4,
              sliderInput("R_reps", "Replicate runs", min = 10, max = 200, value = 40, step = 10),
              actionButton("run_agg", "Compute", class = "btn-primary")
            ),
            column(8,
              plotOutput("fig7", height = 420)
            )
          )
        ),

        tabPanel("Parameter Sweeps",
          div(class = "plot-explanation",
            HTML("<strong>What this shows:</strong> Systematic parameter sweeps across model configurations,
                  averaged over multiple random seeds. Each cell in the heatmap represents the mean outcome
                  for that parameter combination. This reveals which structural features of the research
                  population drive the value of different allocation strategies.")
          ),
          fluidRow(
            column(4,
              selectInput("sweep_choice", "Sweep type",
                choices = setNames(
                  names(SWEEP_CONFIGS),
                  sapply(SWEEP_CONFIGS, `[[`, "name")
                )
              ),
              sliderInput("sweep_seeds", "Number of seeds",
                          min = 3, max = 15, value = 5, step = 1),
              div(class = "param-help",
                  "More seeds = more precise estimates but slower computation."),
              actionButton("run_sweep_btn", "Run Sweep",
                           class = "btn-primary btn-run"),
              hr(),
              uiOutput("sweep_description")
            ),
            column(8,
              plotOutput("sweep_primary_plot", height = 480),
              hr(),
              plotOutput("sweep_secondary_plot", height = 380)
            )
          )
        )
      )
    )
  ),
  
  # Bootstrap collapse functionality
  tags$script(HTML("
    $(document).on('click', '.panel-heading', function() {
      $(this).next('.collapse').slideToggle(200);
      $(this).find('.caret-icon').toggleClass('collapsed');
    });
  "))
)

server <- function(input, output, session) {
  
  # ---- Run both with and without signal for comparison ----
  # ignoreNULL = FALSE ensures it runs on startup with initial values
  res_comparison <- reactiveVal(NULL)
  
  observe({
    # This runs on startup and whenever "Run" is clicked
    input$run  # Take dependency on button
    
    isolate({
      validate(
        need(input$dt > 0, "dt must be > 0."),
        need(input$T > 0, "T must be > 0."),
        need(input$k_shape > 1, "Pareto shape must be > 1."),
        need(input$sd_R > 0, "sd_R must be > 0."),
        need(input$M >= 100, "M should be at least 100.")
      )
      
      # Common parameters
      common_params <- list(
        seed = input$seed,
        n = input$n,
        k_min = input$k_min,
        k_shape = input$k_shape,
        mu_R = input$mu_R,
        sd_R = input$sd_R,
        alpha_out = input$alpha_out,
        eps = input$eps,
        kappa = input$kappa,
        T = input$T,
        dt = input$dt,
        tau = input$tau,
        B = input$B,
        delta = input$delta,
        M = input$M,
        verbose = FALSE
      )
      
      # Run WITH proposal signal
      res_with <- do.call(run_demo_fast, c(common_params, xi = input$xi))
      
      # Run WITHOUT proposal signal (xi = Inf)
      res_without <- do.call(run_demo_fast, c(common_params, xi = Inf))
      
      # Naive allocation (proportional to publications)
      allocate_naive <- function(p, B) {
        w <- p + 1e-6
        B * (w / sum(w))
      }
      g_naive <- allocate_naive(res_with$obs$p, input$B)
      naive_total <- sum(expected_papers_EV_no_decay(
        res_with$truth$K0, 
        res_with$truth$R0 + g_naive, 
        input$T, input$dt, input$eps, input$kappa, input$alpha_out
      ))
      
      res_comparison(list(
        with_signal = res_with,
        without_signal = res_without,
        naive_total = naive_total,
        use_proposal = input$use_proposal
      ))
    })
  })
  
  # ---- Primary result (respects checkbox) ----
  res <- reactive({
    comp <- res_comparison()
    req(comp)  # Wait until comparison is computed
    if (isTRUE(input$use_proposal)) {
      comp$with_signal
    } else {
      comp$without_signal
    }
  })
  
  plots <- reactive({
    make_plots(
      res(),
      base_size = 12,
      contour_bins = 12,
      exemplars_n = input$exemplars_n
    )
  })
  
  # ---- Summary output ----
  output$summary <- renderPrint({
    comp <- res_comparison()
    req(comp)
    r <- res()
    
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("                        KEY RESULTS\n")
    cat("═══════════════════════════════════════════════════════════════\n\n")
    
    cat(sprintf("Population: n = %d researchers | Budget: B = %.0f | Horizon: T = %.0f\n\n",
                r$params$n, r$params$B, r$params$T))
    
    cat("EXPECTED TOTAL OUTPUT (papers over horizon)\n")
    cat("───────────────────────────────────────────\n")
    cat(sprintf("  Baseline (no funding):        %7.1f\n", r$eval$expected_total_baseline))
    cat(sprintf("  Naive (∝ publications):       %7.1f  (+%.1f)\n", 
                comp$naive_total, comp$naive_total - r$eval$expected_total_baseline))
    cat(sprintf("  Optimal (pubs only):          %7.1f  (+%.1f)\n",
                comp$without_signal$eval$expected_total_funded,
                comp$without_signal$eval$expected_gain))
    cat(sprintf("  Optimal (pubs + proposal):    %7.1f  (+%.1f)\n",
                comp$with_signal$eval$expected_total_funded,
                comp$with_signal$eval$expected_gain))
    cat("\n")
    
    signal_value <- comp$with_signal$eval$expected_total_funded - 
                    comp$without_signal$eval$expected_total_funded
    cat(sprintf("VALUE OF PROPOSAL SIGNAL:       %+.1f papers (%.1f%% of gain)\n",
                signal_value,
                100 * signal_value / comp$with_signal$eval$expected_gain))
    cat("\n")
    
    cat(sprintf("Researchers funded: %d / %d\n",
                sum(r$allocation$g > 0), r$params$n))
    cat(sprintf("Signal mode: %s\n", 
                if (isTRUE(comp$use_proposal)) 
                  sprintf("Proposal + publications (ξ = %.1f)", r$params$xi) 
                else 
                  "Publications only"))
  })
  
  # ---- NEW: Signal value comparison plot ----
  output$fig_signal_value <- renderPlot({
    comp <- res_comparison()
    req(comp)
    
    df_compare <- tibble(
      Policy = factor(
        c("Baseline\n(no funding)", 
          "Naive\n(∝ publications)", 
          "Optimal\n(pubs only)", 
          "Optimal\n(pubs + proposal)"),
        levels = c("Baseline\n(no funding)", 
                   "Naive\n(∝ publications)", 
                   "Optimal\n(pubs only)", 
                   "Optimal\n(pubs + proposal)")
      ),
      Total = c(
        comp$with_signal$eval$expected_total_baseline,
        comp$naive_total,
        comp$without_signal$eval$expected_total_funded,
        comp$with_signal$eval$expected_total_funded
      ),
      Type = c("No funding", "Allocation", "Allocation", "Allocation")
    )
    
    baseline <- comp$with_signal$eval$expected_total_baseline
    
    ggplot(df_compare, aes(x = Policy, y = Total, fill = Policy)) +
      geom_col(width = 0.7, alpha = 0.9) +
      geom_hline(yintercept = baseline, linetype = "dashed", alpha = 0.5) +
      geom_text(aes(label = sprintf("%.1f", Total)), 
                vjust = -0.5, size = 4, fontface = "bold") +
      geom_text(aes(label = sprintf("+%.1f", Total - baseline), y = Total + 2),
                vjust = -1.5, size = 3.2, color = "grey40",
                data = . %>% filter(Policy != "Baseline\n(no funding)")) +
      scale_fill_manual(values = c(
        "Baseline\n(no funding)" = "#bdbdbd",
        "Naive\n(∝ publications)" = "#90caf9",
        "Optimal\n(pubs only)" = "#42a5f5",
        "Optimal\n(pubs + proposal)" = "#1565c0"
      )) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        title = "Total expected research output by allocation policy",
        subtitle = sprintf("n = %d researchers, B = %.0f budget, T = %.0f time horizon",
                           comp$with_signal$params$n,
                           comp$with_signal$params$B,
                           comp$with_signal$params$T),
        x = NULL,
        y = "Expected total papers"
      ) +
      theme_sim(base_size = 13) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 11),
        plot.subtitle = element_text(color = "grey40", size = 11)
      )
  }, res = 110)
  
  # ---- Existing plots ----
  output$fig1 <- renderPlot({ plots()$figures$fig1 }, res = 110)
  output$fig2 <- renderPlot({ plots()$figures$fig2 }, res = 110)
  output$fig3 <- renderPlot({ plots()$figures$fig3 }, res = 110)
  output$fig4 <- renderPlot({ plots()$figures$fig4 }, res = 110)
  output$fig5 <- renderPlot({ plots()$figures$fig5 }, res = 110)
  output$fig6a <- renderPlot({ plots()$figures$fig6a }, res = 110)
  output$fig6b <- renderPlot({ plots()$figures$fig6b }, res = 110)
  
  # ---- Aggregate figure ----
  agg <- eventReactive(input$run_agg, {
    ph <- plots()$aggregate_helpers
    ph$make_aggregate_df(R_reps = input$R_reps, seed = input$seed, verbose = FALSE)
  }, ignoreInit = TRUE)
  
  output$fig7 <- renderPlot({
    req(agg())
    ph <- plots()$aggregate_helpers
    ph$plot_aggregate(agg())
  }, res = 110)

  # ---- Parameter Sweeps ----

  output$sweep_description <- renderUI({
    config <- SWEEP_CONFIGS[[input$sweep_choice]]
    req(config)
    div(class = "plot-explanation",
        HTML(paste0("<p>", config$description, "</p>")))
  })

  sweep_results <- eventReactive(input$run_sweep_btn, {
    sweep_name <- input$sweep_choice
    n_seeds    <- input$sweep_seeds

    config      <- SWEEP_CONFIGS[[sweep_name]]
    params_grid <- config$grid_fn()
    total       <- nrow(params_grid) * n_seeds

    withProgress(
      message = paste("Running", config$name, "sweep..."),
      value = 0,
      {
        raw_df <- run_sweep(
          sweep_name  = sweep_name,
          seeds       = seq_len(n_seeds),
          progress_fn = function(i, total) {
            incProgress(1 / total,
                        detail = sprintf("Run %d / %d", i, total))
          }
        )

        summary_df <- summarize_sweep(raw_df, config)
        list(raw = raw_df, summary = summary_df, config = config)
      }
    )
  }, ignoreInit = TRUE)

  output$sweep_primary_plot <- renderPlot({
    res <- sweep_results()
    req(res)
    plot_sweep_from_config(res$summary, res$config, which = "primary")
  }, res = 110)

  output$sweep_secondary_plot <- renderPlot({
    res <- sweep_results()
    req(res)
    plot_sweep_from_config(res$summary, res$config, which = "secondary")
  }, res = 110)
}

shinyApp(ui, server)
