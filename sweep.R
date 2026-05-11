# sweep.R
# ============================================================
# Parameter sweep infrastructure — v5 (rewrite of v3-era sweep)
# ============================================================
#
# What changed vs. the v3 sweep:
#   - run_demo_fast() → run_simulation()
#   - Continuous-time (T, dt) → discrete two-round
#   - kappa removed (harmonic-mean dynamics)
#   - R ~ Normal(mu_R, sd_R) → R ~ Pareto(r_min, r_shape)
#   - Single xi → (tau_k, tau_r); eps → epsilon; alpha_out → gamma
#   - Absolute B → dimensionless b = B / (n · E[R])
#   - "Optimal with/without signal vs naive" frame → 9-strategy comparison
#
# Usage:
#   Rscript -e 'source("app.R"); source("sweep.R"); main_sweep()'
#   Or interactively in R after sourcing app.R.
#
# Each sweep saves a raw .rds and an aggregated .rds under ./sweep_results/.
# Plots are returned by plot_sweep_from_config() — save them with ggsave().
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# Null-coalescing operator (in case not already defined upstream)
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# ============================================================
# DEFAULTS
# ============================================================
#
# Match v5 production defaults (run_simulation arg defaults), with
# smaller n and M for sweep speed. At n=50, M=500, one run takes
# ~0.3-0.6s; a 5×3 grid × 20 trials runs in roughly 5 minutes.

SWEEP_BASE_PARAMS <- list(
  n                   = 50,
  k_min               = 1.0,
  k_shape             = 2.0,
  r_min               = 1.0,
  r_shape             = 2.0,
  rho_kr              = 0,
  gamma               = 1.0,
  epsilon             = 0.3,        # v5 production default
  b                   = 0.5,
  n_steps             = 50,
  tau_r               = 1.0,
  tau_k               = 1.0,
  use_resource_signal = TRUE,
  n_pre_rounds        = 0,
  x_seed              = 0.5,
  M                   = 500,        # vs. production 1500
  strategies          = 1:9
)

# ============================================================
# SWEEP CONFIGURATIONS
# ============================================================
#
# Naming convention for metrics extracted from each run:
#   out_S{n}    = total_expected of strategy n (1..9)
#   alpha_S{n}  = alpha (round-1 share) of strategy n
#
# Derived metrics computed in build_metrics_row():
#   fwd_vs_myo_P      = S7 − S4   (Pubs only:    forward gain)
#   fwd_vs_myo_PG     = S8 − S5   (Pubs+Grant:   forward gain)  ← main hypothesis
#   fwd_vs_myo_PS     = S9 − S6   (Pubs+Seed:    forward gain)
#   signal_myo        = S5 − S4   (Myopic:       grant signal value)
#   signal_fwd        = S8 − S7   (Forward:      grant signal value)
#   seed_myo          = S6 − S4   (Myopic:       seed value)
#   seed_fwd          = S9 − S7   (Forward:      seed value)
#   optimal_vs_naive  = S8 − S3   (Best Bayesian over naive baseline)

SWEEP_CONFIGS <- list(
  
  # ----------------------------------------------------------
  # The headline test: when does forward planning beat myopic?
  # Hypothesis: gap opens at high signal noise (tau_k); ε matters less.
  signal_noise = list(
    name = "Forward vs. myopic across signal noise",
    description = paste(
      "Sweep grant-signal noise tau_k (precision of the proposal/peer-review signal)",
      "against knowledge-growth rate epsilon. Forward planning is hypothesized to",
      "outperform myopic only when grant signals are noisy — round-1 grants then",
      "double as information-acquisition. At sharp tau_k forward and myopic should",
      "converge to within MC noise."
    ),
    grid_fn = function() {
      expand.grid(
        tau_k   = c(0.3, 1.0, 3.0, 10.0),
        epsilon = c(0.1, 0.3, 0.6),
        stringsAsFactors = FALSE
      )
    },
    varied_params = c("tau_k", "epsilon"),
    primary_plot = list(
      type       = "heatmap",
      x_var      = "tau_k",
      y_var      = "epsilon",
      fill_var   = "fwd_vs_myo_PG_mean",
      facet_var  = NULL,
      title      = "Forward minus Myopic (Pubs+Grant)",
      fill_label = "F − M\n(expected output)",
      text_var   = "fwd_vs_myo_PG_mean",
      digits     = 2
    ),
    secondary_plot = list(
      type      = "line",
      x_var     = "tau_k",
      y_var     = "fwd_vs_myo_PG_mean",
      color_var = "epsilon",
      facet_var = NULL,
      title     = "Forward gain over myopic vs. grant-signal noise",
      y_label   = "Forward − Myopic (expected output)"
    )
  ),
  
  # ----------------------------------------------------------
  # Value of the grant signal: in which regimes does observing
  # the grant signal beat using publications alone?
  signal_value = list(
    name = "Value of grant signal",
    description = paste(
      "Sweep grant-signal noise tau_k against the heaviness of the knowledge tail.",
      "More heterogeneous populations (low k_shape) should make the grant signal",
      "more valuable: the gap between funding 'the right' researcher and a typical",
      "one is larger. At low signal noise and high inequality, signal-using strategies",
      "should dominate publication-only strategies; at high noise the signal degrades",
      "and the gap shrinks."
    ),
    grid_fn = function() {
      expand.grid(
        tau_k   = c(0.3, 1.0, 3.0, 10.0),
        k_shape = c(1.3, 1.8, 2.5, 3.5),
        stringsAsFactors = FALSE
      )
    },
    varied_params = c("tau_k", "k_shape"),
    primary_plot = list(
      type       = "heatmap",
      x_var      = "tau_k",
      y_var      = "k_shape",
      fill_var   = "signal_fwd_mean",
      facet_var  = NULL,
      title      = "Forward signal value (S8 − S7) across noise and inequality",
      fill_label = "Signal value\n(expected output)",
      text_var   = "signal_fwd_mean",
      digits     = 2
    ),
    secondary_plot = list(
      type      = "line",
      x_var     = "tau_k",
      y_var     = "signal_fwd_mean",
      color_var = "k_shape",
      facet_var = NULL,
      title     = "Forward signal value across noise levels",
      y_label   = "S8 − S7 (expected output)"
    )
  ),
  
  # ----------------------------------------------------------
  # Where does the planner choose to allocate?
  # Diagnostic: alpha (round-1 share) of forward planner across regimes.
  alpha_regime = list(
    name = "Optimal alpha across regimes",
    description = paste(
      "Where does the forward (CE) planner choose to spend? Sweep tau_k against",
      "epsilon. Force B (knowledge dynamics) pushes alpha above 0.5 (front-load to",
      "compound K). Forces D and E (information value of round-1 grants) also push",
      "alpha up at high tau_k. If both forces are small (low epsilon, sharp signal),",
      "the planner should converge to alpha = 0.5 by symmetry."
    ),
    grid_fn = function() {
      expand.grid(
        tau_k   = c(0.3, 1.0, 3.0, 10.0),
        epsilon = c(0.05, 0.2, 0.5, 0.8),
        stringsAsFactors = FALSE
      )
    },
    varied_params = c("tau_k", "epsilon"),
    primary_plot = list(
      type       = "heatmap",
      x_var      = "tau_k",
      y_var      = "epsilon",
      fill_var   = "alpha_S8_mean",
      facet_var  = NULL,
      title      = "Forward (P+G) round-1 share alpha by regime",
      fill_label = "alpha = sum(g1)/(2B)",
      text_var   = "alpha_S8_mean",
      digits     = 3
    ),
    secondary_plot = list(
      type      = "line",
      x_var     = "tau_k",
      y_var     = "alpha_S8_mean",
      color_var = "epsilon",
      facet_var = NULL,
      title     = "Forward alpha across grant-signal noise",
      y_label   = "alpha (round-1 share)"
    )
  ),
  
  # ----------------------------------------------------------
  # Value of seeding: is uniform seeding ever competitive?
  seed_value = list(
    name = "Value of seed grants",
    description = paste(
      "When does adding a uniform seed grant help? Sweep budget fraction b against",
      "knowledge-tail heaviness. Seeding spreads round-1 budget away from targeting,",
      "which should hurt at large b (signals are accurate enough that targeting wins)",
      "and could help at small b (where targeted grants are too small to overcome",
      "bottlenecks). Inequality interacts: heavy tails make targeting more valuable."
    ),
    grid_fn = function() {
      expand.grid(
        b       = c(0.1, 0.3, 0.5, 1.0, 2.0),
        k_shape = c(1.5, 2.0, 3.0),
        stringsAsFactors = FALSE
      )
    },
    varied_params = c("b", "k_shape"),
    primary_plot = list(
      type       = "heatmap",
      x_var      = "b",
      y_var      = "k_shape",
      fill_var   = "seed_fwd_mean",
      facet_var  = NULL,
      title      = "Seed value under forward planner (S9 − S7)",
      fill_label = "Seed value\n(expected output)",
      text_var   = "seed_fwd_mean",
      digits     = 2
    ),
    secondary_plot = list(
      type      = "line",
      x_var     = "b",
      y_var     = "seed_fwd_mean",
      color_var = "k_shape",
      facet_var = NULL,
      title     = "Seed value across budget levels",
      y_label   = "S9 − S7 (expected output)"
    )
  )
)

# ============================================================
# RUNNER
# ============================================================

# Pull total_expected and alpha for every strategy out of a run result
extract_strategy_metrics <- function(res) {
  out <- list()
  for (s in seq_along(res$strategies)) {
    r <- res$strategies[[s]]
    if (is.null(r)) {
      out[[paste0("out_S",   s)]] <- NA_real_
      out[[paste0("alpha_S", s)]] <- NA_real_
    } else {
      out[[paste0("out_S",   s)]] <- r$total_expected
      out[[paste0("alpha_S", s)]] <- if (is.null(r$alpha)) NA_real_ else r$alpha
    }
  }
  out
}

# Derive named comparison metrics from a single-run row of out_S{n} / alpha_S{n}
build_derived_metrics <- function(row) {
  list(
    fwd_vs_myo_P     = row$out_S7 - row$out_S4,
    fwd_vs_myo_PG    = row$out_S8 - row$out_S5,
    fwd_vs_myo_PS    = row$out_S9 - row$out_S6,
    signal_myo       = row$out_S5 - row$out_S4,
    signal_fwd       = row$out_S8 - row$out_S7,
    seed_myo         = row$out_S6 - row$out_S4,
    seed_fwd         = row$out_S9 - row$out_S7,
    optimal_vs_naive = row$out_S8 - row$out_S3
  )
}

run_sweep <- function(sweep_name,
                      seeds       = 1:20,
                      base_params = NULL,
                      progress    = TRUE) {
  
  config <- SWEEP_CONFIGS[[sweep_name]]
  if (is.null(config)) stop("Unknown sweep: ", sweep_name)
  
  grid <- config$grid_fn()
  bp   <- modifyList(SWEEP_BASE_PARAMS, base_params %||% list())
  
  total_iters <- nrow(grid) * length(seeds)
  results <- vector("list", total_iters)
  counter <- 0L
  t_start <- Sys.time()
  
  for (row_i in seq_len(nrow(grid))) {
    row_vals   <- as.list(grid[row_i, , drop = FALSE])
    run_params <- modifyList(bp, row_vals)
    
    for (s in seeds) {
      counter <- counter + 1L
      run_params$seed <- s
      
      res <- do.call(run_simulation, run_params)
      
      strat_metrics <- extract_strategy_metrics(res)
      derived <- build_derived_metrics(strat_metrics)
      
      result_row <- c(row_vals, list(seed = s), strat_metrics, derived)
      results[[counter]] <- as.data.frame(result_row, stringsAsFactors = FALSE)
      
      if (progress && counter %% 5L == 0L) {
        elapsed <- as.numeric(Sys.time() - t_start, units = "secs")
        eta <- elapsed * (total_iters - counter) / counter
        cat(sprintf("  [%d/%d] elapsed %.0fs, ETA %.0fs\n",
                    counter, total_iters, elapsed, eta))
      }
    }
  }
  
  do.call(rbind, results)
}

# ============================================================
# SUMMARIZE
# ============================================================

# Aggregate raw run-level rows to one row per (varied params) cell:
# mean and SE for every numeric column.
summarize_sweep <- function(raw_df, config) {
  varied <- config$varied_params
  numeric_cols <- setdiff(names(raw_df)[sapply(raw_df, is.numeric)],
                          c(varied, "seed"))
  
  raw_df %>%
    group_by(across(all_of(varied))) %>%
    summarise(
      n_trials = n(),
      across(
        all_of(numeric_cols),
        list(mean = ~mean(.x, na.rm = TRUE),
             se   = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

# ============================================================
# PLOTTING
# ============================================================

theme_sweep <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = base_size + 1, margin = margin(b = 6)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25),
      strip.text = element_text(face = "bold"),
      plot.margin = margin(12, 12, 12, 12)
    )
}

scale_fill_continuous_safe <- function(name = waiver(), ...) {
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    scale_fill_viridis_c(name = name, option = "C", ...)
  } else {
    scale_fill_gradient(name = name, low = "grey90", high = "grey10", ...)
  }
}

plot_sweep_heatmap <- function(summary_df, x_var, y_var, fill_var,
                               facet_var = NULL, title = "",
                               fill_label = "", text_var = NULL,
                               digits = 2) {
  df <- summary_df %>%
    mutate(across(all_of(c(x_var, y_var)), as.factor))
  
  p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]],
                      fill = .data[[fill_var]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_continuous_safe(name = fill_label) +
    labs(title = title, x = x_var, y = y_var) +
    theme_sweep(base_size = 12) +
    theme(legend.position = "right")
  
  if (!is.null(text_var) && text_var %in% names(df)) {
    p <- p + geom_text(
      aes(label = sprintf(paste0("%.", digits, "f"), .data[[text_var]])),
      size = 3.2, color = "white", fontface = "bold"
    )
  }
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)), nrow = 1)
  }
  
  p
}

plot_sweep_line <- function(summary_df, x_var, y_var, color_var,
                            facet_var = NULL, title = "",
                            y_label = "") {
  se_var <- sub("_mean$", "_se", y_var)
  has_se <- se_var %in% names(summary_df)
  
  df <- summary_df %>%
    mutate(across(all_of(color_var), as.factor))
  
  p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]],
                      color = .data[[color_var]],
                      fill  = .data[[color_var]],
                      group = .data[[color_var]]))
  
  if (has_se) {
    p <- p + geom_ribbon(
      aes(ymin = .data[[y_var]] - .data[[se_var]],
          ymax = .data[[y_var]] + .data[[se_var]]),
      alpha = 0.15, color = NA
    )
  }
  
  p <- p +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    labs(title = title, x = x_var, y = y_label,
         color = color_var, fill = color_var) +
    theme_sweep(base_size = 12) +
    theme(legend.position = "bottom")
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  }
  
  p
}

plot_sweep_from_config <- function(summary_df, config, which = "primary") {
  spec <- if (which == "primary") config$primary_plot else config$secondary_plot
  if (is.null(spec)) return(NULL)
  
  if (spec$type == "heatmap") {
    plot_sweep_he