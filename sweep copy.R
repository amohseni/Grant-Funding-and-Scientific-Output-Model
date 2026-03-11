# sweep.R
# ============================================================
# Parameter sweep infrastructure for the grant-funding model
# ============================================================
#
# This file defines:
#   1. Pre-configured parameter sweeps for key research questions
#   2. A generic sweep runner (run_sweep)
#   3. Summary and plotting functions for sweep results
#
# Dependencies:
#   - Model functions from app.R (run_demo_fast, expected_papers_EV_no_decay,
#     lambda_rate, etc.) must be loaded before sourcing this file.
#   - ggplot2, dplyr, tidyr, purrr (loaded by app.R)
#
# Usage:
#   From app.R:  source("sweep.R")  -- called automatically
#   Standalone:  source the model functions first, then source("sweep.R")
# ============================================================

# ----- Utilities -----

allocate_naive_publications <- function(p, B) {
  w <- p + 1e-6
  B * (w / sum(w))
}

scale_fill_continuous_safe <- function(name = waiver(), ...) {
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    scale_fill_viridis_c(name = name, option = "C", ...)
  } else {
    scale_fill_gradient(name = name, low = "grey90", high = "grey10", ...)
  }
}

# ----- Fast defaults for sweeps -----
# These keep each run_demo_fast() call to ~0.1-0.3 seconds
# while preserving qualitative model behavior.

SWEEP_BASE_PARAMS <- list(
  n         = 25,
  k_min     = 1,
  k_shape   = 2.2,
  mu_R      = 10,
  sd_R      = 3,
  alpha_out = 0.12,
  eps       = 0.25,
  kappa     = 2.5,
  T         = 24,
  dt        = 0.5,
  tau       = 6,
  xi        = 3.0,
  B         = 60,
  delta     = 2,
  M         = 800,
  verbose   = FALSE
)

# ============================================================
# SWEEP CONFIGURATIONS
# ============================================================
# Each entry defines a self-contained research question.
# Adding a new sweep = adding one list entry here. No UI changes needed.

SWEEP_CONFIGS <- list(

  # ----------------------------------------------------------
  signal_value = list(
    name = "Signal Value",
    description = paste(
      "When does the proposal signal matter?",
      "This sweep varies the signal noise (\u03be) across levels of",
      "population heterogeneity (knowledge tail weight and resource spread).",
      "The key outcome is what fraction of the total funding gain",
      "is attributable to observing the proposal signal.",
      "We expect the signal to matter most when confounding is severe:",
      "heavy-tailed knowledge (low k_shape) combined with high resource",
      "variation (high sd_R)."
    ),
    grid_fn = function() {
      expand.grid(
        xi      = c(1, 2, 4, 8, 16),
        k_shape = c(1.5, 2.2, 4.0),
        sd_R    = c(1, 3, 8),
        stringsAsFactors = FALSE
      )
    },
    run_without_signal = TRUE,
    varied_params = c("xi", "k_shape", "sd_R"),
    primary_plot = list(
      type       = "heatmap",
      x_var      = "xi",
      y_var      = "k_shape",
      fill_var   = "gain_signal_pct_mean",
      facet_var  = "sd_R",
      title      = "Value of proposal signal (% of total gain)",
      fill_label = "Signal value\n(% of gain)",
      text_var   = "gain_signal_pct_mean",
      digits     = 1
    ),
    secondary_plot = list(
      type      = "line",
      x_var     = "xi",
      y_var     = "gain_signal_mean",
      color_var = "k_shape",
      facet_var = "sd_R",
      title     = "Absolute signal value across noise levels",
      y_label   = "Signal value (expected papers)"
    )
  ),

  # ----------------------------------------------------------
  naive_gap = list(
    name = "Naive vs. Optimal Gap",
    description = paste(
      "When does naive (publication-proportional) allocation suffice?",
      "This sweep varies the population structure: knowledge tail weight",
      "(k_shape) and resource spread (sd_R).",
      "When knowledge is homogeneous (high k_shape) or resources are",
      "homogeneous (low sd_R), the confounding problem is mild and",
      "naive allocation performs nearly as well as optimal.",
      "The gap should be largest when both K and R vary widely and",
      "independently."
    ),
    grid_fn = function() {
      expand.grid(
        k_shape = c(1.3, 1.5, 1.8, 2.2, 2.8, 3.5, 4.0),
        sd_R    = c(1, 2, 4, 6, 8),
        stringsAsFactors = FALSE
      )
    },
    run_without_signal = FALSE,
    varied_params = c("k_shape", "sd_R"),
    primary_plot = list(
      type       = "heatmap",
      x_var      = "k_shape",
      y_var      = "sd_R",
      fill_var   = "gain_over_naive_pct_mean",
      facet_var  = NULL,
      title      = "Advantage of optimal over naive allocation",
      fill_label = "Optimal advantage\n(% of total gain)",
      text_var   = "gain_over_naive_pct_mean",
      digits     = 1
    ),
    secondary_plot = list(
      type      = "line",
      x_var     = "k_shape",
      y_var     = "gain_over_naive_mean",
      color_var = "sd_R",
      facet_var = NULL,
      title     = "Absolute gain of optimal over naive",
      y_label   = "Optimal \u2212 Naive (expected papers)"
    )
  ),

  # ----------------------------------------------------------
  time_horizon = list(
    name = "Time Horizon Effects",
    description = paste(
      "Does a longer time horizon change who gets funded?",
      "Short horizons reward researchers who are productive now.",
      "Long horizons reward those with high growth potential \u2014",
      "high K0 with room to grow (cap = \u03baR far above K0).",
      "This sweep varies the time horizon (T) and the knowledge",
      "growth rate (\u03b5). We look for reversals: does optimal",
      "allocation shift from 'fund the currently productive' to",
      "'fund the high-potential' as T grows? Does naive allocation",
      "fail more at long horizons?"
    ),
    grid_fn = function() {
      expand.grid(
        T   = c(4, 8, 16, 32, 48),
        eps = c(0.05, 0.1, 0.25, 0.5, 0.8),
        stringsAsFactors = FALSE
      )
    },
    run_without_signal = FALSE,
    varied_params = c("T", "eps"),
    primary_plot = list(
      type       = "heatmap",
      x_var      = "T",
      y_var      = "eps",
      fill_var   = "gain_optimal_mean",
      facet_var  = NULL,
      title      = "Total gain from optimal funding",
      fill_label = "Optimal gain\n(expected papers)",
      text_var   = "gain_optimal_mean",
      digits     = 1
    ),
    secondary_plot = list(
      type      = "line",
      x_var     = "T",
      y_var     = "gain_over_naive_pct_mean",
      color_var = "eps",
      facet_var = NULL,
      title     = "Advantage of optimal over naive across horizons",
      y_label   = "Optimal advantage (% of total gain)"
    )
  ),

  # ----------------------------------------------------------
  knowledge_ceiling = list(
    name = "Knowledge Ceiling",
    description = paste(
      "Does the knowledge ceiling (\u03ba) reshape who gets funded?",
      "\u03ba determines how much knowledge a researcher can accumulate",
      "relative to resources (cap = \u03baR). Low \u03ba means funding quickly",
      "hits diminishing returns as knowledge caps out.",
      "High \u03ba means funded researchers keep growing.",
      "This maps to real-world questions: do fields have hard",
      "capacity constraints? We expect \u03ba to interact with T \u2014",
      "the ceiling only binds at long horizons."
    ),
    grid_fn = function() {
      expand.grid(
        kappa = c(0.8, 1.5, 2.5, 4.0, 6.0, 8.0),
        T     = c(8, 16, 24, 48),
        stringsAsFactors = FALSE
      )
    },
    run_without_signal = FALSE,
    varied_params = c("kappa", "T"),
    primary_plot = list(
      type       = "heatmap",
      x_var      = "kappa",
      y_var      = "T",
      fill_var   = "gain_optimal_mean",
      facet_var  = NULL,
      title      = "Total gain from optimal funding by ceiling and horizon",
      fill_label = "Optimal gain\n(expected papers)",
      text_var   = "gain_optimal_mean",
      digits     = 1
    ),
    secondary_plot = list(
      type      = "line",
      x_var     = "kappa",
      y_var     = "gain_over_naive_pct_mean",
      color_var = "T",
      facet_var = NULL,
      title     = "Advantage of optimal over naive by knowledge ceiling",
      y_label   = "Optimal advantage (% of total gain)"
    )
  )
)

# ============================================================
# SWEEP RUNNER
# ============================================================

run_sweep <- function(sweep_name,
                      seeds       = 1:5,
                      base_params = NULL,
                      progress_fn = NULL) {

  config <- SWEEP_CONFIGS[[sweep_name]]
  if (is.null(config)) stop("Unknown sweep: ", sweep_name)

  params_grid  <- config$grid_fn()
  bp           <- modifyList(SWEEP_BASE_PARAMS, base_params %||% list())
  varied_names <- config$varied_params
  total_iters  <- nrow(params_grid) * length(seeds)

  results <- vector("list", total_iters)
  counter <- 0L

  for (row_i in seq_len(nrow(params_grid))) {
    row_vals   <- as.list(params_grid[row_i, , drop = FALSE])
    run_params <- modifyList(bp, row_vals)

    for (s in seeds) {
      counter <- counter + 1L
      run_params$seed <- s

      # --- Primary run (with signal, or whatever xi is set to) ---
      res <- do.call(run_demo_fast, run_params)

      baseline_total <- res$eval$expected_total_baseline
      optimal_total  <- res$eval$expected_total_funded

      # --- Naive allocation ---
      g_naive <- allocate_naive_publications(res$obs$p, run_params$B)
      naive_total <- sum(expected_papers_EV_no_decay(
        res$truth$K0, res$truth$R0 + g_naive,
        run_params$T, run_params$dt, run_params$eps,
        run_params$kappa, run_params$alpha_out
      ))

      # --- Without-signal run (optional) ---
      optimal_no_signal_total <- NA_real_
      if (isTRUE(config$run_without_signal)) {
        run_params_ns      <- modifyList(run_params, list(xi = Inf))
        res_ns             <- do.call(run_demo_fast, run_params_ns)
        optimal_no_signal_total <- res_ns$eval$expected_total_funded
      }

      # --- Derived metrics ---
      gain_naive    <- naive_total - baseline_total
      gain_optimal  <- optimal_total - baseline_total
      gain_over_naive     <- optimal_total - naive_total
      gain_over_naive_pct <- 100 * gain_over_naive /
                             pmax(gain_optimal, 1e-6)
      gain_signal     <- optimal_total - optimal_no_signal_total
      gain_signal_pct <- if (!is.na(optimal_no_signal_total)) {
        100 * gain_signal / pmax(gain_optimal, 1e-6)
      } else {
        NA_real_
      }

      # --- Assemble row ---
      result_row <- c(
        row_vals,
        list(
          seed                    = s,
          baseline_total          = baseline_total,
          naive_total             = naive_total,
          optimal_total           = optimal_total,
          optimal_no_signal_total = optimal_no_signal_total,
          gain_naive              = gain_naive,
          gain_optimal            = gain_optimal,
          gain_over_naive         = gain_over_naive,
          gain_over_naive_pct     = gain_over_naive_pct,
          gain_signal             = gain_signal,
          gain_signal_pct         = gain_signal_pct
        )
      )
      results[[counter]] <- as.data.frame(result_row, stringsAsFactors = FALSE)

      if (!is.null(progress_fn)) progress_fn(counter, total_iters)
    }
  }

  do.call(rbind, results)
}

# ============================================================
# SUMMARIZE SWEEP
# ============================================================

summarize_sweep <- function(raw_df, config) {
  varied <- config$varied_params
  outcome_cols <- c("baseline_total", "naive_total", "optimal_total",
                    "optimal_no_signal_total",
                    "gain_naive", "gain_optimal",
                    "gain_over_naive", "gain_over_naive_pct",
                    "gain_signal", "gain_signal_pct")

  # Group by varied parameters, compute mean and SE for each outcome
  raw_df %>%
    group_by(across(all_of(varied))) %>%
    summarise(
      n_seeds = n(),
      across(
        all_of(outcome_cols),
        list(
          mean = ~mean(.x, na.rm = TRUE),
          se   = ~sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

plot_sweep_heatmap <- function(summary_df, x_var, y_var, fill_var,
                               facet_var = NULL, title = "",
                               fill_label = "", text_var = NULL,
                               digits = 1) {
  # Convert axes to factors for discrete tiles
  df <- summary_df %>%
    mutate(
      across(all_of(c(x_var, y_var)), as.factor)
    )
  if (!is.null(facet_var)) {
    df <- df %>% mutate(across(all_of(facet_var), ~paste0(facet_var, " = ", .)))
  }

  p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]],
                       fill = .data[[fill_var]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_continuous_safe(name = fill_label) +
    labs(title = title, x = x_var, y = y_var) +
    theme_sim(base_size = 12) +
    theme(
      axis.text    = element_text(size = 10),
      legend.position = "right"
    )

  # Add text labels
  if (!is.null(text_var) && text_var %in% names(df)) {
    p <- p + geom_text(
      aes(label = sprintf(paste0("%.", digits, "f"), .data[[text_var]])),
      size = 3.2, color = "white", fontface = "bold"
    )
  }

  # Facet if needed
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)), nrow = 1)
  }

  p
}

plot_sweep_line <- function(summary_df, x_var, y_var, color_var,
                            facet_var = NULL, title = "",
                            y_label = "") {
  # SE variable name
  se_var <- sub("_mean$", "_se", y_var)
  has_se <- se_var %in% names(summary_df)

  df <- summary_df %>%
    mutate(across(all_of(color_var), as.factor))

  p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]],
                       color = .data[[color_var]],
                       fill  = .data[[color_var]],
                       group = .data[[color_var]]))

  # SE ribbon

  if (has_se) {
    p <- p + geom_ribbon(
      aes(ymin = .data[[y_var]] - .data[[se_var]],
          ymax = .data[[y_var]] + .data[[se_var]]),
      alpha = 0.15, color = NA
    )
  }

  p <- p +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    labs(title = title, x = x_var, y = y_label,
         color = color_var, fill = color_var) +
    theme_sim(base_size = 12) +
    theme(legend.position = "bottom")

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  }

  p
}

# Dispatch: read plot spec from config and call the appropriate function
plot_sweep_from_config <- function(summary_df, config, which = "primary") {
  spec <- if (which == "primary") config$primary_plot else config$secondary_plot
  if (is.null(spec)) return(NULL)

  if (spec$type == "heatmap") {
    plot_sweep_heatmap(
      summary_df,
      x_var      = spec$x_var,
      y_var      = spec$y_var,
      fill_var   = spec$fill_var,
      facet_var  = spec$facet_var,
      title      = spec$title,
      fill_label = spec$fill_label,
      text_var   = spec$text_var,
      digits     = spec$digits %||% 1
    )
  } else if (spec$type == "line") {
    plot_sweep_line(
      summary_df,
      x_var     = spec$x_var,
      y_var     = spec$y_var,
      color_var = spec$color_var,
      facet_var = spec$facet_var,
      title     = spec$title,
      y_label   = spec$y_label
    )
  }
}
