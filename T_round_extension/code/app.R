# app.R
# ============================================================
# Grant-funding / scientific-output model — Shiny explorer
# ------------------------------------------------------------
# Interactive front-end for the model. ALL simulation logic lives in
# model.R (sourced below); this file is the presentation layer only:
# the UI controls, the server that calls run_simulation_T(), and the plots.
#
# The model is the T-round Bayesian grant-allocation model (T in 1..5):
# a funder allocates a fixed budget across n researchers and T rounds under
# 9 strategies (no-funding -> naive -> myopic -> forward-looking CE planner),
# inferring each researcher's knowledge/resources from noisy publications and
# two signals. Knowledge compounds between rounds. The classic case is T = 2.
#
# See README.md for the project map, and T_round_extension/ for the analysis
# package (data, validation, results) and its STATE_OF_PLAY.md history.
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
# MODEL ENGINE
# ============================================================
# All simulation logic lives in model.R (pure base R): it defines
# STRATEGY_NAMES, run_simulation() (the 2-round v5 reference) and
# run_simulation_T() (the default T-round model). This file (app.R)
# is the Shiny presentation layer only.
source("model.R")

# ============================================================
# PRESENTATION CONSTANTS
# ============================================================
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
           T sequential funding rounds when researchers differ in knowledge
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
                        sliderInput("T_rounds", "Funding rounds (T)", min = 1, max = 5, value = 2, step = 1),
                        div(class = "param-help", "Number of sequential grant rounds. The total budget is split evenly across the T rounds, and knowledge compounds between them. T = 2 is the classic case; larger T lets the forward-looking planner front-load to compound knowledge."),

                        sliderInput("n", "Number of researchers", min = 10, max = 500, value = 100, step = 10),
                        div(class = "param-help", "Population size."),

                        sliderInput("b", "Budget fraction (b = B / n·E[R])", min = 0.0, max = 2.0, value = 0.5, step = 0.05),
                        div(class = "param-help", "Sets the total budget \\(B_{total} = 2\\,b\\,n\\,E[R]\\), split evenly across the T rounds; b scales it relative to the population's aggregate expected resource."),
                        
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
                        div(class = "param-help", "Rate of knowledge growth between rounds. After each round K grows in proportion to that round's output rate: \\(K_{t+1} = K_t + \\varepsilon \\cdot K_t R_t / (K_t + R_t)\\).")
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
                  (sum of \\(\\lambda\\) over all rounds and all researchers) for
                  each selected strategy. Higher bars = more total output.")
                           ),
                           plotOutput("fig_strategy_comparison", height = 520)
                  ),
                  
                  tabPanel("Funding Effects",
                           div(class = "plot-explanation",
                               HTML("<strong>What this shows:</strong> How the T funding rounds
                  reshape the \\((K, R)\\) distribution for a selected strategy. The first
                  panel is the state entering round 1 (post-pre-rounds); each subsequent
                  panel is the effective \\((K, R)\\) during that round — knowledge grows
                  round to round, and \\(R\\) includes that round's grant. Arrows on the
                  final panel show each researcher's trajectory from start to the last round.")
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
                  measures across researchers at initial conditions and after each of
                  the T rounds (grey → deepening blue). Select a strategy to view.")
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
    updateSliderInput(session, "T_rounds", value = 2)
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
          all_runs[[i]] <- run_simulation_T(
            seed = base_seed + (i - 1L),
            T_rounds = input$T_rounds,
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
            strategies = strats,
            detail = (i == 1L)   # trajectories only for run 1 (the state plots use it)
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
            bidxs <- vapply(all_runs, function(r) {
              v <- r$strategies[[s]]$b_idx;            if (is.null(v)) NA_real_ else v
            }, numeric(1))

            res$strategies[[s]]$total_output     <- mean(outs, na.rm = TRUE)
            res$strategies[[s]]$se_total_output  <- sd(outs,   na.rm = TRUE) / sqrt(sum(!is.na(outs)))
            res$strategies[[s]]$total_expected   <- mean(exps, na.rm = TRUE)
            res$strategies[[s]]$se_total_expected<- sd(exps,   na.rm = TRUE) / sqrt(sum(!is.na(exps)))
            res$strategies[[s]]$alpha            <- mean(alphas, na.rm = TRUE)
            res$strategies[[s]]$se_alpha         <- sd(alphas, na.rm = TRUE) / sqrt(sum(!is.na(alphas)))
            res$strategies[[s]]$b_idx            <- mean(bidxs, na.rm = TRUE)
            res$strategies[[s]]$se_b_idx         <- sd(bidxs, na.rm = TRUE) / sqrt(sum(!is.na(bidxs)))
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
      cat(sprintf("Trials: %d (seeds %d–%d) | T = %d rounds | n = %d | b = %.2f (B_total = %.1f) | Pre-rounds: %d\n",
                  res$n_trials, res$seed_range[1], res$seed_range[2],
                  res$params$T_rounds, res$params$n, res$params$b,
                  res$params$B_total, res$params$n_pre_rounds))
    } else {
      cat(sprintf("Seed: %d | T = %d rounds | n = %d | b = %.2f (B_total = %.1f) | Pre-rounds: %d\n",
                  res$params$seed, res$params$T_rounds, res$params$n, res$params$b,
                  res$params$B_total, res$params$n_pre_rounds))
    }
    cat(sprintf("K ~ Pareto(min=%.1f, shape=%.1f) | R ~ Pareto(min=%.1f, shape=%.1f) | rho = %.2f\n",
                res$params$k_min, res$params$k_shape,
                res$params$r_min, res$params$r_shape,
                res$params$rho_kr))
    cat(sprintf("K-R correlation: initial = %.3f, post-pre-rounds = %.3f\n\n",
                cor(res$initial_state$K, res$initial_state$R),
                cor(res$post_preround_state$K, res$post_preround_state$R)))
    
    if (multi) {
      cat(sprintf("  %-30s %16s   %5s  %6s\n", "Strategy", "Output (mean±SE)", "alpha", "b_idx"))
    } else {
      cat(sprintf("  %-30s %8s   %5s  %6s\n", "Strategy", "Output", "alpha", "b_idx"))
    }
    for (s in seq_along(res$strategies)) {
      r <- res$strategies[[s]]
      if (is.null(r)) next
      alpha_str <- if (is.na(r$alpha)) "  —  " else sprintf("%.3f", r$alpha)
      bidx_str  <- if (is.null(r$b_idx) || is.na(r$b_idx)) "  —  " else sprintf("%.3f", r$b_idx)
      if (multi) {
        se <- if (is.null(r$se_total_expected) || is.na(r$se_total_expected)) 0 else r$se_total_expected
        cat(sprintf("  %-30s %8.1f ± %4.2f   %5s  %6s\n", r$name, r$total_expected, se, alpha_str, bidx_str))
      } else {
        cat(sprintf("  %-30s %8.1f   %5s  %6s\n", r$name, r$total_expected, alpha_str, bidx_str))
      }
    }
    cat("\n  alpha: round-1 share of total spend.   b_idx: schedule center-of-mass\n")
    cat("  (0.5 = even across rounds, >0.5 = front-loaded). Equal-tranche strategies give b_idx = 0.5.\n")
  })
  
  output$fig_strategy_comparison <- renderPlot({
    res <- sim_result()
    req(res)
    
    multi <- !is.null(res$n_trials) && res$n_trials > 1
    
    df <- tibble(Strategy = character(), Total = numeric(), SE = numeric())
    for (s in seq_along(res$strategies)) {
      r <- res$strategies[[s]]
      if (is.null(r)) next
      se_val <- if (!is.null(r$se_total_expected) && !is.na(r$se_total_expected)) r$se_total_expected else 0
      df <- bind_rows(df, tibble(Strategy = r$name, Total = r$total_expected, SE = se_val))
    }
    df$Strategy <- factor(df$Strategy, levels = df$Strategy)
    
    baseline_val <- df$Total[df$Strategy == "No funding"]
    if (length(baseline_val) == 0) baseline_val <- 0
    
    subtitle_text <- if (multi) {
      sprintf("T = %d, n = %d, b = %.2f (B_total = %.1f), gamma = %.1f | mean ± SE over %d trials",
              res$params$T_rounds, res$params$n, res$params$b, res$params$B_total, res$params$gamma, res$n_trials)
    } else {
      sprintf("T = %d, n = %d, b = %.2f (B_total = %.1f), gamma = %.1f",
              res$params$T_rounds, res$params$n, res$params$b, res$params$B_total, res$params$gamma)
    }
    
    p <- ggplot(df, aes(x = Strategy, y = Total, fill = Strategy)) +
      geom_col(width = 0.7, alpha = 0.9) +
      {if (baseline_val > 0) geom_hline(yintercept = baseline_val, linetype = "dashed", alpha = 0.5)}
    
    if (multi) {
      p <- p + geom_errorbar(aes(ymin = Total - SE, ymax = Total + SE),
                             width = 0.25, alpha = 0.7, linewidth = 0.5)
    }
    
    p +
      geom_text(aes(label = sprintf("%.1f", Total)),
                vjust = -0.5, size = 3.5, fontface = "bold",
                nudge_y = if (multi) max(df$SE) * 1.1 else 0) +
      scale_fill_manual(values = STRATEGY_COLORS, drop = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(
        title = sprintf("Total expected output by strategy (%d round%s)",
                        res$params$T_rounds, if (res$params$T_rounds == 1) "" else "s"),
        subtitle = subtitle_text,
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
  
  output$fig_preround_effects <- renderPlot({
    res <- sim_result()
    req(res)
    
    df_init <- tibble(
      K = res$initial_state$K,
      R = res$initial_state$R,
      D = (res$initial_state$K - res$initial_state$R) /
        (res$initial_state$K + res$initial_state$R),
      Phase = "Initial (raw draws)"
    )
    cor_init <- cor(res$initial_state$K, res$initial_state$R)
    
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
    
    max_val <- max(c(df_all$K, df_all$R), na.rm = TRUE) * 1.05
    
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
  
  output$fig_funding_effects <- renderPlot({
    res <- sim_result()
    req(res)
    
    s_idx <- as.integer(input$fe_strategy)
    strat <- res$strategies[[s_idx]]
    req(strat, !is.null(strat$K_rounds))

    n  <- res$params$n
    Tr <- res$params$T_rounds
    stage_levels <- c("Before funding", sprintf("Round %d", seq_len(Tr)))

    df_start <- tibble(
      K = res$K_at_start,
      R = res$R0_at_start,
      D = (res$K_at_start - res$R0_at_start) /
        (res$K_at_start + res$R0_at_start),
      id = seq_len(n),
      Stage = "Before funding"
    )

    # one panel per round: knowledge entering round t vs resources that round
    df_rounds <- bind_rows(lapply(seq_len(Tr), function(t) {
      Kt <- strat$K_rounds[[t]]; Rt <- strat$R_rounds[[t]]
      tibble(K = Kt, R = Rt, D = (Kt - Rt) / (Kt + Rt),
             id = seq_len(n), Stage = sprintf("Round %d", t))
    }))

    df_all <- bind_rows(df_start, df_rounds)
    df_all$Stage <- factor(df_all$Stage, levels = stage_levels)

    max_val <- max(c(df_all$K, df_all$R), na.rm = TRUE) * 1.05

    ann <- tibble(
      Stage = factor(stage_levels, levels = stage_levels),
      label = c(sprintf("r = %.3f", cor(res$K_at_start, res$R0_at_start)),
                vapply(seq_len(Tr), function(t)
                  sprintf("r = %.3f", cor(strat$K_rounds[[t]], strat$R_rounds[[t]])),
                  character(1)))
    )

    # arrows: initial state -> final round, drawn on the last-round facet
    df_arrows <- tibble(
      x_start = res$R0_at_start,
      y_start = res$K_at_start,
      x_end = strat$R_rounds[[Tr]],
      y_end = strat$K_rounds[[Tr]],
      moved = abs(strat$R_rounds[[Tr]] - res$R0_at_start) > 0.01 |
              abs(strat$K_rounds[[Tr]] - res$K_at_start) > 0.01,
      Stage = factor(sprintf("Round %d", Tr), levels = stage_levels)
    ) %>% filter(moved)

    ggplot(df_all, aes(x = R, y = K, color = D)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.4) +
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
                inherit.aes = FALSE, hjust = 0, size = 3.4, fontface = "bold") +
      facet_wrap(~ Stage, ncol = 3) +
      coord_cartesian(xlim = c(0, max_val), ylim = c(0, max_val)) +
      labs(
        title = sprintf("Funding effects on (K, R) distribution: %s", strat$name),
        subtitle = sprintf("T = %d | b = %.2f (B_total = %.1f) | alpha = %s | Arrows: initial state → final round",
                           Tr, res$params$b, res$params$B_total,
                           if (is.na(strat$alpha)) "NA" else sprintf("%.3f", strat$alpha)),
        x = "Resources (R)",
        y = "Knowledge (K)"
      ) +
      theme_sim(base_size = 12) +
      theme(plot.subtitle = element_text(color = "grey40", size = 10))
  }, res = 110)
  
  output$fig_bottleneck <- renderPlot({
    res <- sim_result()
    req(res)
    
    s_idx <- as.integer(input$bn_strategy)
    strat <- res$strategies[[s_idx]]
    req(strat, !is.null(strat$K_rounds))

    Tr <- res$params$T_rounds
    stage_levels <- c("Initial", sprintf("After Round %d", seq_len(Tr)))

    bn0 <- compute_bottleneck(res$K_at_start, res$R0_at_start)
    df_bn <- bind_rows(
      tibble(D = bn0$D, S = bn0$S, Stage = "Initial"),
      bind_rows(lapply(seq_len(Tr), function(t) {
        bn <- compute_bottleneck(strat$K_rounds[[t]], strat$R_rounds[[t]])
        tibble(D = bn$D, S = bn$S, Stage = sprintf("After Round %d", t))
      }))
    )
    df_bn$Stage <- factor(df_bn$Stage, levels = stage_levels)

    df_d <- df_bn %>% select(Stage, value = D) %>% mutate(Measure = "Direction (D)")
    df_s <- df_bn %>% select(Stage, value = S) %>% mutate(Measure = "Severity (S)")

    df_long <- bind_rows(df_d, df_s)
    df_long$Measure <- factor(df_long$Measure,
                              levels = c("Direction (D)", "Severity (S)"))

    # colour ramp: grey Initial → deepening blue across the rounds
    stage_cols <- setNames(c("#bdbdbd", colorRampPalette(c("#90caf9", "#0d47a1"))(Tr)),
                           stage_levels)

    ggplot(df_long, aes(x = value, fill = Stage)) +
      geom_density(alpha = 0.4, linewidth = 0.4) +
      facet_wrap(~ Measure, scales = "free") +
      scale_fill_manual(values = stage_cols) +
      labs(
        title = sprintf("Bottleneck measures: %s", strat$name),
        x = "Value",
        y = "Density",
        fill = "Stage"
      ) +
      theme_sim(base_size = 12) +
      theme(legend.position = "bottom")
  }, res = 110)
  
  output$fig_signal_value <- renderPlot({
    res <- sim_result()
    req(res)
    
    comparisons <- tibble(
      Comparison = character(),
      Without = numeric(), With = numeric(), Gain = numeric(),
      SE_Without = numeric(), SE_With = numeric()
    )
    
    s4 <- res$strategies[[4]]; s5 <- res$strategies[[5]]
    if (!is.null(s4) && !is.null(s5)) {
      comparisons <- bind_rows(comparisons, tibble(
        Comparison = "Myopic:\ngrant signal",
        Without = s4$total_expected, With = s5$total_expected,
        Gain = s5$total_expected - s4$total_expected,
        SE_Without = se_or_zero(s4), SE_With = se_or_zero(s5)
      ))
    }
    
    s7 <- res$strategies[[7]]; s8 <- res$strategies[[8]]
    if (!is.null(s7) && !is.null(s8)) {
      comparisons <- bind_rows(comparisons, tibble(
        Comparison = "Forward:\ngrant signal",
        Without = s7$total_expected, With = s8$total_expected,
        Gain = s8$total_expected - s7$total_expected,
        SE_Without = se_or_zero(s7), SE_With = se_or_zero(s8)
      ))
    }
    
    if (nrow(comparisons) == 0) {
      plot.new()
      text(0.5, 0.5, "Run strategies 4, 5, 7, 8 to see signal value comparisons",
           cex = 1.2)
      return(invisible(NULL))
    }
    
    comparisons$Comparison <- factor(comparisons$Comparison, levels = comparisons$Comparison)
    make_value_plot(
      comparisons,
      title = "Value of grant signals",
      subtitle = "Gain from adding grant signal to pubs-only strategies",
      baseline_label = "Without signal", treatment_label = "With signal",
      baseline_color = "#90caf9",         treatment_color = "#1565c0"
    )
  }, res = 110)
  
  output$fig_forward_value <- renderPlot({
    res <- sim_result()
    req(res)
    
    comparisons <- tibble(
      Comparison = character(),
      Without = numeric(), With = numeric(), Gain = numeric(),
      SE_Without = numeric(), SE_With = numeric()
    )
    
    # Pairs at each information / intervention setting
    pairs <- list(
      list(label = "Pubs only",      m = 4, f = 7),
      list(label = "Pubs + grant",   m = 5, f = 8),
      list(label = "Pubs + seed",    m = 6, f = 9)
    )
    for (p in pairs) {
      sm <- res$strategies[[p$m]]; sf <- res$strategies[[p$f]]
      if (!is.null(sm) && !is.null(sf)) {
        comparisons <- bind_rows(comparisons, tibble(
          Comparison = p$label,
          Without = sm$total_expected, With = sf$total_expected,
          Gain = sf$total_expected - sm$total_expected,
          SE_Without = se_or_zero(sm), SE_With = se_or_zero(sf)
        ))
      }
    }
    
    if (nrow(comparisons) == 0) {
      plot.new()
      text(0.5, 0.5, "Run at least one matched (myopic, forward) pair: 4↔7, 5↔8, or 6↔9",
           cex = 1.1)
      return(invisible(NULL))
    }
    
    comparisons$Comparison <- factor(comparisons$Comparison, levels = comparisons$Comparison)
    make_value_plot(
      comparisons,
      title    = "Value of forward-looking planning",
      subtitle = "Gain from upgrading myopic to forward (CE) planner, holding signal and seed fixed",
      baseline_label = "Myopic",  treatment_label = "Forward",
      baseline_color = "#1e88e5", treatment_color = "#c62828"
    )
  }, res = 110)
  
  output$fig_seed_value <- renderPlot({
    res <- sim_result()
    req(res)
    
    comparisons <- tibble(
      Comparison = character(),
      Without = numeric(), With = numeric(), Gain = numeric(),
      SE_Without = numeric(), SE_With = numeric()
    )
    
    # Pairs at each planner type (no grant signal in either side, so the
    # comparison isolates the seed intervention)
    pairs <- list(
      list(label = "Myopic",  no = 4, ye = 6),
      list(label = "Forward", no = 7, ye = 9)
    )
    for (p in pairs) {
      s_no <- res$strategies[[p$no]]; s_ye <- res$strategies[[p$ye]]
      if (!is.null(s_no) && !is.null(s_ye)) {
        comparisons <- bind_rows(comparisons, tibble(
          Comparison = p$label,
          Without = s_no$total_expected, With = s_ye$total_expected,
          Gain = s_ye$total_expected - s_no$total_expected,
          SE_Without = se_or_zero(s_no), SE_With = se_or_zero(s_ye)
        ))
      }
    }
    
    if (nrow(comparisons) == 0) {
      plot.new()
      text(0.5, 0.5, "Run at least one matched (no-seed, seed) pair: 4↔6 or 7↔9",
           cex = 1.1)
      return(invisible(NULL))
    }
    
    comparisons$Comparison <- factor(comparisons$Comparison, levels = comparisons$Comparison)
    make_value_plot(
      comparisons,
      title    = "Value of seed grants",
      subtitle = "Gain from adding uniform seed to pubs-only strategies (grant signal off in both)",
      baseline_label = "No seed",  treatment_label = "With seed",
      baseline_color = "#42a5f5",  treatment_color = "#7e57c2"
    )
  }, res = 110)
}

shinyApp(ui, server)
