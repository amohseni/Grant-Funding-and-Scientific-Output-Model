# figures/make_figures.R
# ============================================================
# Regenerates the publication figure set from the canonical sweep data
# (T_round_extension/data/). Each figure is written as both PNG (300 dpi)
# and PDF into figures/. Run from the project root:  Rscript figures/make_figures.R
# ============================================================
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr); library(scales)})

DATA <- "T_round_extension/data"
OUT  <- "figures"; dir.create(OUT, showWarnings = FALSE)
G    <- function(nm) as.data.frame(readRDS(file.path(DATA, paste0(nm, "_summary.rds"))))
have <- function(nm) file.exists(file.path(DATA, paste0(nm, "_summary.rds")))

# --- shared theme + palette -------------------------------------------------
theme_pub <- function(base = 12) theme_minimal(base_size = base) +
  theme(plot.title = element_text(face = "bold"),
        plot.title.position = "plot",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
        strip.text = element_text(face = "bold"),
        legend.position = "right")
TEAL <- "#12727c"; RUST <- "#bb4a30"
save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
  cat("  wrote", name, "(png+pdf)\n")
}
heat <- function(df, x, y, fill, lab, title, subtitle = NULL, digits = 2, diverging = FALSE) {
  df[[x]] <- factor(df[[x]]); df[[y]] <- factor(df[[y]])
  p <- ggplot(df, aes(.data[[x]], .data[[y]], fill = .data[[fill]])) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf(paste0("%.", digits, "f"), .data[[fill]])),
              size = 3.3, colour = "grey15") +
    labs(title = title, subtitle = subtitle, x = x, y = y, fill = lab) + theme_pub()
  if (diverging) p + scale_fill_gradient2(low = RUST, mid = "grey92", high = TEAL, midpoint = 0)
  else p + scale_fill_gradient(low = "#eaf3f2", high = "#0d5c64")
}

cat("Generating figures ->", OUT, "\n")

# --- Fig 1: horizon phase diagram (forward vs myopic) -----------------------
if (have("horizon_growth")) {
  s <- G("horizon_growth")
  save_fig(heat(s, "T_rounds", "epsilon", "fwd_vs_myo_PG_mean", "S8 - S5",
                "Forward vs. myopic advantage across horizon and growth",
                "fwd_vs_myo_PG = S8 - S5 (expected output); positive = forward wins",
                digits = 2, diverging = TRUE),
           "fig1_horizon_phase_diagram", 7.5, 5)
}

# --- Fig 2: signal-value law (inequality x precision) -----------------------
if (have("signal_value")) {
  s <- G("signal_value")
  save_fig(heat(s, "tau_k", "k_shape", "signal_fwd_mean", "S8 - S7",
                "Grant-signal value by inequality and precision",
                "signal_fwd = S8 - S7; k_shape smaller = heavier knowledge tail; tau_k smaller = sharper signal",
                digits = 1),
           "fig2_signal_value_law", 7, 5)
}

# --- Fig 3: signal precision decay ------------------------------------------
if (have("signal_precision")) {
  s <- G("signal_precision") |> arrange(tau_k)
  p <- ggplot(s, aes(tau_k, signal_fwd_mean)) +
    geom_ribbon(aes(ymin = signal_fwd_mean - signal_fwd_se, ymax = signal_fwd_mean + signal_fwd_se),
                fill = TEAL, alpha = 0.15) +
    geom_line(colour = TEAL, linewidth = 0.9) + geom_point(colour = TEAL, size = 2) +
    scale_x_log10() +
    labs(title = "Signal value decays with signal noise",
         subtitle = "signal_fwd = S8 - S7 vs grant-signal noise tau_k (log scale)",
         x = "tau_k (signal noise, log)", y = "S8 - S7") + theme_pub()
  save_fig(p, "fig3_signal_precision", 7, 4.5)
}

# --- Fig 4: correlation robustness ------------------------------------------
if (have("correlation")) {
  s <- G("correlation")
  p <- ggplot(s, aes(factor(rho_kr), signal_fwd_mean, fill = factor(k_shape))) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = signal_fwd_mean - signal_fwd_se, ymax = signal_fwd_mean + signal_fwd_se),
                  position = position_dodge(0.8), width = 0.2) +
    scale_fill_manual(values = c("#0a4f57", "#2c8a8f", "#95c2c1"), name = "k_shape (tail)") +
    labs(title = "Signal value survives K–R correlation",
         subtitle = "signal_fwd = S8 - S7 vs Gaussian-copula correlation rho_c (tau_k = 0.3)",
         x = "rho_c (K–R correlation)", y = "S8 - S7") + theme_pub()
  save_fig(p, "fig4_correlation_robustness", 7, 4.5)
}

# --- Fig 5: seeding never helps ---------------------------------------------
if (have("seed_value")) {
  s <- G("seed_value")
  save_fig(heat(s, "x_seed", "b", "seed_fwd_mean", "S9 - S7",
                "Uniform seed floors never help",
                "seed_fwd = S9 - S7; negative everywhere", digits = 2, diverging = TRUE),
           "fig5_seed_value", 7, 5)
}

# --- Fig 6: tail_map (knowledge x resource inequality) ----------------------
if (have("tail_map")) {
  s <- G("tail_map")
  save_fig(heat(s, "r_shape", "k_shape", "signal_fwd_mean", "S8 - S7",
                "Signal value by knowledge and resource inequality",
                "signal_fwd = S8 - S7 (tau_k = 0.3); smaller shape = heavier tail. NOTE: budget co-varies with r_shape",
                digits = 1),
           "fig6_tail_map", 7, 5)
}

# --- Fig 7: information channel (eps -> 0) -----------------------------------
if (have("info_value")) {
  s <- G("info_value")
  p <- ggplot(s, aes(T_rounds, fwd_vs_myo_PG_mean, colour = factor(tau_k))) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    geom_ribbon(aes(ymin = fwd_vs_myo_PG_mean - fwd_vs_myo_PG_se,
                    ymax = fwd_vs_myo_PG_mean + fwd_vs_myo_PG_se, fill = factor(tau_k)),
                alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9) + geom_point(size = 2) +
    labs(title = "Information channel isolated (no compounding, eps -> 0)",
         subtitle = "fwd_vs_myo_PG = S8 - S5 at eps ~ 0: forward's residual edge is pure information value",
         x = "Horizon T", y = "S8 - S5 (information value)", colour = "tau_k", fill = "tau_k") +
    theme_pub()
  save_fig(p, "fig7_information_value", 7, 4.5)
}

# --- Fig 8: horizon_long saturation (corrected granularity) -----------------
if (have("horizon_long")) {
  s <- G("horizon_long")
  p <- ggplot(s, aes(T_rounds, fwd_vs_myo_PG_mean, colour = factor(epsilon))) +
    geom_ribbon(aes(ymin = fwd_vs_myo_PG_mean - fwd_vs_myo_PG_se,
                    ymax = fwd_vs_myo_PG_mean + fwd_vs_myo_PG_se, fill = factor(epsilon)),
                alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.9) + geom_point(size = 2) +
    labs(title = "Forward advantage at long horizons",
         subtitle = "fwd_vs_myo_PG = S8 - S5, T = 5..10 (re-run at converged granularity)",
         x = "Horizon T", y = "S8 - S5", colour = "epsilon", fill = "epsilon") + theme_pub()
  save_fig(p, "fig8_horizon_long", 7, 4.5)
} else cat("  (skipping fig8: horizon_long not yet in canonical data)\n")

cat("Done.\n")
