# Plot the resource_regime reversal map. Run after the sweep finishes.
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
indir  <- "sweep_results/T_run_smooth_supplement"
outdir <- "figures/resource_regime"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
s <- readRDS(file.path(indir, "resource_regime_summary.rds"))

# round-1 share of forward planner (alpha_S8_mean) if present
has_r1 <- "alpha_S8_mean" %in% names(s)
s$eps_lab <- factor(sprintf("%g", s$epsilon), levels = sprintf("%g", sort(unique(s$epsilon))))
s$rmin_lab <- factor(sprintf("%g", s$r_min),  levels = sprintf("%g", sort(unique(s$r_min))))

# ---- (1) Heatmap of b_idx_S8 with 0.5 boundary ----
p1 <- ggplot(s, aes(eps_lab, rmin_lab, fill = b_idx_S8_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", b_idx_S8_mean),
                color = b_idx_S8_mean < 0.5), size = 3) +
  scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_gradient2(midpoint = 0.5, low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       name = "b_idx (S8)") +
  labs(x = expression(epsilon~"(compounding rate)"),
       y = expression(r[min]~"(baseline resources)"),
       title = "Forward-planner schedule: front-load (blue, <0.5) vs back-load (red, >0.5)",
       subtitle = "Decoupled purse (budget_ref=K), T=5, smooth allocator. Boundary = b_idx=0.5.") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "resource_regime_heatmap.png"), p1, width = 8.5, height = 5, dpi = 130)

# ---- (2) Line: b_idx_S8 vs epsilon by r_min, with 0.5 line ----
p2 <- ggplot(s, aes(epsilon, b_idx_S8_mean, color = rmin_lab, group = rmin_lab)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_x_log10() +
  labs(x = expression(epsilon~"(log scale)"), y = "b_idx (S8): schedule center-of-mass",
       color = expression(r[min]),
       title = "Reversal is governed by compounding, not baseline resources",
       subtitle = "b_idx<0.5 front-load; >0.5 back-load. Curves shift up with r_min but the crossing is set by epsilon.") +
  annotate("text", x = min(s$epsilon), y = 0.503, label = "even split", hjust = 0, size = 3, color = "grey40") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "resource_regime_bidx_vs_eps.png"), p2, width = 8.5, height = 5, dpi = 130)

# ---- (3) Round-1 share vs r_min at each epsilon (shows scarcity pushes spend forward) ----
if (has_r1) {
  p3 <- ggplot(s, aes(r_min, alpha_S8_mean, color = eps_lab, group = eps_lab)) +
    geom_hline(yintercept = 1/5, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 0.8) + geom_point(size = 2) + scale_x_log10() +
    labs(x = expression(r[min]~"(log scale)"), y = "round-1 budget share (S8)",
         color = expression(epsilon),
         title = "Scarcer baseline resources push spending toward round 1",
         subtitle = "Even split = 1/T = 0.2 (dashed). Directional shift, but rarely crosses into strict front-loading.") +
    theme_minimal(base_size = 12)
  ggsave(file.path(outdir, "resource_regime_r1share.png"), p3, width = 8.5, height = 5, dpi = 130)
}

# ---- console summary + the 0.5 crossing per r_min ----
cat("\n=== resource_regime summary (b_idx_S8_mean) ===\n")
tab <- s %>% select(r_min, epsilon, b_idx_S8_mean,
                    any_of(c("alpha_S8_mean","fwd_vs_myo_PG_mean","signal_fwd_mean"))) %>%
  arrange(r_min, epsilon)
print(as.data.frame(tab), digits = 4, row.names = FALSE)

cat("\n=== front-load cells (b_idx_S8 < 0.5) ===\n")
fl <- tab %>% filter(b_idx_S8_mean < 0.5)
if (nrow(fl)) print(as.data.frame(fl), digits = 4, row.names = FALSE) else cat("(none)\n")

cat("\nfigures written to", outdir, "\n")
