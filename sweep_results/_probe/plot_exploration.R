# Figures + console digest for the three exploration sweeps (paid-information regime).
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
indir  <- "sweep_results/T_run_smooth_supplement"
outdir <- "figures/resource_regime"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

digest <- function(nm, keys) {
  s <- readRDS(file.path(indir, paste0(nm, "_summary.rds")))
  s$vs_unif <- s$out_S8_mean - s$out_S2_mean
  d <- s[do.call(order, s[keys]), c(keys, "b_idx_S8_mean","alpha_S8_mean","fwd_vs_myo_PG_mean","vs_unif")]
  names(d) <- c(keys, "b_idx","r1_share","PG","vs_unif")
  cat(sprintf("\n===== %s =====\n", nm)); print(as.data.frame(d), digits = 3, row.names = FALSE)
  s
}

sc <- digest("exploration_corner",  c("tau_k","epsilon"))
sp <- digest("exploration_poverty", c("tau_k","r_min"))
sd_ <- digest("exploration_depth",  c("tau_k","b"))

# ---- Depth figure: the decisive test ----
d <- sd_ %>% mutate(vs_unif = out_S8_mean - out_S2_mean,
                    signal = factor(ifelse(tau_k < 10, "free signal on (tau_k=1)", "no free signal (tau_k=100)")))
p1 <- ggplot(d, aes(b, b_idx_S8_mean, color = signal)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  labs(x = "budget scale b (per-capita depth = 2bE[K]/T)", y = "b_idx (S8)",
       title = "Does funding depth tilt the schedule early? (poverty, eps~0, T=6)",
       subtitle = "b_idx<0.5 = money front-loaded; 0.5 = even; >0.5 = back-loaded (seed-and-harvest)") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(outdir, "exploration_depth_bidx.png"), p1, width = 8, height = 5, dpi = 130)

p2 <- ggplot(d, aes(b, vs_unif, color = signal)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  labs(x = "budget scale b", y = "out(S8) - out(S2)",
       title = "Does deep funding make paid information work? (discrimination value vs uniform)",
       subtitle = "At tau_k=100 all information must come from funded output; positive = the funder learned who to fund") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(outdir, "exploration_depth_vsunif.png"), p2, width = 8, height = 5, dpi = 130)

# ---- Poverty figure: schedule tilt requires learnable output ----
d2 <- sp %>% mutate(signal = factor(ifelse(tau_k < 10, "signal on", "no signal")))
p3 <- ggplot(d2, aes(r_min, alpha_S8_mean, color = signal)) +
  geom_hline(yintercept = 1/6, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) + scale_x_log10() +
  labs(x = "r_min (baseline resources, log)", y = "round-1 budget share (S8)",
       title = "Money follows information: round-1 share drops only where output is talent-revealing",
       subtitle = "eps~0, T=6, heavy tail. Dashed = even split (1/6). Poverty flattens the tilt; wealth + no signal defers spending.") +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave(file.path(outdir, "exploration_poverty_r1share.png"), p3, width = 8, height = 5, dpi = 130)

# ---- S8 full schedules (rawlong) for depth cells, no-signal ----
rl <- readRDS(file.path(indir, "exploration_depth_rawlong.rds"))
g <- sd_ %>% select(cell_id = any_of("cell_id"))  # cell ids live in raw; map via grid order
grid <- expand.grid(b = c(0.5, 1.5, 3, 6), tau_k = c(1, 100))
grid$cell_id <- seq_len(nrow(grid))
sched <- rl %>% filter(strategy == 8) %>% group_by(cell_id, round) %>%
  summarise(alpha_t = mean(alpha_t), .groups = "drop") %>% left_join(grid, by = "cell_id")
p4 <- ggplot(sched %>% filter(tau_k == 100), aes(round, alpha_t, color = factor(b), group = b)) +
  geom_hline(yintercept = 1/6, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  labs(x = "round", y = "budget share alpha_t", color = "b",
       title = "Forward planner's full schedule, no free signal (poverty, eps~0)",
       subtitle = "Even split = 1/6 dashed. Seed-and-harvest = below-even early, above-even late.") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "exploration_depth_schedules.png"), p4, width = 8, height = 5, dpi = 130)

cat("\nfigures written to", outdir, "\n")
