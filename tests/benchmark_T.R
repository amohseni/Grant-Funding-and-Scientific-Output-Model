# ============================================================
# tests/benchmark_T.R  —  §8 Step 1 micro-benchmark + projection
# ------------------------------------------------------------
# Measures per-run cost vs T (all 9 strategies), plus stress cells
# (large n, heavy tail), then projects full-manifest wall-clock and
# writes benchmark_report.txt with a go/no-go.
# ============================================================
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages(library(parallel))
source("model.R"); source("sweep.R"); source("sweep_T.R")

NCORES <- SWEEP_CORES
time_cell <- function(seeds, ...) {
  t0 <- Sys.time()
  for (sd in seeds) run_simulation_T(seed = sd, M = 200, strategies = 1:9, ...)
  as.numeric(Sys.time() - t0, units = "secs") / length(seeds)
}

seeds <- 1:10
cat("=== per-run cost vs T (default cell: n=50, M=200, all 9 strat) ===\n")
costT <- numeric(5)
for (Tr in 1:5) {
  costT[Tr] <- time_cell(seeds, T_rounds = Tr, n = 50, epsilon = 0.1, b = 0.5, k_shape = 2.0)
  cat(sprintf("  T=%d : %.4f s/run\n", Tr, costT[Tr]))
}
# stress cells
c_bign  <- time_cell(seeds, T_rounds = 2, n = 200, epsilon = 0.1, b = 0.5, k_shape = 2.0)
c_htail <- time_cell(seeds, T_rounds = 5, n = 50,  epsilon = 0.1, b = 0.5, k_shape = 1.3)
c_bigb  <- time_cell(seeds, T_rounds = 5, n = 50,  epsilon = 0.1, b = 1.0, k_shape = 2.0)
cat(sprintf("  stress n=200,T=2: %.4f | heavy-tail k=1.3,T=5: %.4f | b=1.0,T=5: %.4f\n",
            c_bign, c_htail, c_bigb))

# ---- project full manifest ----
# Tier-1 cells per T: horizon_growth 5 + horizon_noise 5 + horizon_scale 4 = 14 at each T in 1..5.
tier1_per_T <- 14
tier234_cells <- sum(vapply(SWEEP_CONFIGS_T, function(cf)
  if (cf$tier > 1) nrow(cf$grid_fn()) else 0, numeric(1)))
SEEDS_FULL <- 200
# serial seconds
serial <- SEEDS_FULL * (tier1_per_T * sum(costT) +      # Tier-1: 14 cells at each T=1..5
                        tier234_cells * costT[2])        # Tier 2-4 all at T=2
# pop_size & heavy-tail cells cost more than the default costT[2]; add a rough surcharge
# (pop_size n=100/200 cells and low-alpha_K cells). Use a 1.3x safety factor.
serial_safe <- serial * 1.3
parallel_hr <- serial_safe / NCORES / 3600

total_cells <- tier1_per_T * 5 + tier234_cells
lines <- c(
  "=== Full-manifest wall-clock projection (§8) ===",
  sprintf("cores available            : %d", NCORES),
  sprintf("total cells (post b-trim)  : %d", total_cells),
  sprintf("seeds per cell             : %d", SEEDS_FULL),
  sprintf("total runs                 : %d", total_cells * SEEDS_FULL),
  "",
  sprintf("per-run cost T=1..5 (s)    : %s", paste(sprintf("%.4f", costT), collapse=" ")),
  sprintf("O(T^2) check cost5/cost2   : %.2f (vs 6.25 if pure T^2)", costT[5]/costT[2]),
  sprintf("stress: n200/T2=%.3f htail/T5=%.3f b1/T5=%.3f", c_bign, c_htail, c_bigb),
  "",
  sprintf("serial estimate            : %.0f s (%.2f h)", serial, serial/3600),
  sprintf("serial + 1.3x safety       : %.0f s (%.2f h)", serial_safe, serial_safe/3600),
  sprintf("PARALLEL (%d cores)         : %.2f h", NCORES, parallel_hr),
  "",
  sprintf("GO/NO-GO (< 20h budget)    : %s", if (parallel_hr < 20) "GO — launch full manifest, no dials needed" else "NO-GO — apply §8 dials"),
  sprintf("generated (this machine)   : cost anchored to %d-seed timing", length(seeds))
)
writeLines(lines, "benchmark_report.txt")
cat("\n"); cat(lines, sep = "\n"); cat("\n\nWrote benchmark_report.txt\n")
