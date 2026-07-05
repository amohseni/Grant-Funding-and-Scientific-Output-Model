# reproduce.R
# ============================================================
# One-command reproduction of the entire study, from seeds.
# Run from the project root:  Rscript reproduce.R
#
# Regenerates the full parameter-sweep dataset and re-runs the
# validation suite. Everything is deterministic given the seeds
# (each run does set.seed() internally), so a fresh checkout
# reproduces the published results bit-for-bit.
#
# Runtime: ~45 min on ~8 cores (the 16-sweep manifest ~30 min +
# validation). Requires: R >= 4.4, and the packages listed in
# docs/ENVIRONMENT.md (shiny not needed for reproduction).
# ============================================================

# Run this from the PROJECT ROOT (the directory containing model.R).
if (!file.exists("model.R"))
  stop("Run reproduce.R from the project root (the folder containing model.R).")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(scales); library(parallel)
})

cat("== Environment ==\n")
cat(R.version.string, "|", R.version$platform, "| cores:", detectCores(), "\n\n")

source("model.R")     # the model engine (base R)
source("sweep.R")     # sweep infrastructure
source("sweep_T.R")   # the 16-sweep manifest (SWEEP_CONFIGS_T)

# ---- Stage 1: validation gate (fast) --------------------------------------
cat("== Stage 1: T=2 regression anchor (must be bit-identical to v5) ==\n")
system("Rscript tests/test_T2_reduction.R")

# ---- Stage 2: full manifest (main 13 sweeps + 3 supplements = 16) ---------
cat("\n== Stage 2: full parameter-sweep manifest (200 trials/cell) ==\n")
t0 <- Sys.time()
main_sweep_T(seeds = 1:200, out_dir = "sweep_results/T_run_fixed")   # resume=TRUE: skips any already present
cat(sprintf("   manifest wall-clock: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))

# ---- Stage 3: validity assertions -----------------------------------------
cat("\n== Stage 3: §6 validity assertions ==\n")
system("Rscript tests/assertions_T.R")

cat("\nDone. Data in sweep_results/T_run_fixed/ ; see T_round_extension/DATA_DICTIONARY.md.\n")
cat("Assertions log: docs/assertions_log.txt (or repo root if regenerated).\n")
