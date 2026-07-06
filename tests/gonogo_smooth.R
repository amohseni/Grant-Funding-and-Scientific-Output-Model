# ============================================================
# DECISIVE go/no-go test for the smooth allocator (Approach A, Phase 2)
# ============================================================
# Reproduces the fwd_vs_myo_P (S7-S4) contrast that hosted the b=0.3 resonance,
# comparing THREE allocators head-to-head over a b-grid:
#   smooth      : continuous water-fill + projected-gradient (NO n_steps)
#   greedy ns50 : the buggy baseline that shows the resonance comb
#   greedy ns800: the converged fine-greedy reference (the "truth")
#
# PASS criteria (docs/SMOOTH_ALLOCATOR_HANDOFF.md sec.6):
#   A) eps=0  : smooth S7-S4 ~ 0 across ALL b (no spike at b=0.30), i.e. it lands
#      on greedy-ns800's converged ~0 row WITHOUT needing fine granularity.
#   B) eps>0  : smooth S7-S4 matches greedy-ns800 (real forward advantage
#      preserved), and does NOT reproduce the greedy-ns50 comb.

setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages(library(parallel)); source("model.R")
NC <- max(1, detectCores() - 1)
SEEDS <- 1:100
BGRID <- c(0.20, 0.24, 0.28, 0.30, 0.32, 0.36, 0.50)

# S7-S4 contrast for a (eps, b, allocator, n_steps) cell, averaged over seeds.
PGp <- function(eps, b, allocator, ns) {
  v <- unlist(mclapply(SEEDS, function(s) {
    r <- run_simulation_T(seed = s, T_rounds = 2, n = 50, epsilon = eps, b = b,
                          tau_k = 1, k_shape = 2, r_shape = 2, gamma = 1,
                          x_seed = 0.25, M = 400, n_steps = ns,
                          strategies = c(4, 7), allocator = allocator)
    r$strategies[[7]]$total_expected - r$strategies[[4]]$total_expected
  }, mc.cores = NC))
  c(mean = mean(v), se = sd(v) / sqrt(length(v)))
}

row <- function(eps, b) {
  sm <- PGp(eps, b, "smooth", 50)     # n_steps ignored by smooth
  g5 <- PGp(eps, b, "greedy", 50)
  g8 <- PGp(eps, b, "greedy", 800)
  cat(sprintf("  b=%.2f | smooth %+.4f (z%+5.1f) | greedy50 %+.4f (z%+5.1f) | greedy800 %+.4f (z%+5.1f)\n",
              b, sm[1], sm[1]/sm[2], g5[1], g5[1]/g5[2], g8[1], g8[1]/g8[2]))
  c(smooth = sm[1], greedy50 = g5[1], greedy800 = g8[1])
}

t0 <- proc.time()[3]
cat(sprintf("go/no-go: %d seeds, %d cores, b-grid={%s}\n\n",
            length(SEEDS), NC, paste(BGRID, collapse=",")))

cat("=== A) eps=0 : smooth S7-S4 must be ~0 across b (greedy50 shows the comb) ===\n")
A <- t(sapply(BGRID, function(b) row(0.0, b)))

cat("\n=== B) eps=0.1 : smooth must MATCH greedy800 (not the greedy50 comb) ===\n")
B <- t(sapply(BGRID, function(b) row(0.1, b)))

cat(sprintf("\n--- VERDICT (%.0fs) ---\n", proc.time()[3] - t0))
maxA_sm <- max(abs(A[, "smooth.mean"])); maxA_g5 <- max(abs(A[, "greedy50.mean"]))
cat(sprintf("A) eps=0  : max|smooth S7-S4| = %.4f   (greedy50 max = %.4f)\n", maxA_sm, maxA_g5))
d_sm_g8 <- max(abs(B[, "smooth.mean"] - B[, "greedy800.mean"]))
d_g5_g8 <- max(abs(B[, "greedy50.mean"] - B[, "greedy800.mean"]))
cat(sprintf("B) eps=.1 : max|smooth - greedy800| = %.4f   (greedy50 vs greedy800 = %.4f)\n",
            d_sm_g8, d_g5_g8))
passA <- maxA_sm < 0.02
passB <- d_sm_g8 < 0.02
cat(sprintf("VERDICT: A %s, B %s => %s\n",
            ifelse(passA, "PASS", "FAIL"), ifelse(passB, "PASS", "FAIL"),
            ifelse(passA && passB, "GO (integrate + re-run sweeps)", "NO-GO (use fallback)")))
