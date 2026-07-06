# ============================================================
# Supplementary validation of the smooth allocator, answering:
#   "Is smooth fully valid + free of the greedy artifact?" beyond the T=2 / S7-S4
#   go/no-go. Two gaps the go/no-go left open:
#     CHECK 1 — the ORIGINAL diligence spike was the GRANT contrast S8-S5 (PG),
#               not S7-S4 (P). Does smooth kill the PG comb too?
#     CHECK 2 — a SEPARATE greedy artifact lives at HIGHER T (horizon_long). Does
#               the projected-gradient planner stay well-behaved at H=5 (dominate
#               fine greedy, satisfy KKT, no stall / non-concavity failure)?
# ============================================================
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages(library(parallel)); source("model.R")
NC <- max(1, detectCores() - 1)

# ---- CHECK 1: grant-signal contrast S8-S5 (PG), the original z24 anomaly ----
cat("=== CHECK 1: PG contrast S8-S5, eps=0.1, resonance b-grid ===\n")
cat("  greedy50 hosted the original comb; smooth must be flat & match greedy800.\n")
SEEDS <- 1:80
PGg <- function(b, allocator, ns) {
  v <- unlist(mclapply(SEEDS, function(s) {
    r <- run_simulation_T(seed = s, T_rounds = 2, n = 50, epsilon = 0.1, b = b,
                          tau_k = 1, k_shape = 2, r_shape = 2, gamma = 1,
                          x_seed = 0.25, M = 400, n_steps = ns,
                          strategies = c(5, 8), allocator = allocator)
    r$strategies[[8]]$total_expected - r$strategies[[5]]$total_expected
  }, mc.cores = NC))
  c(mean(v), sd(v) / sqrt(length(v)))
}
maxsm <- 0; maxdiff <- 0
for (b in c(0.24, 0.28, 0.30, 0.32, 0.36)) {
  sm <- PGg(b, "smooth", 50); g5 <- PGg(b, "greedy", 50); g8 <- PGg(b, "greedy", 800)
  maxsm <- max(maxsm, abs(sm[1])); maxdiff <- max(maxdiff, abs(sm[1] - g8[1]))
  cat(sprintf("  b=%.2f | smooth %+.4f (z%+5.1f) | greedy50 %+.4f (z%+5.1f) | greedy800 %+.4f\n",
              b, sm[1], sm[1]/sm[2], g5[1], g5[1]/g5[2], g8[1]))
}
cat(sprintf("  -> max|smooth PG| = %.4f ; max|smooth - greedy800| = %.4f\n\n", maxsm, maxdiff))

# ---- CHECK 2: projected-gradient robustness at H=5 (single-plan diagnostics) ----
cat("=== CHECK 2: forward planner at H=5 (dominance + KKT, no stall) ===\n")
plan_obj <- function(posts, plan, gamma, epsilon) {
  v <- 0
  for (i in seq_along(posts)) {
    bar1 <- ce_reweight_posterior(posts[[i]], plan[i, 1], gamma)
    v <- v + fwd_researcher_value(posts[[i]], bar1, numeric(0), plan[i, ], gamma, epsilon)
  }
  v
}
kkt_spread <- function(posts, plan, gamma, epsilon, dg = 1e-6, thresh = 1e-6) {
  m <- c(); H <- ncol(plan)
  for (i in seq_along(posts)) {
    bar1 <- ce_reweight_posterior(posts[[i]], plan[i, 1], gamma)
    for (s in seq_len(H)) if (plan[i, s] > thresh)
      m <- c(m, fwd_marginal(posts[[i]], bar1, numeric(0), plan[i, ], s, dg, gamma, epsilon))
  }
  if (length(m) < 2) 0 else diff(range(m))
}
for (eps in c(0.0, 0.1)) for (b in c(0.3, 0.5)) {
  set.seed(11); n <- 40; H <- 5
  E_R <- 2; B_rem <- 2 * b * n * E_R
  pop <- draw_initial_population(n, 1, 2, 1, 2, 0)
  p0  <- rpois(n, pmax(lambda_rate(pop$K0, pop$R0, 1), 1e-12))
  sig <- draw_signals(pop$K0, pop$R0, 1, 1)
  posts <- build_posteriors_hist(p0, sig$sigma_r, sig$sigma_k, 400, 1, 2, 1, 2, 1, 1, 1, TRUE, FALSE)
  t <- system.time(psm <- plan_forward_smooth(posts, B_rem, 1, eps, H))[3]
  pgr <- plan_forward_ce(posts, B_rem, B_rem / H / 400, 1, eps, H)   # fine greedy
  o_s <- plan_obj(posts, psm, 1, eps); o_g <- plan_obj(posts, pgr, 1, eps)
  cat(sprintf("  H=5 eps=%.1f b=%.2f | smooth obj=%.5f spend=%.4f (t=%.2fs) | smooth-greedy400=%+.2e | KKT=%.2e\n",
              eps, b, o_s, sum(psm), t, o_s - o_g, kkt_spread(posts, psm, 1, eps)))
}
cat("\nDONE\n")
