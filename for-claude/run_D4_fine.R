# D-4 (fine grid + alpha_K = 3.5): gap-rule convergence and value of review.
# Same cell definition as sweep_results/_probe/run_D4.R (hash 21b0d9a); run in the
# cloud container session 2026-09-08 on the staged model.R, validated bit-identical
# to the canonical D_gap_convergence cell (tau_k=1, k_shape=2) before this run.
suppressPackageStartupMessages(library(parallel))
source("model.R")  # run from the repo root
gap_rule <- function(K, R, budget) {
  f <- function(cc) sum(pmax(cc * K - R, 0)) - budget
  cc <- uniroot(f, c(1e-9, 1e6), tol = 1e-10)$root
  pmax(cc * K - R, 0)
}
cell <- function(tk, ks, seeds = 1:50) {
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks,
                           allocator = "smooth", strategies = c(1, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks,
                           allocator = "smooth", strategies = c(5), oracle = TRUE)
    K <- rp$K_at_start; R <- rp$R0_at_start
    tranche <- rp$params$B_total / 2
    gstar <- gap_rule(K, R, tranche)
    g5 <- rp$strategies[[5]]$g_rounds[[1]]
    g4 <- rp$strategies[[4]]$g_rounds[[1]]
    c(corr_S5 = suppressWarnings(cor(g5, gstar)),
      corr_S4 = suppressWarnings(cor(g4, gstar)),
      gap_S5  = ro$strategies[[5]]$total_expected - rp$strategies[[5]]$total_expected,
      gap_S4  = ro$strategies[[5]]$total_expected - rp$strategies[[4]]$total_expected,
      out_S1  = rp$strategies[[1]]$total_expected)
  }, mc.cores = max(1L, detectCores() - 1L))
  d <- do.call(rbind, rows)
  data.frame(tau_k = tk, k_shape = ks,
             corr_S5 = mean(d[,"corr_S5"], na.rm=TRUE), corr_S5_se = sd(d[,"corr_S5"], na.rm=TRUE)/sqrt(nrow(d)),
             corr_S4 = mean(d[,"corr_S4"], na.rm=TRUE),
             gap_S5 = mean(d[,"gap_S5"]), gap_S4 = mean(d[,"gap_S4"]),
             gap_S5_pctS1 = 100*mean(d[,"gap_S5"]/d[,"out_S1"]),
             gap_S4_pctS1 = 100*mean(d[,"gap_S4"]/d[,"out_S1"]))
}
grid <- expand.grid(tau_k = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 14, 20),
                    k_shape = c(1.3, 2, 3.5))
t0 <- Sys.time()
tab <- do.call(rbind, Map(cell, grid$tau_k, grid$k_shape))
write.csv(tab, "gap_convergence_fine.csv", row.names = FALSE)
cat(sprintf("DONE %d cells in %.0fs\n", nrow(tab), as.numeric(Sys.time()-t0, units="secs")))
