# D-4 (addendum): allocation convergence to the gap rule g* = max(cK - R, 0).
# corr(g_realized_round1, g*) as tau_k varies; output gap to the full-information oracle.
suppressPackageStartupMessages(library(parallel))
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R")
outdir <- "sweep_results/D_gap_convergence"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
grid <- expand.grid(tau_k = c(0.05, 0.3, 1, 3, 10, 20), k_shape = c(1.3, 2))
seeds <- 1:50
cores <- max(1L, detectCores() - 1L)

gap_rule <- function(K, R, budget) {
  # water-fill c: sum_i max(c*K_i - R_i, 0) = budget
  f <- function(cc) sum(pmax(cc * K - R, 0)) - budget
  cc <- uniroot(f, c(1e-9, 1e6), tol = 1e-10)$root
  pmax(cc * K - R, 0)
}

cell <- function(tk, ks) {
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks,
                           allocator = "smooth", strategies = c(1, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks,
                           allocator = "smooth", strategies = c(5), oracle = TRUE)
    K <- rp$K_at_start; R <- rp$R0_at_start   # truth from the detail=TRUE run (same seed/population)
    tranche <- rp$params$B_total / 2
    gstar <- gap_rule(K, R, tranche)
    g5 <- rp$strategies[[5]]$g_rounds[[1]]
    g4 <- rp$strategies[[4]]$g_rounds[[1]]
    c(corr_S5 = suppressWarnings(cor(g5, gstar)),
      corr_S4 = suppressWarnings(cor(g4, gstar)),
      gap_S5  = ro$strategies[[5]]$total_expected - rp$strategies[[5]]$total_expected,
      gap_S4  = ro$strategies[[5]]$total_expected - rp$strategies[[4]]$total_expected,
      out_S1  = rp$strategies[[1]]$total_expected)
  }, mc.cores = cores)
  d <- do.call(rbind, rows)
  data.frame(tau_k = tk, k_shape = ks,
             corr_S5 = mean(d[,"corr_S5"], na.rm=TRUE), corr_S5_se = sd(d[,"corr_S5"], na.rm=TRUE)/sqrt(nrow(d)),
             corr_S4 = mean(d[,"corr_S4"], na.rm=TRUE),
             gap_S5 = mean(d[,"gap_S5"]), gap_S4 = mean(d[,"gap_S4"]),
             gap_S5_pctS1 = 100*mean(d[,"gap_S5"]/d[,"out_S1"]))
}
tab <- do.call(rbind, Map(cell, grid$tau_k, grid$k_shape))
write.csv(tab, file.path(outdir, "gap_convergence.csv"), row.names = FALSE)
print(tab, digits = 3, row.names = FALSE)
writeLines(c(paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
             paste("date:", format(Sys.time())), "seeds: 1:50, oracle mode"),
           file.path(outdir, "RUN_INFO.txt"))
cat(sprintf("D-4 DONE in %.0fs\n", as.numeric(Sys.time()-t0, units="secs")))
