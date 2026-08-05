suppressPackageStartupMessages({library(parallel)})
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R")
outdir <- "sweep_results/_probe"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
t0 <- Sys.time()
grid <- expand.grid(r_min = c(0.001, 0.01, 0.1, 1),
                    epsilon = c(1e-4, 0.05, 0.1, 0.2))
seeds <- 1:12; Tr <- 5; ks <- 1.3   # heavy talent tail -> max info value
cores <- max(1L, detectCores() - 1L)
cell <- function(rm, ep) {
  rows <- mclapply(seeds, function(sd) {
    res <- run_simulation_T(seed=sd, T_rounds=Tr, n=50, epsilon=ep, b=0.5,
      M=150, tau_k=1.0, k_shape=ks, r_min=rm, budget_ref="K",
      allocator="smooth", strategies=c(5,7,8))
    s <- res$strategies
    c(b8=s[[8]]$b_idx, a1=s[[8]]$alpha_t[1],
      PG=s[[8]]$total_expected-s[[5]]$total_expected,
      sig=s[[8]]$total_expected-s[[7]]$total_expected)
  }, mc.cores=cores)
  d <- do.call(rbind, rows)
  data.frame(r_min=rm, epsilon=ep, b_idx_S8=mean(d[,"b8"]),
             r1_share=mean(d[,"a1"]), PG=mean(d[,"PG"]), signal=mean(d[,"sig"]))
}
tab <- do.call(rbind, Map(cell, grid$r_min, grid$epsilon))
tab$load <- ifelse(tab$b_idx_S8 < 0.5, "FRONT", "back")
write.csv(tab, file.path(outdir,"probe_loweps.csv"), row.names=FALSE)
cat(sprintf("\n== T=%d, k_shape=%.1f (heavy tail), tau_k=1, budget_ref=K, smooth, %d seeds ==\n", Tr, ks, length(seeds)))
cat(sprintf("even split round-1 share = 1/T = %.3f\n", 1/Tr))
for (ep in unique(tab$epsilon)) {
  cat(sprintf("\n-- epsilon = %g --\n", ep))
  print(tab[tab$epsilon==ep, c("r_min","b_idx_S8","r1_share","PG","signal","load")], digits=4, row.names=FALSE)
}
cat(sprintf("\nb_idx: <0.5 = FRONT-load (invest early); >0.5 = back-load\n"))
cat(sprintf("elapsed %.0fs\n", as.numeric(Sys.time()-t0,units="secs")))
