suppressPackageStartupMessages({library(parallel)})
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R")
outdir <- "sweep_results/_probe"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
t0 <- Sys.time()
grid <- expand.grid(r_min = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1),
                    tau_k = c(0.3, 10))
seeds <- 1:12; eps <- 0.6; Tr <- 6
cores <- max(1L, detectCores() - 1L)
cell <- function(rm, tk) {
  rows <- mclapply(seeds, function(sd) {
    res <- run_simulation_T(seed=sd, T_rounds=Tr, n=50, epsilon=eps, b=0.5,
      M=120, tau_k=tk, r_min=rm, budget_ref="K", allocator="smooth",
      strategies=c(5,7,8))
    s <- res$strategies
    c(b8=s[[8]]$b_idx, PG=s[[8]]$total_expected-s[[5]]$total_expected,
      sig=s[[8]]$total_expected-s[[7]]$total_expected)
  }, mc.cores=cores)
  d <- do.call(rbind, rows)
  data.frame(r_min=rm, tau_k=tk, b_idx_S8=mean(d[,"b8"]),
             PG=mean(d[,"PG"]), signal=mean(d[,"sig"]))
}
tab <- do.call(rbind, Map(cell, grid$r_min, grid$tau_k))
tab$load <- ifelse(tab$b_idx_S8 < 0.5, "FRONT", "back")
write.csv(tab, file.path(outdir,"probe_consolidated.csv"), row.names=FALSE)
cat(sprintf("\n== eps=%.2f T=%d, budget_ref=K (fixed purse), smooth, %d seeds ==\n", eps, Tr, length(seeds)))
for (tk in unique(tab$tau_k)) {
  cat(sprintf("\n-- tau_k=%.1f (%s) --\n", tk, ifelse(tk<1,"sharp signal","noisy signal")))
  print(tab[tab$tau_k==tk, c("r_min","b_idx_S8","PG","signal","load")], digits=4, row.names=FALSE)
}
cat(sprintf("\nb_idx: <0.5 = FRONT-load (invest early); >0.5 = back-load\n"))
cat(sprintf("elapsed %.0fs\n", as.numeric(Sys.time()-t0,units="secs")))
