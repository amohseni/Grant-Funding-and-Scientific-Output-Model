# Probe the pure-exploration corner: no free signal (tau_k huge), poverty (r_min~0),
# heavy talent tail. Question: does the forward planner front-load money (b_idx<0.5),
# or seed-and-harvest (positive round-1 spend, mass late)? Full alpha_t schedules printed.
suppressPackageStartupMessages({library(parallel)})
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R")
outdir <- "sweep_results/_probe"
t0 <- Sys.time()
cells <- list(
  list(lab="CORNER: no signal, poor, eps~0",  tau_k=100, epsilon=1e-4, r_min=0.001),
  list(lab="signal on, poor, eps~0",          tau_k=1,   epsilon=1e-4, r_min=0.001),
  list(lab="no signal, poor, eps=0.3",        tau_k=100, epsilon=0.3,  r_min=0.001),
  list(lab="no signal, RICH, eps~0",          tau_k=100, epsilon=1e-4, r_min=1)
)
seeds <- 1:12; Tr <- 6; ks <- 1.3
cores <- max(1L, detectCores() - 1L)
res_rows <- list(); sched_rows <- list()
for (ci in seq_along(cells)) {
  cl <- cells[[ci]]
  rows <- mclapply(seeds, function(sd) {
    res <- run_simulation_T(seed=sd, T_rounds=Tr, n=50, epsilon=cl$epsilon, b=0.5,
      M=200, tau_k=cl$tau_k, k_shape=ks, r_min=cl$r_min, budget_ref="K",
      allocator="smooth", strategies=c(2,5,7,8))
    s <- res$strategies
    list(m=c(b8=s[[8]]$b_idx, a1=s[[8]]$alpha_t[1],
             PG=s[[8]]$total_expected-s[[5]]$total_expected,
             sig=s[[8]]$total_expected-s[[7]]$total_expected,
             vs_unif=s[[8]]$total_expected-s[[2]]$total_expected,
             out8=s[[8]]$total_expected),
         sched=s[[8]]$alpha_t)
  }, mc.cores=cores)
  m <- do.call(rbind, lapply(rows, `[[`, "m"))
  sch <- colMeans(do.call(rbind, lapply(rows, `[[`, "sched")))
  res_rows[[ci]] <- data.frame(cell=cl$lab, t(colMeans(m)))
  sched_rows[[ci]] <- data.frame(cell=cl$lab, round=seq_len(Tr), alpha_t=sch)
}
tab <- do.call(rbind, res_rows)
write.csv(tab, file.path(outdir,"probe_exploration.csv"), row.names=FALSE)
write.csv(do.call(rbind, sched_rows), file.path(outdir,"probe_exploration_sched.csv"), row.names=FALSE)
cat(sprintf("\n== T=%d, k_shape=%.1f, budget_ref=K, smooth, %d seeds. even share=%.3f ==\n\n",
            Tr, ks, length(seeds), 1/Tr))
print(tab, digits=4, row.names=FALSE)
cat("\n== mean S8 budget-share schedule alpha_t by round ==\n")
for (ci in seq_along(cells)) {
  s <- sched_rows[[ci]]
  cat(sprintf("%-34s : %s\n", cells[[ci]]$lab, paste(sprintf("%.3f", s$alpha_t), collapse=" ")))
}
cat(sprintf("\nelapsed %.0fs\n", as.numeric(Sys.time()-t0,units="secs")))
