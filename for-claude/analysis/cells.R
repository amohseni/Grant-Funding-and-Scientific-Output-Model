suppressPackageStartupMessages(library(parallel))
source("model.R")
cell <- function(tk, ks, seeds = 1:50) {
  rows <- mclapply(seeds, function(sd) {
    rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks,
                           allocator = "smooth", strategies = c(1, 4, 5), detail = TRUE)
    ro <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks,
                           allocator = "smooth", strategies = c(5), oracle = TRUE)
    c(S1 = rp$strategies[[1]]$total_expected, S4 = rp$strategies[[4]]$total_expected,
      S5 = rp$strategies[[5]]$total_expected, OR = ro$strategies[[5]]$total_expected)
  }, mc.cores = 2)
  d <- do.call(rbind, rows)
  data.frame(tau_k = tk, k_shape = ks,
    gap_S4_pctS1 = 100*mean((d[,"OR"]-d[,"S4"])/d[,"S1"]),
    gap_S4_pctgain = 100*mean((d[,"OR"]-d[,"S4"])/(d[,"OR"]-d[,"S1"])),
    gap_S5_pctS1 = 100*mean((d[,"OR"]-d[,"S5"])/d[,"S1"]),
    gap_S5_pctgain = 100*mean((d[,"OR"]-d[,"S5"])/(d[,"OR"]-d[,"S1"])),
    rev_pctS1 = 100*mean((d[,"S5"]-d[,"S4"])/d[,"S1"]),
    rev_pctgain5 = 100*mean((d[,"S5"]-d[,"S4"])/(d[,"S5"]-d[,"S1"])),
    oracle_gain_pctS1 = 100*mean((d[,"OR"]-d[,"S1"])/d[,"S1"]))
}
grid <- expand.grid(tau_k = c(0.05, 1, 3, 20), k_shape = c(1.3, 2, 3.5))
t0 <- Sys.time()
tab <- do.call(rbind, Map(cell, grid$tau_k, grid$k_shape))
print(tab, digits=3)
write.csv(tab, "review_cells.csv", row.names = FALSE)
cat(sprintf("DONE in %.0fs\n", as.numeric(Sys.time()-t0, units="secs")))
