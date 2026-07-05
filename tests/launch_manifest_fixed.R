# Full manifest RE-RUN on fixed app.R (deterministic ce_reweight), M=400 -> T_run_fixed/
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
cat("RE-RUN (fixed): 200 seeds, M=400,", SWEEP_CORES, "cores\n")
t0 <- Sys.time()
main_sweep_T(seeds = 1:200, out_dir = "sweep_results/T_run_fixed", cores = SWEEP_CORES)
cat(sprintf("TOTAL wall: %.1f min\n", as.numeric(Sys.time()-t0, units="mins")))
