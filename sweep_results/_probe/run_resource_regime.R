# Run the resource_regime sweep (front-load reversal map) — smooth allocator, decoupled purse.
# Focal strategies only (5,7,8) for speed; b_idx_S8, fwd_vs_myo_PG (S8-S5), signal_fwd (S8-S7).
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
t0 <- Sys.time()
sweep_one_T("resource_regime",
            seeds = 1:64,
            base_params = list(allocator = "smooth", strategies = c(5, 7, 8), M = 200),
            out_dir = "sweep_results/T_run_smooth_supplement",
            resume = FALSE)
cat(sprintf("TOTAL elapsed %.0fs\n", as.numeric(Sys.time() - t0, units = "secs")))
