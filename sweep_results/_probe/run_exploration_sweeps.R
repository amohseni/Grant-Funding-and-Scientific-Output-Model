# Run the three exploration sweeps (paid-information regime). Strategies include S2
# (uniform) so out_S8 - out_S2 measures whether the planner discriminates at all.
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
t0 <- Sys.time()
for (nm in c("exploration_corner", "exploration_poverty", "exploration_depth")) {
  sweep_one_T(nm, seeds = 1:64,
              base_params = list(allocator = "smooth", strategies = c(2, 5, 7, 8), M = 200),
              out_dir = "sweep_results/T_run_smooth_supplement", resume = TRUE)
}
cat(sprintf("TOTAL elapsed %.0fs\n", as.numeric(Sys.time() - t0, units = "secs")))
