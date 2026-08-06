# E-2 (addendum): headline surface at scale. (a) horizon_growth at n=200;
# (b) T x b at eps=0.85. 100 seeds each.
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
t0 <- Sys.time()
SWEEP_CONFIGS_T$E2_horizon_scale_hieps <- list(
  name = "E2_horizon_scale_hieps", tier = 1,
  description = "T x b at eps=0.85 (budget-independence of the forward advantage at high growth).",
  grid_fn = function() cbind(expand.grid(T_rounds = c(2,3,5), b = c(0.1, 0.5, 1)), epsilon = 0.85),
  varied_params = c("T_rounds", "b"),
  primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                      color_var = "b", title = "Forward advantage vs horizon by budget (eps=0.85)",
                      y_label = "S8 - S5"),
  secondary_plot = NULL)
sweep_one_T("E2_horizon_scale_hieps", seeds = 1:100,
            base_params = list(allocator = "smooth", strategies = c(1,5,8)),
            out_dir = "sweep_results/E2_headline_scale", resume = TRUE)
cat(sprintf("[E2b done] %.0fs\n", as.numeric(Sys.time()-t0, units="secs")))
sweep_one_T("horizon_growth", seeds = 1:100,
            base_params = list(allocator = "smooth", n = 200, strategies = c(1,5,8)),
            out_dir = "sweep_results/E2_headline_scale_n200", resume = TRUE)
for (d in c("sweep_results/E2_headline_scale", "sweep_results/E2_headline_scale_n200"))
  writeLines(c(paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
               paste("date:", format(Sys.time())), "seeds: 1:100"), file.path(d, "RUN_INFO.txt"))
cat(sprintf("E-2 DONE in %.0f min\n", as.numeric(Sys.time()-t0, units="mins")))
