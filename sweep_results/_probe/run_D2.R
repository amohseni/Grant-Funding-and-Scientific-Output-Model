setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
SWEEP_CONFIGS_T$D_misspecified_trust <- list(
  name = "D_misspecified_trust", tier = 1,
  description = "Misspecified trust: true grant-signal noise x believed noise x talent tail, T=2 smooth.",
  grid_fn = function() expand.grid(tau_k_true = c(0.3, 1, 3),
                                   tau_k_belief = c(0.1, 0.3, 1, 3, 10),
                                   k_shape = c(1.3, 2)),
  varied_params = c("tau_k_true", "tau_k_belief", "k_shape"),
  primary_plot = list(type = "line", x_var = "tau_k_belief", y_var = "signal_myo_mean",
                      color_var = "tau_k_true", title = "Signal value under misweighted trust (S5-S4)",
                      y_label = "S5 - S4"),
  secondary_plot = list(type = "line", x_var = "tau_k_belief", y_var = "signal_fwd_mean",
                        color_var = "tau_k_true", title = "S8 - S7 under misweighted trust",
                        y_label = "S8 - S7"))
sweep_one_T("D_misspecified_trust", seeds = 1:200,
            base_params = list(allocator = "smooth", strategies = c(1,2,4,5,7,8)),
            out_dir = "sweep_results/D_misspecified_trust", resume = TRUE)
writeLines(c(paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
             paste("date:", format(Sys.time())), "seeds: 1:200"),
           "sweep_results/D_misspecified_trust/RUN_INFO.txt")
cat("D-2 DONE\n")
