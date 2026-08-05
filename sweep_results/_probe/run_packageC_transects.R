# Package C2 refinement: 5-point transects across the b_idx=0.5 contour at the two
# most informative crossings (eps_paid = 0.3 and 0.85), 200 seeds (handoff rule).
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
SWEEP_CONFIGS_T$bload_transect <- list(
  name = "bload_transect", tier = 1,
  description = "Transects across the b_idx=0.5 contour of the decoupled-growth grid.",
  grid_fn = function() rbind(
    cbind(expand.grid(eps_free = c(0.0125, 0.025, 0.05, 0.1, 0.2)),  eps_paid = 0.3),
    cbind(expand.grid(eps_free = c(0.05, 0.1, 0.15, 0.2, 0.3)),      eps_paid = 0.85)) |>
    transform(T_rounds = 5, epsilon = 0.3),
  varied_params = c("eps_free", "eps_paid"),
  primary_plot = list(type = "line", x_var = "eps_free", y_var = "b_idx_S8_mean",
                      color_var = "eps_paid", title = "b_idx across the contour transects",
                      y_label = "b_idx (S8)"),
  secondary_plot = NULL)
sweep_one_T("bload_transect", seeds = 1:200,
            base_params = list(allocator = "smooth", strategies = c(1, 5, 8)),
            out_dir = "sweep_results/bload_decouple", resume = TRUE)
cat("TRANSECTS DONE\n")
