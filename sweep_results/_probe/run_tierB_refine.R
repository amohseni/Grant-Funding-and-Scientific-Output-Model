# Tier B refinement rule: bracket within-cell reversals of gini vs ces_gamma
# with two extra points at 400 seeds (handoff B3). Reversal detected at
# (k_shape=3.5, b=1.0) [interior max near -6]; bracket with -4 and -9.
# Also run (k_shape=3.5, b=0.5) as the adjacent borderline cell.
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
SWEEP_CONFIGS_T$sigma_tierB_refine <- list(
  name = "sigma_tierB_refine", tier = 1,
  description = "Refinement: bracket the interior gini maximum (ces_gamma -4, -9) at k_shape=3.5, b in {0.5, 1.0}, 400 seeds.",
  grid_fn = function() cbind(
    expand.grid(ces_gamma = c(-4, -9), k_shape = 3.5, b = c(0.5, 1.0)),
    T_rounds = 1),
  varied_params = c("ces_gamma", "k_shape", "b"),
  primary_plot = list(type = "line", x_var = "ces_gamma", y_var = "gini_g1_S5_mean",
                      color_var = "b", title = "Refined gini near the interior max",
                      y_label = "Gini (S5)"),
  secondary_plot = NULL)
sweep_one_T("sigma_tierB_refine", seeds = 1:400,
            base_params = list(allocator = "greedy", n_steps = 400, strategies = c(1, 5)),
            out_dir = "sweep_results/sigma_tierB_concentration", resume = TRUE)
cat("TIER B REFINEMENT DONE\n")
