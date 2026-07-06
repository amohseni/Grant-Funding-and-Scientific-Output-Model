# Full manifest RE-RUN under the SMOOTH allocator (Approach A), M=400 -> T_run_smooth/
# Runs the entire 16-sweep T-round manifest with allocator="smooth" (continuous
# water-fill + projected-gradient forward planner; no greedy step delta). Writes to
# a NEW dir so sweep_results/T_run_fixed/ (greedy, superseded) is preserved for
# provenance. base_params routes allocator through modifyList -> run_simulation_T.
# The model.R default stays "greedy" (v5 T2 anchor preserved) — smooth is opt-in here.
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
cat("RE-RUN (smooth): 200 seeds, M=400,", SWEEP_CORES, "cores\n")
t0 <- Sys.time()
main_sweep_T(seeds = 1:200, base_params = list(allocator = "smooth"),
             out_dir = "sweep_results/T_run_smooth", cores = SWEEP_CORES)
cat(sprintf("TOTAL wall: %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))
