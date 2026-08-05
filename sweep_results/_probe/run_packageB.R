# Package B driver (SWEEP_HANDOFF_2026-08-05): Tier B concentration gate, then
# Tier A CES robustness re-runs, then Leontief boundary checks. Sequential to
# avoid core contention; each sweep checkpoints via sweep_one_T resume logic.
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")

write_run_info <- function(dir, seeds_txt) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(c(
    paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
    paste("date:", format(Sys.time())),
    paste("seeds:", seeds_txt),
    capture.output(sessionInfo())
  ), file.path(dir, "RUN_INFO.txt"))
}

t0 <- Sys.time()

# ---- Tier B: the concentration gate (72 cells, T=1, greedy-400, 200 seeds) ----
sweep_one_T("sigma_tierB", seeds = 1:200,
            base_params = list(allocator = "greedy", n_steps = 400,
                               strategies = c(1, 5)),
            out_dir = "sweep_results/sigma_tierB_concentration", resume = TRUE)
write_run_info("sweep_results/sigma_tierB_concentration", "1:200")
cat(sprintf("[tierB done] %.0fs elapsed\n", as.numeric(Sys.time() - t0, units = "secs")))

# ---- Tier B Leontief endpoint (separate; n_steps=800, degeneracy caveat) ----
# reuse the sigma_tierB grid restricted to ces_gamma = -Inf via base_params override:
SWEEP_CONFIGS_T$sigma_tierB_leontief <- modifyList(SWEEP_CONFIGS_T$sigma_tierB, list(
  name = "sigma_tierB_leontief",
  description = "Leontief endpoint of the concentration gate (ces_gamma=-Inf, n_steps=800).",
  grid_fn = function() cbind(
    expand.grid(ces_gamma = -Inf, k_shape = c(1.3, 2, 3.5), b = c(0.1, 0.5, 1.0)),
    T_rounds = 1)))
sweep_one_T("sigma_tierB_leontief", seeds = 1:200,
            base_params = list(allocator = "greedy", n_steps = 800,
                               strategies = c(1, 5)),
            out_dir = "sweep_results/sigma_tierB_concentration", resume = TRUE)
cat(sprintf("[tierB leontief done] %.0fs elapsed\n", as.numeric(Sys.time() - t0, units = "secs")))

# ---- Tier A: four headline sweeps x ces_gamma in {0, -3}, greedy-400, 100 seeds ----
tierA_sweeps <- c("signal_value", "horizon_growth", "seed_value", "correlation")
for (gc in c(0, -3)) {
  dir <- sprintf("sweep_results/sigma_tierA_gc%s", ifelse(gc == 0, "0", "m3"))
  for (nm in tierA_sweeps) {
    sweep_one_T(nm, seeds = 1:100,
                base_params = list(allocator = "greedy", n_steps = 400, ces_gamma = gc),
                out_dir = dir, resume = TRUE)
    cat(sprintf("[tierA %s gc=%g done] %.0fs elapsed\n", nm, gc,
                as.numeric(Sys.time() - t0, units = "secs")))
  }
  write_run_info(dir, "1:100")
}

# ---- Tier A Leontief boundary check: signal_value only, n_steps=800, 50 seeds ----
sweep_one_T("signal_value", seeds = 1:50,
            base_params = list(allocator = "greedy", n_steps = 800, ces_gamma = -Inf),
            out_dir = "sweep_results/sigma_tierA_leontief", resume = TRUE)
write_run_info("sweep_results/sigma_tierA_leontief", "1:50")

cat(sprintf("PACKAGE B DONE in %.0f min\n", as.numeric(Sys.time() - t0, units = "mins")))
