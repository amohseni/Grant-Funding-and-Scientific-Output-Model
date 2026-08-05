# Package A sweeps (SWEEP_HANDOFF_2026-08-05): D3 seed-on-signal, D4 persistent floor.
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")

write_run_info <- function(dir) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(c(
    paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
    paste("date:", format(Sys.time())),
    "seeds: 1:200",
    capture.output(sessionInfo())
  ), file.path(dir, "RUN_INFO.txt"))
}

t0 <- Sys.time()
sweep_one_T("D3_seed_signal", seeds = 1:200,
            base_params = list(allocator = "smooth", strategies = c(1,2,4,5,7,8,10,11)),
            out_dir = "sweep_results/D3_seed_signal", resume = TRUE)
write_run_info("sweep_results/D3_seed_signal")

sweep_one_T("D4_seed_persistent", seeds = 1:200,
            base_params = list(allocator = "smooth", strategies = c(1,2,4,6),
                               seed_persistent = TRUE),
            out_dir = "sweep_results/D4_seed_persistent", resume = TRUE)
write_run_info("sweep_results/D4_seed_persistent")
cat(sprintf("PACKAGE A SWEEPS DONE in %.0fs\n", as.numeric(Sys.time() - t0, units = "secs")))
