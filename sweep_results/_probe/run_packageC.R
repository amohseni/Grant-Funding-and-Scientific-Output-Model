# Package C2 (SWEEP_HANDOFF_2026-08-05): back-loading attribution sweep.
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
t0 <- Sys.time()
sweep_one_T("bload_decouple", seeds = 1:200,
            base_params = list(allocator = "smooth", strategies = c(1, 5, 8)),
            out_dir = "sweep_results/bload_decouple", resume = TRUE)
writeLines(c(
  paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
  paste("date:", format(Sys.time())),
  "seeds: 1:200",
  capture.output(sessionInfo())
), "sweep_results/bload_decouple/RUN_INFO.txt")
cat(sprintf("PACKAGE C2 DONE in %.0f min\n", as.numeric(Sys.time() - t0, units = "mins")))
