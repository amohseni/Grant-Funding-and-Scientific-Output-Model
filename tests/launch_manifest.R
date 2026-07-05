# Full T-round manifest launch (§5), 200 seeds, checkpointed to sweep_results/T_run/.
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
src <- readLines("app.R", warn=FALSE); end <- grep("^shinyApp\\(", src)[1]-1L
eval(parse(text=paste(src[1:end],collapse="\n")), envir=globalenv())
source("simulate_T.R"); source("sweep.R"); source("sweep_T.R")
cat("Launch: 200 seeds,", SWEEP_CORES, "cores, ", length(SWEEP_CONFIGS_T), "sweeps\n")
t0 <- Sys.time()
main_sweep_T(seeds = 1:200, out_dir = "sweep_results/T_run", cores = SWEEP_CORES)
cat(sprintf("TOTAL wall: %.1f min\n", as.numeric(Sys.time()-t0, units="mins")))
