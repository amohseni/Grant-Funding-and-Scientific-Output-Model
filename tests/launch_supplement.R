# Supplementary sweeps (Roadmap Stage 1.1): tail_map, info_value, horizon_long
# 200 seeds, M=400 (matches main manifest) -> sweep_results/T_run_supplement/
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2)})
source("model.R"); source("sweep.R"); source("sweep_T.R")
t0 <- Sys.time()
for (nm in c("tail_map","info_value","horizon_long"))
  sweep_one_T(nm, seeds = 1:200, out_dir = "sweep_results/T_run_supplement", cores = SWEEP_CORES)
cat(sprintf("SUPPLEMENT DONE — %.1f min\n", as.numeric(Sys.time()-t0, units="mins")))
