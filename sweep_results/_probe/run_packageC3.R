# Package C3 (optional corroboration): resource_regime at 200 seeds -> sweep_results/resource_regime/
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
sweep_one_T("resource_regime", seeds = 1:200,
            base_params = list(allocator = "smooth", strategies = c(1, 5, 7, 8), M = 200),
            out_dir = "sweep_results/resource_regime", resume = TRUE)
writeLines(c(paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
             paste("date:", format(Sys.time())), "seeds: 1:200"),
           "sweep_results/resource_regime/RUN_INFO.txt")
cat("C3 DONE\n")
