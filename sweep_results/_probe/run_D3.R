setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
sweep_one_T("D_resource_ablation", seeds = 1:200,
            base_params = list(allocator = "smooth", strategies = c(1,4,5,7,8)),
            out_dir = "sweep_results/D_resource_ablation", resume = TRUE)
writeLines(c(paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
             paste("date:", format(Sys.time())), "seeds: 1:200"),
           "sweep_results/D_resource_ablation/RUN_INFO.txt")
cat("D-3 DONE\n")
