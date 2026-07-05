# Find the n_steps at which PG(S8-S5) converges for high-T cells, so horizon_long
# can be re-run at trustworthy granularity. Reports per-round increments = 2*n_steps/T.
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages(library(parallel)); source("model.R")
NC<-max(1,detectCores()-1)
PG<-function(Tr,eps,ns,seeds){v<-unlist(mclapply(seeds,function(sd){r<-run_simulation_T(seed=sd,T_rounds=Tr,n=50,epsilon=eps,b=0.5,tau_k=1,k_shape=2,x_seed=0.25,M=400,n_steps=ns,strategies=c(5,8));r$strategies[[8]]$total_expected-r$strategies[[5]]$total_expected},mc.cores=NC));mean(v)}
seeds<-1:24; NS<-c(50,100,200,400,800)
cat("PG(S8-S5) vs n_steps (per-round increments in parens), 24 seeds:\n\n")
cat(sprintf("%3s %5s |",'T','eps')); for(ns in NS) cat(sprintf(" ns=%-4d",ns)); cat("\n")
for(cell in list(c(5,0.3),c(6,0.3),c(8,0.3),c(10,0.3),c(8,0.85),c(10,0.85))){
  Tr<-cell[1];eps<-cell[2]; cat(sprintf("%3d %5.2f |",Tr,eps))
  for(ns in NS) cat(sprintf(" %+6.3f",PG(Tr,eps,ns,seeds)))
  cat(sprintf("   incr/round @ns=800: %d\n", round(2*800/Tr)))
}
