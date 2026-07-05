# Does the erratic high-T PG(n) dependence vanish as greedy granularity (n_steps) rises?
# If PG(T=5) across n=40/50/60 converges as n_steps grows, the high-T decline is a
# discretization artifact of the forward planner, not a structural model property.
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages(library(parallel))
source("model.R")
NC <- max(1, detectCores()-1)
PG <- function(n,Tr,eps,ns,seeds){ v<-unlist(mclapply(seeds,function(sd){
  r<-run_simulation_T(seed=sd,T_rounds=Tr,n=n,epsilon=eps,b=0.5,tau_k=1,k_shape=2,
                      x_seed=0.25,M=200,n_steps=ns,strategies=c(5,8))
  r$strategies[[8]]$total_expected-r$strategies[[5]]$total_expected},mc.cores=NC));c(mean(v),sd(v)/sqrt(length(v)))}
cat("PG(T=5, eps=0.3) vs n_steps (granularity), 80 seeds. spread across n = artifact size:\n\n")
cat(sprintf("%8s | %11s %11s %11s | %s\n","n_steps","n=40","n=50","n=60","spread(max-min)"))
for(ns in c(50,100,200,400)){
  r40<-PG(40,5,0.3,ns,1:80); r50<-PG(50,5,0.3,ns,1:80); r60<-PG(60,5,0.3,ns,1:80)
  sp<-max(r40[1],r50[1],r60[1])-min(r40[1],r50[1],r60[1])
  cat(sprintf("%8d | %+6.3f±%.2f %+6.3f±%.2f %+6.3f±%.2f | %.3f\n",
      ns,r40[1],r40[2],r50[1],r50[2],r60[1],r60[2],sp))
}
cat("\nAlso T=2 (bit-identical to v5, should be granularity-stable) as a control:\n")
for(ns in c(50,200)){ r<-PG(50,2,0.3,ns,1:80); cat(sprintf("  n=50 T=2 n_steps=%3d: PG=%+.3f±%.3f\n",ns,r[1],r[2])) }
