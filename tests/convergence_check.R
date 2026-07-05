# Does high-T forward underperformance converge as granularity->fine? And is high-eps
# (the "positive" phase-diagram region) also contaminated? Report PG and the raw
# S7(fwd-pubs) & S4(myo-pubs) output levels to check for planner pathology.
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages(library(parallel))
src <- readLines("app.R", warn=FALSE); end <- grep("^shinyApp\\(", src)[1]-1L
eval(parse(text=paste(src[1:end],collapse="\n")), envir=globalenv()); source("simulate_T.R")
NC <- max(1, detectCores()-1)
stat <- function(n,Tr,eps,ns,seeds){ M<-do.call(rbind,mclapply(seeds,function(sd){
  r<-run_simulation_T(seed=sd,T_rounds=Tr,n=n,epsilon=eps,b=0.5,tau_k=1,k_shape=2,
                      x_seed=0.25,M=200,n_steps=ns,strategies=c(4,5,7,8))
  o<-sapply(c(4,5,7,8),function(s) r$strategies[[s]]$total_expected)
  c(PG=o[4]-o[2], PGp=o[3]-o[1], S8=o[4], S5=o[2])},mc.cores=NC))
  c(PG=mean(M[,"PG"]),se=sd(M[,"PG"])/sqrt(nrow(M)),PGp=mean(M[,"PGp"]),S8=mean(M[,"S8"]),S5=mean(M[,"S5"]))}
sd_ <- 1:30
cat("Convergence of PG(T=5) vs n_steps, 30 seeds:\n")
cat(sprintf("%5s %8s | %13s %13s | %8s %8s\n","eps","n_steps","PG=S8-S5","PGpubs=S7-S4","S8","S5"))
for(eps in c(0.3,0.85)) for(ns in c(50,200,800)){
  r<-stat(50,5,eps,ns,sd_)
  cat(sprintf("%5.2f %8d | %+7.3f±%.2f %+7.3f      | %8.2f %8.2f\n",eps,ns,r["PG"],r["se"],r["PGp"],r["S8"],r["S5"]))
}
cat("\nT=3 for reference (eps=0.3):\n")
for(ns in c(50,200,800)){ r<-stat(50,3,0.3,ns,sd_); cat(sprintf("  n_steps=%3d: PG=%+7.3f±%.2f\n",ns,r["PG"],r["se"])) }
