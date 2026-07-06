# Diagnose the b=0.3 fwd_vs_myo_P spike: greedy-discretization resonance?
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages(library(parallel)); source("model.R")
NC <- max(1, detectCores()-1); SEEDS <- 1:150; E_R <- 2
runcell <- function(seed, Tr, eps, b, ns, strat) run_simulation_T(seed=seed, T_rounds=Tr, n=50,
  epsilon=eps, b=b, tau_k=1, k_shape=2, r_shape=2, gamma=1, x_seed=0.25, M=400, n_steps=ns, strategies=strat, detail=("d" %in% strat))
PGp <- function(eps,b,ns,seeds){ v<-unlist(mclapply(seeds,function(s){r<-runcell(s,2,eps,b,ns,c(4,7));r$strategies[[7]]$total_expected-r$strategies[[4]]$total_expected},mc.cores=NC)); c(mean(v),sd(v)/sqrt(length(v))) }
PGg <- function(eps,b,ns,seeds){ v<-unlist(mclapply(seeds,function(s){r<-runcell(s,2,eps,b,ns,c(5,8));r$strategies[[8]]$total_expected-r$strategies[[5]]$total_expected},mc.cores=NC)); c(mean(v),sd(v)/sqrt(length(v))) }

cat("=== TEST 1: eps=0, T=2, b=0.3, pubs-only. fwd_vs_myo_P (S7-S4) vs n_steps ===\n")
cat("  (predict: ~0.12 at ns=50, monotone -> 0)\n")
for(ns in c(50,100,200,400,800,1600)){ r<-PGp(0,0.3,ns,SEEDS); cat(sprintf("  n_steps=%-4d  S7-S4=%+.4f (se %.4f, z%+.1f)\n",ns,r[1],r[2],r[1]/r[2])) }

cat("\n=== TEST 2: eps=0.1, T=2, fine b-grid, S7-S4 [and S8-S5] at n_steps 50 vs 800 ===\n")
cat(sprintf("  %6s | %-22s | %-22s\n","b","n_steps=50 (P; PG)","n_steps=800 (P; PG)"))
for(b in c(0.24,0.26,0.28,0.30,0.32,0.34,0.36)){
  p50<-PGp(0.1,b,50,SEEDS); g50<-PGg(0.1,b,50,SEEDS); p800<-PGp(0.1,b,800,SEEDS); g800<-PGg(0.1,b,800,SEEDS)
  cat(sprintf("  %6.2f | %+.3f ; %+.3f        | %+.3f ; %+.3f\n",b,p50[1],g50[1],p800[1],g800[1]))
}

cat("\n=== TEST 3: eps=0, T=2, n_steps=50. L1 allocation distance ||g1_S7 - g1_S4||_1 / B ===\n")
cat("  (predict: spikes at b=0.3, ~0 elsewhere)\n")
for(b in c(0.1,0.2,0.3,0.5,1.0)){
  B <- b*50*E_R
  d <- unlist(mclapply(SEEDS,function(s){ r<-runcell(s,2,0,b,50,1:9); sum(abs(r$strategies[[7]]$g_rounds[[1]]-r$strategies[[4]]$g_rounds[[1]]))/B },mc.cores=NC))
  cat(sprintf("  b=%-4g  L1/B = %.4f (se %.4f)\n",b,mean(d),sd(d)/sqrt(length(d))))
}
cat("\nDONE\n")
