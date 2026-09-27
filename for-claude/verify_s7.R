# S7 (timing) claims re-derivation, mirroring verify_all_claims.R sections 1, 2, 3(A5), 4(B1 scope), 5.
U <- "sweep_results"
ck <- function(lab, val, fmt = "%.4g") cat(sprintf("  %-58s %s\n", lab, sprintf(fmt, val)))

cat("== Result 1: front/back boundary in eps (resource_regime) ==\n")
s <- readRDS(file.path(U, "resource_regime/resource_regime_summary.rds"))
fl <- s[s$b_idx_S8_mean < 0.5, ]
ck("front-load cells (b_idx<0.5)", nrow(fl), "%d")
ck("max eps among front-load cells", max(fl$epsilon))
d <- s[s$epsilon == 0.85, ]; d <- d[order(d$r_min), ]
cat("  b_idx at eps=0.85, r_min 0.001..3:", sprintf("%.3f", d$b_idx_S8_mean), "\n")
cat("  round-1 share same rows:", sprintf("%.3f", d$alpha_S8_mean), "\n")

cat("== Results 2-3: depth (exploration_200) ==\n")
dd <- readRDS(file.path(U, "exploration_200/exploration_depth_summary.rds"))
dd <- dd[dd$tau_k == 100, ]; dd <- dd[order(dd$b), ]
cat("  b:", dd$b, "| S8-S2:", sprintf("%.1f", dd$out_S8_mean - dd$out_S2_mean),
    "| r1 share:", sprintf("%.3f", dd$alpha_S8_mean), "| b_idx:", sprintf("%.3f", dd$b_idx_S8_mean), "\n")
co <- readRDS(file.path(U, "exploration_200/exploration_corner_summary.rds"))
cat("  corner (b=0.5): S8-S2 by tau_k:", paste(sprintf("tau=%g: %.1f", co$tau_k, co$out_S8_mean - co$out_S2_mean), collapse=", "), "\n")
allb <- do.call(rbind, lapply(c("exploration_corner","exploration_poverty","exploration_depth"),
  function(nm) readRDS(file.path(U, sprintf("exploration_200/%s_summary.rds", nm)))[, c("b_idx_S8_mean","b_idx_S8_se")]))
ck("min b_idx_S8 across exploration cells", min(allb$b_idx_S8_mean))
ck("n exploration cells", nrow(allb), "%d")
pv <- readRDS(file.path(U, "exploration_200/exploration_poverty_summary.rds"))
cat("  poverty r1 share by r_min:", paste(sprintf("r_min=%g: %.3f", pv$r_min, pv$alpha_S8_mean), collapse=", "), "\n")

cat("== Result 4: attribution (bload_decouple) ==\n")
bd <- readRDS(file.path(U, "bload_decouple/bload_decouple_summary.rds"))
pf <- bd[bd$eps_paid==0 & bd$eps_free>0, ]; pf <- pf[order(pf$eps_free), ]
cat("  pure-free b_idx (ef .1/.3/.85):", sprintf("%.3f", pf$b_idx_S8_mean), "\n")
pp <- bd[bd$eps_free==0 & bd$eps_paid>0, ]; pp <- pp[order(pp$eps_paid), ]
cat("  pure-paid b_idx (ep .1/.3/.85):", sprintf("%.4f", pp$b_idx_S8_mean), "| max SE:", sprintf("%.4f", max(pp$b_idx_S8_se)), "\n")
ck("PG at pure-paid eps_paid=.85", bd$fwd_vs_myo_PG_mean[bd$eps_free==0 & bd$eps_paid==0.85], "%.3f")
tr <- readRDS(file.path(U, "bload_decouple/bload_transect_summary.rds"))
for (ep in c(0.3, 0.85)) {
  t1 <- tr[tr$eps_paid==ep, ]; t1 <- t1[order(t1$eps_free), ]
  lo <- max(t1$eps_free[t1$b_idx_S8_mean < 0.5]); hi <- min(t1$eps_free[t1$b_idx_S8_mean >= 0.5])
  cat(sprintf("  contour @eps_paid=%.2f: eps_free in (%.3g, %.3g), ratio %.2f-%.2f\n", ep, lo, hi, lo/ep, hi/ep))
}

cat("== Result 5 + calibration: PG magnitudes, T exponent ==\n")
hg <- readRDS(file.path(U, "T_run_smooth/horizon_growth_summary.rds"))
hl <- readRDS(file.path(U, "T_run_smooth/horizon_long_summary.rds"))
ddx <- rbind(hg[,c("T_rounds","epsilon","fwd_vs_myo_PG_mean","signal_fwd_mean","out_S1_mean")],
             hl[,c("T_rounds","epsilon","fwd_vs_myo_PG_mean","signal_fwd_mean","out_S1_mean")])
for (ep in c(0.3, 0.85)) {
  d5 <- ddx[ddx$epsilon==ep & ddx$fwd_vs_myo_PG_mean>0 & ddx$T_rounds>=2 & ddx$T_rounds<=5,]
  f <- lm(log(fwd_vs_myo_PG_mean) ~ log(T_rounds), d5)
  dl <- ddx[ddx$epsilon==ep & ddx$fwd_vs_myo_PG_mean>0 & ddx$T_rounds>=5,]
  f2 <- lm(log(fwd_vs_myo_PG_mean) ~ log(T_rounds), dl)
  cat(sprintf("  exponent eps=%.2f: T<=5: %.2f ; T=5..10: %.2f\n", ep, coef(f)[2], coef(f2)[2]))
}
t5 <- hg[hg$T_rounds==5, ]
cat("  T=5 PG vs signal value by eps:", paste(sprintf("eps=%g: PG=%.2f sig=%.2f (ratio %.3f)",
    t5$epsilon, t5$fwd_vs_myo_PG_mean, t5$signal_fwd_mean, t5$fwd_vs_myo_PG_mean/t5$signal_fwd_mean), collapse="; "), "\n")

cat("== Cobb-Douglas scope (sigma_tierA) ==\n")
g0 <- readRDS(file.path(U, "sigma_tierA_gc0/horizon_growth_summary.rds"))
ck("CD: max PG(T=5) over eps", max(g0$fwd_vs_myo_PG_mean[g0$T_rounds==5]))
ck("CD: max b_idx (T=5)", max(g0$b_idx_S8_mean[g0$T_rounds==5]), "%.3f")
g3 <- readRDS(file.path(U, "sigma_tierA_gcm3/horizon_growth_summary.rds"))
ck("gc=-3: PG(T=5, eps=.85)", g3$fwd_vs_myo_PG_mean[g3$T_rounds==5 & g3$epsilon==0.85], "%.2f")
ck("gc=-3: b_idx(T=5, eps=.85)", g3$b_idx_S8_mean[g3$T_rounds==5 & g3$epsilon==0.85], "%.3f")

cat("== Bootstrap verifications ==\n")
h0 <- read.csv(file.path(U, "bootstrap_verify/honest_schedules_eps0.csv"))
h3b <- read.csv(file.path(U, "bootstrap_verify/honest_schedules_eps03.csv"))
ev0 <- h0$out[h0$schedule=="even"]
xs <- suppressWarnings(as.numeric(sub("x([0-9.]+)_.*","\\1", h0$schedule)))
ks <- suppressWarnings(as.numeric(sub(".*_k([0-9]+)","\\1", h0$schedule)))
em <- which(!is.na(xs) & xs > ks/6 + 1e-9)
ck("eps~0: even output", ev0, "%.1f")
ck("eps~0: best early-mass minus even", max(h0$out[em]) - ev0, "%.1f")
ck("eps=0.3: best early-mass minus even", max(h3b$out[em]) - h3b$out[h3b$schedule=="even"], "%.1f")
cc <- read.csv(file.path(U, "bootstrap_verify/ce_self_consistency.csv"))
ck("CE(S8)-CE(S5) mean", mean(cc$ce_S8 - cc$ce_S5), "%.2f")
ck("seeds with CE(S8)>=CE(S5)", sum(cc$ce_S8 >= cc$ce_S5), "%d")
U <- "sweep_results"
# (a) PG / signal ratio across ALL tested horizon cells (T>=2), and PG as % of S1
hg <- readRDS(file.path(U,"T_run_smooth/horizon_growth_summary.rds"))
hl <- readRDS(file.path(U,"T_run_smooth/horizon_long_summary.rds"))
cols <- c("T_rounds","epsilon","fwd_vs_myo_PG_mean","signal_fwd_mean","out_S1_mean")
d <- rbind(hg[,cols], hl[,cols]); d <- d[d$T_rounds>=2,]
d$ratio <- d$fwd_vs_myo_PG_mean / d$signal_fwd_mean
d$pg_pct <- 100*d$fwd_vs_myo_PG_mean/d$out_S1_mean
d$sig_pct <- 100*d$signal_fwd_mean/d$out_S1_mean
cat(sprintf("cells: %d; ratio range: %.4f .. %.4f (max at T=%d eps=%g)\n",
    nrow(d), min(d$ratio), max(d$ratio), d$T_rounds[which.max(d$ratio)], d$epsilon[which.max(d$ratio)]))
cat(sprintf("ratio at default eps=0.3: %s\n", paste(sprintf("T=%d: %.3f", d$T_rounds[d$epsilon==0.3], d$ratio[d$epsilon==0.3]), collapse=" ")))
cat(sprintf("PG %% of S1 range: %.3f .. %.2f; signal %% of S1 range: %.1f .. %.1f\n",
    min(d$pg_pct), max(d$pg_pct), min(d$sig_pct), max(d$sig_pct)))
# (b) T3 percents
h0 <- read.csv(file.path(U,"bootstrap_verify/honest_schedules_eps0.csv"))
h3 <- read.csv(file.path(U,"bootstrap_verify/honest_schedules_eps03.csv"))
xs <- suppressWarnings(as.numeric(sub("x([0-9.]+)_.*","\\1", h0$schedule)))
ks <- suppressWarnings(as.numeric(sub(".*_k([0-9]+)","\\1", h0$schedule)))
em <- which(!is.na(xs) & xs > ks/6 + 1e-9)
ev0 <- h0$out[h0$schedule=="even"]; ev3 <- h3$out[h3$schedule=="even"]
cat(sprintf("T3: eps~0 loss %.1f of even %.1f = %.2f%%; eps=.3 loss %.1f of even %.1f = %.2f%%\n",
    ev0-max(h0$out[em]), ev0, 100*(ev0-max(h0$out[em]))/ev0,
    ev3-max(h3$out[em]), ev3, 100*(ev3-max(h3$out[em]))/ev3))
# (c) T6 percents: depth gains over uniform, % of uniform output
dd <- readRDS(file.path(U,"exploration_200/exploration_depth_summary.rds"))
dd <- dd[dd$tau_k==100,]; dd <- dd[order(dd$b),]
cat("T6 depth: b:", dd$b, "| S8-S2 as % of S2:", sprintf("%.2f", 100*(dd$out_S8_mean-dd$out_S2_mean)/dd$out_S2_mean),
    "| r1 share:", sprintf("%.3f", dd$alpha_S8_mean), "\n")
co <- readRDS(file.path(U,"exploration_200/exploration_corner_summary.rds"))
c1 <- co[co$tau_k==1,]; c100 <- co[co$tau_k==100,]
cat("corner b=0.5 S8-S2 %S2, tau=1 cells:", sprintf("%.2f", 100*(c1$out_S8_mean-c1$out_S2_mean)/c1$out_S2_mean),
    "| tau=100 cells:", sprintf("%.2f", 100*(c100$out_S8_mean-c100$out_S2_mean)/c100$out_S2_mean), "\n")
# (d) overtrust as fractions
dt <- readRDS(file.path(U,"D_misspecified_trust/D_misspecified_trust_summary.rds"))
for (ks_ in c(2, 1.3)) {
  cal <- dt$signal_myo_mean[dt$k_shape==ks_ & dt$tau_k_true==3 & dt$tau_k_belief==3]
  s1 <- dt$out_S1_mean[dt$k_shape==ks_ & dt$tau_k_true==3 & dt$tau_k_belief==3]
  ot <- dt[dt$k_shape==ks_ & dt$tau_k_true==3 & dt$tau_k_belief<=1,]
  cat(sprintf("k=%s: calibrated %.2f (%.2f%% of S1); overtrust cells %% of S1: %s\n", ks_, cal, 100*cal/s1,
      paste(sprintf("%.2f", 100*ot$signal_myo_mean/ot$out_S1_mean), collapse=" ")))
}
ut <- dt[dt$k_shape==2 & dt$tau_k_true==0.3,]
cat(sprintf("undertrust k=2 true=.3: calibrated %.2f; belief=3: %.2f (%.0f%%); belief=10: %.2f\n",
    ut$signal_myo_mean[ut$tau_k_belief==0.3], ut$signal_myo_mean[ut$tau_k_belief==3],
    100*ut$signal_myo_mean[ut$tau_k_belief==3]/ut$signal_myo_mean[ut$tau_k_belief==0.3],
    ut$signal_myo_mean[ut$tau_k_belief==10]))
# (e) boundary figure data: b_idx vs eps by r_min
rr <- readRDS(file.path(U,"resource_regime/resource_regime_summary.rds"))
cat("resource_regime dims: eps values:", sort(unique(rr$epsilon)), "| r_min:", sort(unique(rr$r_min)), "\n")
U <- "sweep_results"
# Overtrust figure data: value as % of S1 vs belief, true tau=3, both tails
dt <- readRDS(file.path(U,"D_misspecified_trust/D_misspecified_trust_summary.rds"))
d3 <- dt[dt$tau_k_true==3,]; d3 <- d3[order(d3$k_shape, d3$tau_k_belief),]
d3$pct <- 100*d3$signal_myo_mean/d3$out_S1_mean
d3$pct_se <- 100*d3$signal_myo_se/d3$out_S1_mean
for (r in seq_len(nrow(d3))) cat(sprintf("k=%.1f belief=%g pct=%.3f se=%.3f\n", d3$k_shape[r], d3$tau_k_belief[r], d3$pct[r], d3$pct_se[r]))
# Boundary figure data: b_idx vs eps by r_min
rr <- readRDS(file.path(U,"resource_regime/resource_regime_summary.rds"))
rr <- rr[order(rr$r_min, rr$epsilon),]
for (rm in unique(rr$r_min)) {
  s <- rr[rr$r_min==rm,]
  cat(sprintf("r_min=%g: eps %s -> b_idx %s\n", rm, paste(s$epsilon, collapse=","), paste(sprintf("%.3f", s$b_idx_S8_mean), collapse=",")))
}
# Ratio figure data: PG/signal vs T by eps
hg <- readRDS(file.path(U,"T_run_smooth/horizon_growth_summary.rds"))
hl <- readRDS(file.path(U,"T_run_smooth/horizon_long_summary.rds"))
cols <- c("T_rounds","epsilon","fwd_vs_myo_PG_mean","signal_fwd_mean")
d <- unique(rbind(hg[,cols], hl[,cols])); d <- d[d$T_rounds>=2,]
d$ratio <- pmax(d$fwd_vs_myo_PG_mean,0)/d$signal_fwd_mean
d <- d[order(d$epsilon, d$T_rounds),]
for (ep in unique(d$epsilon)) {
  s <- d[d$epsilon==ep,]
  cat(sprintf("eps=%g: T %s -> ratio %s\n", ep, paste(s$T_rounds, collapse=","), paste(sprintf("%.4f", s$ratio), collapse=",")))
}
# CD gain as % of S1
g0 <- readRDS(file.path(U,"sigma_tierA_gc0/horizon_growth_summary.rds"))
cat(sprintf("CD: max PG(T=5) %% of S1: %.4f\n", max(100*g0$fwd_vs_myo_PG_mean[g0$T_rounds==5]/g0$out_S1_mean[g0$T_rounds==5])))
g3 <- readRDS(file.path(U,"sigma_tierA_gcm3/horizon_growth_summary.rds"))
i <- g3$T_rounds==5 & g3$epsilon==0.85
cat(sprintf("gc=-3: PG(T=5,eps=.85) %% of S1: %.3f\n", 100*g3$fwd_vs_myo_PG_mean[i]/g3$out_S1_mean[i]))
# corner config check
cfgline <- grep("exploration_corner", readLines("/mnt/user-data/uploads/Grant-Funding-and-Scientific-Output-Model/sweep_T.R"))
cat(paste(readLines("/mnt/user-data/uploads/Grant-Funding-and-Scientific-Output-Model/sweep_T.R")[cfgline[1]:(cfgline[1]+8)], collapse="\n"), "\n")
suppressPackageStartupMessages(library(parallel))
source("model.R")
cores <- max(1L, detectCores() - 1L)
shares <- function(b_, sd) {
  rp <- run_simulation_T(seed = sd, T_rounds = 6, b = b_, budget_ref = "K", epsilon = 1e-4,
                         r_min = 0.001, k_shape = 1.3, tau_k = 100, allocator = "smooth",
                         strategies = c(8), M = 200)
  sp <- vapply(rp$strategies[[8]]$g_rounds, sum, numeric(1)); sp/sum(sp)
}
m <- do.call(rbind, mclapply(1:24, function(sd) shares(3, sd), mc.cores = 4))
cat(sprintf("b=3 shares: %s\nSEs: %s\n", paste(sprintf("%.4f", colMeans(m)), collapse=" "),
    paste(sprintf("%.4f", apply(m,2,sd)/sqrt(nrow(m))), collapse=" ")))
suppressPackageStartupMessages(library(parallel))
source("model.R")
gap_rule <- function(K, R, budget) {
  f <- function(cc) sum(pmax(cc * K - R, 0)) - budget
  cc <- uniroot(f, c(1e-9, 1e6), tol = 1e-10)$root
  pmax(cc * K - R, 0)
}
Tn <- 20
run1 <- function(sd) {
  rp <- run_simulation_T(seed = sd, T_rounds = Tn, tau_k = 1, k_shape = 2,
                         allocator = "smooth", strategies = c(4, 5), detail = TRUE)
  tranche <- rp$params$B_total / Tn
  out <- c()
  for (S in c(4, 5)) {
    st <- rp$strategies[[S]]
    out <- c(out, sapply(1:Tn, function(t)
      suppressWarnings(cor(st$g_rounds[[t]], gap_rule(st$K_rounds[[t]], rp$R0_at_start, tranche)))))
  }
  out
}
mm <- do.call(rbind, mclapply(1:50, function(sd) run1(sd), mc.cores = max(1L, detectCores()-1L)))
mu <- colMeans(mm, na.rm = TRUE)
cat("S4:", sprintf("%.4f", mu[1:Tn]), "\n")
cat("S5:", sprintf("%.4f", mu[Tn+1:Tn]), "\n")
