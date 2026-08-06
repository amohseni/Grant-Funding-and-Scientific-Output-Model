# Read-only verification pass: recompute every headline number cited in
# docs/PAPER_HANDOFF_2026-08-06.md directly from the canonical data files.
# Rerun anytime: Rscript sweep_results/_probe/verify_all_claims.R
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
ck <- function(lab, val, fmt = "%.4g") cat(sprintf("  %-58s %s\n", lab, sprintf(fmt, val)))

cat("== 1. resource_regime (200 seeds) ==\n")
s <- readRDS("sweep_results/resource_regime/resource_regime_summary.rds")
fl <- s[s$b_idx_S8_mean < 0.5, ]
ck("front-load cells (b_idx<0.5)", nrow(fl), "%d")
ck("max eps among front-load cells", max(fl$epsilon))
d <- s[s$epsilon == 0.85, ]; d <- d[order(d$r_min), ]
cat("  b_idx at eps=0.85, r_min 0.001..3:", sprintf("%.3f", d$b_idx_S8_mean), "\n")
cat("  round-1 share same row:", sprintf("%.3f", d$alpha_S8_mean), "\n")

cat("== 2. exploration_200: depth corner ==\n")
d <- readRDS("sweep_results/exploration_200/exploration_depth_summary.rds")
dd <- d[d$tau_k == 100, ]; dd <- dd[order(dd$b), ]
cat("  b:", dd$b, "| r1 share:", sprintf("%.3f", dd$alpha_S8_mean),
    "| S8-S2:", sprintf("%.1f", dd$out_S8_mean - dd$out_S2_mean), "\n")
allb <- do.call(rbind, lapply(c("exploration_corner","exploration_poverty","exploration_depth"),
  function(nm) readRDS(sprintf("sweep_results/exploration_200/%s_summary.rds", nm))[, c("b_idx_S8_mean","b_idx_S8_se")]))
ck("min b_idx_S8 across all 32 exploration cells", min(allb$b_idx_S8_mean))
ck("cells below 0.5 by >2 SE", sum(allb$b_idx_S8_mean < 0.5 - 2*allb$b_idx_S8_se), "%d")

cat("== 3. Package A ==\n")
s <- readRDS("sweep_results/D3_seed_signal/D3_seed_signal_summary.rds")
key <- s[s$k_shape == 1.3 & s$b == 0.5 & s$x_seed == 0.75, ]
ck("P-A1 focal cost S8 floor (raw)", -key$seed_sig_fwd_mean)
ck("P-A1 focal cost as % of S1", -100*key$seed_sig_fwd_mean/key$out_S1_mean)
kap <- -s$seed_sig_fwd_mean / ((s$out_S5_mean - s$out_S2_mean) * s$x_seed / 2)
cat("  kappa range:", sprintf("%.2f", range(kap)), "\n")
p <- readRDS("sweep_results/D4_seed_persistent/D4_seed_persistent_summary.rds")
id <- p[p$x_seed == 1, ]
ck("P-A2 max |S6-S2| at x_seed=1 (must be 0)", max(abs(id$out_S6_mean - id$out_S2_mean)))
hg <- readRDS("sweep_results/T_run_smooth/horizon_growth_summary.rds")
hl <- readRDS("sweep_results/T_run_smooth/horizon_long_summary.rds")
dd <- rbind(hg[,c("T_rounds","epsilon","fwd_vs_myo_PG_mean")], hl[,c("T_rounds","epsilon","fwd_vs_myo_PG_mean")])
for (ep in c(0.3, 0.85)) {
  d5 <- dd[dd$epsilon==ep & dd$fwd_vs_myo_PG_mean>0 & dd$T_rounds>=2 & dd$T_rounds<=5,]
  f <- lm(log(fwd_vs_myo_PG_mean) ~ log(T_rounds), d5)
  ck(sprintf("A5 exponent T<=5, eps=%.2f", ep), coef(f)[2], "%.2f")
  dl <- dd[dd$epsilon==ep & dd$fwd_vs_myo_PG_mean>0 & dd$T_rounds>=5,]
  f2 <- lm(log(fwd_vs_myo_PG_mean) ~ log(T_rounds), dl)
  ck(sprintf("A5 exponent T=5..10, eps=%.2f", ep), coef(f2)[2], "%.2f")
}

cat("== 4. Package B ==\n")
tb <- readRDS("sweep_results/sigma_tierB_concentration/sigma_tierB_summary.rds")
g21 <- tb[tb$k_shape==2 & tb$b==1, ]; g21 <- g21[order(-g21$ces_gamma), ]
cat("  gini k2 b1 (gc 0..-12):", sprintf("%.3f", g21$gini_g1_S5_mean), "\n")
lf <- readRDS("sweep_results/sigma_tierB_concentration/sigma_tierB_leontief_summary.rds")
ck("gini k2 b1 Leontief", lf$gini_g1_S5_mean[lf$k_shape==2 & lf$b==1], "%.3f")
g13 <- tb[tb$k_shape==1.3 & tb$b==0.1, ]; g13 <- g13[order(-g13$ces_gamma), ]
cat("  gini k1.3 b0.1 (gc 0..-12):", sprintf("%.3f", range(g13$gini_g1_S5_mean)),
    "| Leontief:", sprintf("%.3f", lf$gini_g1_S5_mean[lf$k_shape==1.3 & lf$b==0.1]), "\n")
rf <- readRDS("sweep_results/sigma_tierB_concentration/sigma_tierB_refine_summary.rds")
r35 <- rbind(tb[tb$k_shape==3.5 & tb$b==1, c("ces_gamma","gini_g1_S5_mean","gini_g1_S5_se")],
             rf[rf$b==1, c("ces_gamma","gini_g1_S5_mean","gini_g1_S5_se")])
pk <- r35[which.max(r35$gini_g1_S5_mean), ]
m12 <- r35[r35$ces_gamma == -12, ]
ck("interior max k3.5 b1: peak gc", pk$ces_gamma, "%g")
ck("z of peak vs gc=-12 drop", (pk$gini_g1_S5_mean - m12$gini_g1_S5_mean)/sqrt(pk$gini_g1_S5_se^2 + m12$gini_g1_S5_se^2), "%.1f")
g0 <- readRDS("sweep_results/sigma_tierA_gc0/horizon_growth_summary.rds")
ck("P-B1 CD: max PG(T=5) over eps", max(g0$fwd_vs_myo_PG_mean[g0$T_rounds==5]))
ck("P-B1 CD: b_idx range T=5", max(g0$b_idx_S8_mean[g0$T_rounds==5]), "%.3f")
g3 <- readRDS("sweep_results/sigma_tierA_gcm3/horizon_growth_summary.rds")
ck("gc=-3: PG(T=5, eps=.85)", g3$fwd_vs_myo_PG_mean[g3$T_rounds==5 & g3$epsilon==0.85], "%.2f")
ck("gc=-3: b_idx(T=5, eps=.85)", g3$b_idx_S8_mean[g3$T_rounds==5 & g3$epsilon==0.85], "%.3f")

cat("== 5. Package C ==\n")
bd <- readRDS("sweep_results/bload_decouple/bload_decouple_summary.rds")
pf <- bd[bd$eps_paid==0 & bd$eps_free>0, ]; pf <- pf[order(pf$eps_free), ]
cat("  pure-free b_idx (ef .1/.3/.85):", sprintf("%.3f", pf$b_idx_S8_mean), "\n")
pp <- bd[bd$eps_free==0 & bd$eps_paid>0, ]; pp <- pp[order(pp$eps_paid), ]
cat("  pure-paid b_idx (ep .1/.3/.85):", sprintf("%.4f", pp$b_idx_S8_mean),
    "| max SE:", sprintf("%.4f", max(pp$b_idx_S8_se)), "\n")
ck("PG at pure-paid eps_paid=.85", bd$fwd_vs_myo_PG_mean[bd$eps_free==0 & bd$eps_paid==0.85], "%.3f")
dg <- bd[bd$eps_free==bd$eps_paid & bd$eps_free %in% c(0.3, 0.85), ]
cat("  diagonal b_idx (.3,.85):", sprintf("%.4f", dg$b_idx_S8_mean[order(dg$eps_free)]),
    " vs coupled horizon_growth T=5:", sprintf("%.4f", hg$b_idx_S8_mean[hg$T_rounds==5 & hg$epsilon %in% c(0.3,0.85)]), "\n")
tr <- readRDS("sweep_results/bload_decouple/bload_transect_summary.rds")
for (ep in c(0.3, 0.85)) {
  t1 <- tr[tr$eps_paid==ep, ]; t1 <- t1[order(t1$eps_free), ]
  lo <- max(t1$eps_free[t1$b_idx_S8_mean < 0.5]); hi <- min(t1$eps_free[t1$b_idx_S8_mean >= 0.5])
  cat(sprintf("  contour @eps_paid=%.2f: eps_free in (%.3g, %.3g), ratio %.2f-%.2f\n", ep, lo, hi, lo/ep, hi/ep))
}

cat("== 6-9. Package D ==\n")
imp <- read.csv("sweep_results/D_auc_calibration/implied_tau.csv")
cat("  D-1 implied tau rows (only light tail resolvable <=20):\n"); print(imp, row.names=FALSE)
auc <- read.csv("sweep_results/D_auc_calibration/auc_grid.csv")
ck("D-1 min AUC at aK=1.3 (tau=20)", min(auc$auc[auc$alpha_K==1.3]))
dt <- readRDS("sweep_results/D_misspecified_trust/D_misspecified_trust_summary.rds")
ot <- dt[dt$k_shape==2 & dt$tau_k_true==3 & dt$tau_k_belief<=1, ]
cat("  D-2 overtrust cells k2 (belief .1/.3/1):", sprintf("%.2f", ot$signal_myo_mean[order(ot$tau_k_belief)]), "\n")
ck("D-2 heavy-tail overtrust (true3, belief .1)", dt$signal_myo_mean[dt$k_shape==1.3 & dt$tau_k_true==3 & dt$tau_k_belief==0.1], "%.2f")
ck("D-2 heavy-tail calibrated (true3, belief3)", dt$signal_myo_mean[dt$k_shape==1.3 & dt$tau_k_true==3 & dt$tau_k_belief==3], "%.2f")
ra <- readRDS("sweep_results/D_resource_ablation/D_resource_ablation_summary.rds")
w <- merge(ra[ra$use_resource_signal==TRUE,], ra[ra$use_resource_signal==FALSE,],
           by=c("rho_kr","k_shape"), suffixes=c(".T",".F"))
ck("D-3 max own-contribution % of S1", max(100*(w$out_S8_mean.T - w$out_S8_mean.F)/w$out_S1_mean.T))
gc <- read.csv("sweep_results/D_gap_convergence/gap_convergence.csv")
cat("  D-4 corr range k1.3:", sprintf("%.3f", range(gc$corr_S5[gc$k_shape==1.3])),
    "| k2:", sprintf("%.3f", range(gc$corr_S5[gc$k_shape==2])), "\n")
cat("  D-4 oracle gap %S1 (k1.3 tau20 vs tau.05):",
    sprintf("%.2f", gc$gap_S5_pctS1[gc$k_shape==1.3 & gc$tau_k %in% c(20, 0.05)]), "\n")

cat("== 10-12. Package E ==\n")
rt <- read.csv("sweep_results/E1_heterogeneity/E1a_rule_tracking.csv")
cat("  E-1a rule-tracking corr range:", sprintf("%.3f", range(rt$corr)), "\n")
eb <- read.csv("sweep_results/E1_heterogeneity/E1b_reduction.csv")
h3 <- eb[eb$tail=="heavy" & eb$tau_k==0.3, ]
cat("  E-1b signal_fwd heavy/0.3, s=1/.7/.4:", sprintf("%.2f", h3$signal_fwd[order(-h3$s)]),
    sprintf("(drop %.0f%%)", 100*(1 - min(h3$signal_fwd)/max(h3$signal_fwd))), "\n")
rc <- read.csv("sweep_results/E1_heterogeneity/E1b_rankcorr.csv")
cat("  E-1b rank-corr range:", sprintf("%.3f", range(rc$rank_corr_vs_s1)), "\n")
n2 <- readRDS("sweep_results/E2_headline_scale_n200/horizon_growth_summary.rds")
m <- merge(n2[,c("T_rounds","epsilon","fwd_vs_myo_PG_mean","b_idx_S8_mean","out_S1_mean")],
           hg[,c("T_rounds","epsilon","fwd_vs_myo_PG_mean","b_idx_S8_mean","out_S1_mean")],
           by=c("T_rounds","epsilon"), suffixes=c(".200",".50"))
ck("E-2 corr PG% n200 vs n50",
   cor(100*m$fwd_vs_myo_PG_mean.200/m$out_S1_mean.200, 100*m$fwd_vs_myo_PG_mean.50/m$out_S1_mean.50), "%.3f")
ck("E-2 corr b_idx", cor(m$b_idx_S8_mean.200, m$b_idx_S8_mean.50), "%.3f")

cat("== 13. Bootstrap ==\n")
h0 <- read.csv("sweep_results/bootstrap_verify/honest_schedules_eps0.csv")
h3b <- read.csv("sweep_results/bootstrap_verify/honest_schedules_eps03.csv")
ev0 <- h0$out[h0$schedule=="even"]
ck("item1 eps~0: even output", ev0, "%.1f")
xs <- suppressWarnings(as.numeric(sub("x([0-9.]+)_.*","\\1", h0$schedule)))
ks <- suppressWarnings(as.numeric(sub(".*_k([0-9]+)","\\1", h0$schedule)))
em <- which(!is.na(xs) & xs > ks/6 + 1e-9)   # early-mass: round-1..k share above even
ck("item1 eps~0: best early-mass minus even", max(h0$out[em]) - ev0, "%.1f")
ev3 <- h3b$out[h3b$schedule=="even"]
ck("item1 eps=0.3: best early-mass minus even", max(h3b$out[em]) - ev3, "%.1f")
cc <- read.csv("sweep_results/bootstrap_verify/ce_self_consistency.csv")
ck("item2 CE(S8)-CE(S5) mean (n=20)", mean(cc$ce_S8 - cc$ce_S5), "%.2f")
ck("item2 seeds with CE(S8)>=CE(S5)", sum(cc$ce_S8 >= cc$ce_S5), "%d")
for (nm in c("exploration_corner","exploration_poverty","exploration_depth")) {
  a <- readRDS(sprintf("sweep_results/exploration_200/%s_summary.rds", nm))
  b <- readRDS(sprintf("sweep_results/T_run_smooth_supplement/%s_summary.rds", nm))
  key <- intersect(names(a), c("tau_k","epsilon","r_min","b"))
  mm <- merge(a[,c(key,"b_idx_S8_mean")], b[,c(key,"b_idx_S8_mean")], by=key)
  ck(sprintf("item3 %s max |d b_idx| 200 vs 64", nm), max(abs(mm[,ncol(mm)]-mm[,ncol(mm)-1])))
}
cat("\nVERIFICATION PASS COMPLETE\n")
