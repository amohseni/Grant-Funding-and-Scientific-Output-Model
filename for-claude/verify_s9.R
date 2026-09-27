U <- "sweep_results"
cat("== RUN_INFOs ==\n")
for (f in c("sigma_tierB_concentration/RUN_INFO.txt","sigma_tierA_gc0/RUN_INFO.txt","sigma_tierA_leontief/RUN_INFO.txt"))
  { cat("--", f, "\n"); cat(readLines(file.path(U,f)), sep="\n"); cat("\n") }
cat("\n== Tier B concentration (T=1, S5): Gini of optimal-ish allocation by gamma ==\n")
tb <- readRDS(file.path(U,"sigma_tierB_concentration/sigma_tierB_summary.rds"))
cat("cols:", paste(head(names(tb), 20), collapse=","), "...\n")
cat("grid: gamma", sort(unique(tb$ces_gamma)), "| k_shape", sort(unique(tb$k_shape)), "| b", sort(unique(tb$b)), "\n")
lf <- readRDS(file.path(U,"sigma_tierB_concentration/sigma_tierB_leontief_summary.rds"))
rf <- readRDS(file.path(U,"sigma_tierB_concentration/sigma_tierB_refine_summary.rds"))
cat("refine grid: gamma", sort(unique(rf$ces_gamma)), "| k", sort(unique(rf$k_shape)), "| b", sort(unique(rf$b)), "\n")
for (ks in sort(unique(tb$k_shape))) for (b_ in sort(unique(tb$b))) {
  r <- tb[tb$k_shape==ks & tb$b==b_,]; r <- r[order(-r$ces_gamma),]
  lrow <- lf[lf$k_shape==ks & lf$b==b_,]
  cat(sprintf("k=%.1f b=%.1f: gc %s -> gini %s | Leontief %.3f (se %.3f)\n", ks, b_,
      paste(r$ces_gamma, collapse=","), paste(sprintf("%.3f", r$gini_g1_S5_mean), collapse=","),
      lrow$gini_g1_S5_mean, lrow$gini_g1_S5_se))
}
cat("refine rows (k=3.5):\n")
rr <- rf[order(rf$b, -rf$ces_gamma),]
for (i in seq_len(nrow(rr))) cat(sprintf("  b=%.1f gc=%g gini=%.4f se=%.4f\n", rr$b[i], rr$ces_gamma[i], rr$gini_g1_S5_mean[i], rr$gini_g1_S5_se[i]))
cat("\n== Tier A signal value across the family ==\n")
for (nm in c("sigma_tierA_gc0","sigma_tierA_gcm3","sigma_tierA_leontief")) {
  sv <- readRDS(file.path(U, nm, "signal_value_summary.rds"))
  cat("--", nm, ": cols", paste(intersect(c("k_shape","b","tau_k","ces_gamma","signal_myo_mean","out_S1_mean"), names(sv)), collapse=","), "\n")
  sv <- sv[order(sv$k_shape, sv$tau_k),]
  for (i in seq_len(nrow(sv))) cat(sprintf("  k=%.1f tau=%g b=%s: signal %%S1 = %.2f (se %.2f)\n",
      sv$k_shape[i], sv$tau_k[i], if("b" %in% names(sv)) sv$b[i] else "?",
      100*sv$signal_myo_mean[i]/sv$out_S1_mean[i], 100*sv$signal_myo_se[i]/sv$out_S1_mean[i]))
}
cat("\n== Tier A seed cost (no-signal family) gc0 ==\n")
sd0 <- readRDS(file.path(U,"sigma_tierA_gc0/seed_value_summary.rds"))
cat("cols:", paste(intersect(c("k_shape","b","x_seed","seed_myo_mean","out_S1_mean"), names(sd0)), collapse=","), "\n")
cat(sprintf("worst seed cost %%S1: %.3f\n", min(100*sd0$seed_myo_mean/sd0$out_S1_mean)))
