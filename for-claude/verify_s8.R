U <- "sweep_results"
ck <- function(lab, val, fmt = "%.4g") cat(sprintf("  %-58s %s\n", lab, sprintf(fmt, val)))
cat("== Package A: seed floors ==\n")
s <- readRDS(file.path(U,"D3_seed_signal/D3_seed_signal_summary.rds"))
cat("cols:", paste(names(s), collapse=","), "\n")
cat("grid: k_shape", sort(unique(s$k_shape)), "| b", sort(unique(s$b)), "| x_seed", sort(unique(s$x_seed)), "\n")
key <- s[s$k_shape == 1.3 & s$b == 0.5 & s$x_seed == 0.75, ]
ck("P-A1 focal cost (S11 vs S8, raw)", -key$seed_sig_fwd_mean)
ck("P-A1 focal cost as % of S1", -100*key$seed_sig_fwd_mean/key$out_S1_mean)
kap <- -s$seed_sig_fwd_mean / ((s$out_S5_mean - s$out_S2_mean) * s$x_seed / 2)
cat("  kappa range:", sprintf("%.2f", range(kap)), "\n")
# cost as % of S1 across cells
s$cost_pct <- -100*s$seed_sig_fwd_mean/s$out_S1_mean
for (ks in sort(unique(s$k_shape))) for (b_ in sort(unique(s$b))) {
  r <- s[s$k_shape==ks & s$b==b_,]; r <- r[order(r$x_seed),]
  cat(sprintf("  k=%.1f b=%.1f: x_seed %s -> cost %%S1 %s | S5-S2 %%S1 %s\n", ks, b_,
      paste(r$x_seed, collapse=","), paste(sprintf("%.2f", r$cost_pct), collapse=","),
      paste(sprintf("%.2f", 100*(r$out_S5_mean-r$out_S2_mean)/r$out_S1_mean), collapse=",")))
}
p <- readRDS(file.path(U,"D4_seed_persistent/D4_seed_persistent_summary.rds"))
cat("D4 cols:", paste(names(p), collapse=","), "\n")
id <- p[p$x_seed == 1, ]
ck("P-A2 max |S6-S2| at x_seed=1 (must be 0)", max(abs(id$out_S6_mean - id$out_S2_mean)))
pp <- p[order(p$k_shape, p$x_seed),]
for (ks in sort(unique(p$k_shape))) {
  r <- pp[pp$k_shape==ks,]
  cat(sprintf("  D4 k=%.1f: x_seed %s -> (S6-S2)/S1 %% %s\n", ks, paste(r$x_seed, collapse=","),
      paste(sprintf("%.2f", 100*(r$out_S6_mean-r$out_S2_mean)/r$out_S1_mean), collapse=",")))
}
U <- "sweep_results"
s <- readRDS(file.path(U,"D3_seed_signal/D3_seed_signal_summary.rds"))
s <- s[s$x_seed == 0.1,]  # one row per (k_shape, b); strategy outputs identical across x_seed
for (ks in sort(unique(s$k_shape))) { r <- s[s$k_shape==ks,]; r <- r[order(r$b),]
  cat(sprintf("k=%.1f: b %s | uniform gain (S2-S1)/S1 %% %s | targeting value rel to uniform gain (S5-S2)/(S2-S1) %s | (S4-S2)/S1 %% %s\n",
    ks, paste(r$b, collapse=","),
    paste(sprintf("%.1f", 100*(r$out_S2_mean-r$out_S1_mean)/r$out_S1_mean), collapse=","),
    paste(sprintf("%.2f", (r$out_S5_mean-r$out_S2_mean)/(r$out_S2_mean-r$out_S1_mean)), collapse=","),
    paste(sprintf("%.2f", 100*(r$out_S4_mean-r$out_S2_mean)/r$out_S1_mean), collapse=",")))
}
cat("tau_k in D3:", unique(s$tau_k_mean), "\n")
