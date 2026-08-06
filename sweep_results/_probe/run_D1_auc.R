# D-1 (addendum): model-implied AUC of the peer-review signal for identifying
# top-q researchers by K. Pure Monte Carlo, no simulation engine.
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
outdir <- "sweep_results/D_auc_calibration"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(20260805)
N <- 1e6
grid <- expand.grid(alpha_K = c(1.3, 2, 3.5),
                    tau_k   = c(0.05, 0.1, 0.3, 0.5, 1, 1.5, 2.5, 5, 10, 20),
                    q       = c(0.1, 0.2))
auc_one <- function(aK, tk, q) {
  K <- 1 / runif(N)^(1 / aK)                 # Pareto(k_min=1, shape=aK)
  sig <- K + rnorm(N, 0, tk)
  pos <- K >= quantile(K, 1 - q)
  r <- rank(sig)                             # AUC via rank-sum
  n1 <- sum(pos); n0 <- N - n1
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
grid$auc <- mapply(auc_one, grid$alpha_K, grid$tau_k, grid$q)
write.csv(grid, file.path(outdir, "auc_grid.csv"), row.names = FALSE)

# implied tau_k at AUC = 0.54, per (alpha_K, q): log-linear interpolation + bracket
cat("== implied tau_k at AUC = 0.54 ==\n")
imp <- list()
for (aK in unique(grid$alpha_K)) for (qq in unique(grid$q)) {
  d <- grid[grid$alpha_K == aK & grid$q == qq, ]; d <- d[order(d$tau_k), ]
  above <- d$auc > 0.54
  if (all(above)) { cat(sprintf("aK=%.1f q=%.1f: AUC>0.54 at all tau<=20 (implied tau > 20)\n", aK, qq)); next }
  i <- max(which(above))
  f <- (0.54 - d$auc[i]) / (d$auc[i+1] - d$auc[i])
  tk_imp <- exp(log(d$tau_k[i]) + f * (log(d$tau_k[i+1]) - log(d$tau_k[i])))
  imp[[length(imp)+1]] <- data.frame(alpha_K = aK, q = qq, tau_implied = tk_imp,
                                     bracket_lo = d$tau_k[i], bracket_hi = d$tau_k[i+1])
  cat(sprintf("aK=%.1f q=%.1f: implied tau_k = %.1f (bracket %.1f-%.1f)\n",
              aK, qq, tk_imp, d$tau_k[i], d$tau_k[i+1]))
}
imp <- do.call(rbind, imp)
write.csv(imp, file.path(outdir, "implied_tau.csv"), row.names = FALSE)

# figure
suppressPackageStartupMessages(library(ggplot2))
grid$lab <- sprintf("alpha_K=%.1f, q=%.1f", grid$alpha_K, grid$q)
p <- ggplot(grid, aes(tau_k, auc, color = factor(alpha_K), linetype = factor(q))) +
  geom_hline(yintercept = 0.54, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 2.5, linetype = "dotted", color = "grey60") +
  annotate("text", x = 2.5, y = 0.95, label = "signal-value half-decay (elbow)", size = 3, hjust = -0.05, color = "grey40") +
  annotate("text", x = 0.05, y = 0.545, label = "Fang et al. AUC = 0.54", size = 3, hjust = 0, vjust = -0.5, color = "grey40") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.5) + scale_x_log10() +
  labs(x = expression(tau[K]~"(log)"), y = "AUC (top-q identification)",
       color = expression(alpha[K]), linetype = "q",
       title = "Model-implied AUC of the peer-review signal vs noise",
       subtitle = "Where does Fang et al.'s AUC = 0.54 place real NIH review on the precision axis?") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "auc_vs_tau.png"), p, width = 8.5, height = 5.5, dpi = 130)
writeLines(c(paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
             paste("date:", format(Sys.time())), "MC: N=1e6/cell, seed 20260805"),
           file.path(outdir, "RUN_INFO.txt"))
cat("D-1 done\n")
