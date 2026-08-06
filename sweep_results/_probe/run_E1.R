# E-1 (addendum): heterogeneous productivity.
# E-1a: A observable — direction sweep + scaled-gap-rule tracking.
# E-1b: A latent — the product-reduction test (variance split s of log T).
suppressPackageStartupMessages(library(parallel))
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
t0 <- Sys.time()
cores <- max(1L, detectCores() - 1L)

# ---------- E-1a directions (manifest machinery) ----------
SWEEP_CONFIGS_T$E1a_hetA <- list(
  name = "E1a_hetA", tier = 1,
  description = "Observable heterogeneous productivity A_i (lognormal, E[A]=1): do headline directions survive?",
  grid_fn = function() expand.grid(A_obs_sdlog = c(0, 0.25, 0.5),
                                   k_shape = c(1.3, 3.5), tau_k = c(0.3, 10)),
  varied_params = c("A_obs_sdlog", "k_shape", "tau_k"),
  primary_plot = list(type = "line", x_var = "A_obs_sdlog", y_var = "signal_fwd_mean",
                      color_var = "k_shape", title = "Signal value under heterogeneous productivity",
                      y_label = "S8 - S7"),
  secondary_plot = list(type = "line", x_var = "A_obs_sdlog", y_var = "seed_fwd_mean",
                        color_var = "k_shape", title = "Seed cost under heterogeneous productivity",
                        y_label = "S9 - S7"))
sweep_one_T("E1a_hetA", seeds = 1:100,
            base_params = list(allocator = "smooth", strategies = c(1,2,5,7,8,9)),
            out_dir = "sweep_results/E1_heterogeneity", resume = TRUE)
cat(sprintf("[E1a directions done] %.0fs\n", as.numeric(Sys.time()-t0, units="secs")))

# ---------- E-1a rule tracking: corr(g1_S5, scaled gap rule) ----------
scaled_rule <- function(K, R, gam_i, budget) {
  f <- function(nu) sum(pmax(K * (sqrt(gam_i / nu) - 1) - R, 0)) - budget
  hi <- max(gam_i * (K / (K + R))^2) * 0.999999
  nu <- uniroot(f, c(hi * 1e-12, hi), tol = 1e-12)$root
  pmax(K * (sqrt(gam_i / nu) - 1) - R, 0)
}
cells <- expand.grid(sdlog = c(0.25, 0.5), k_shape = c(1.3, 2), tau_k = c(0.3, 1))
track <- do.call(rbind, lapply(seq_len(nrow(cells)), function(ci) {
  cl <- cells[ci, ]
  rows <- mclapply(1:30, function(sd) {
    r <- run_simulation_T(seed = sd, T_rounds = 2, k_shape = cl$k_shape, tau_k = cl$tau_k,
                          A_obs_sdlog = cl$sdlog, allocator = "smooth",
                          strategies = c(5), detail = TRUE)
    K <- r$K_at_start; R <- r$R0_at_start; gam_i <- r$params$gamma * r$A_i
    gstar <- scaled_rule(K, R, gam_i, r$params$B_total / 2)
    suppressWarnings(cor(r$strategies[[5]]$g_rounds[[1]], gstar))
  }, mc.cores = cores)
  data.frame(cl, corr = mean(unlist(rows), na.rm = TRUE), corr_se = sd(unlist(rows), na.rm=TRUE)/sqrt(30))
}))
write.csv(track, "sweep_results/E1_heterogeneity/E1a_rule_tracking.csv", row.names = FALSE)
cat("== E-1a rule tracking ==\n"); print(track, digits = 3, row.names = FALSE)

# ---------- E-1b: the product-reduction test ----------
tails <- list(heavy = list(V = (1/1.3)^2, meanlogT = 1/1.3),
              moderate = list(V = (1/2)^2, meanlogT = 1/2))
grid <- expand.grid(s = c(1, 0.7, 0.4), tail = names(tails), tau_k = c(0.3, 1),
                    stringsAsFactors = FALSE)
res_rows <- list(); gvecs <- list()
for (ci in seq_len(nrow(grid))) {
  g <- grid[ci, ]; tl <- tails[[g$tail]]
  rows <- mclapply(1:100, function(sd) {
    r <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = g$tau_k, allocator = "smooth",
                          strategies = c(1,2,5,7,8),
                          latentA_cfg = list(split_s = g$s, V = tl$V, meanlogT = tl$meanlogT))
    o <- function(x) r$strategies[[x]]$total_expected
    list(m = c(out_S8 = o(8), signal_fwd = o(8) - o(7), vs_unif = o(8) - o(2),
               gini_g1 = gini(r$strategies[[8]]$g_rounds[[1]]), out_S1 = o(1)),
         g1 = r$strategies[[8]]$g_rounds[[1]])
  }, mc.cores = cores)
  m <- do.call(rbind, lapply(rows, `[[`, "m"))
  res_rows[[ci]] <- data.frame(g, t(colMeans(m)),
                               out_S8_se = sd(m[,"out_S8"])/10, signal_fwd_se = sd(m[,"signal_fwd"])/10)
  gvecs[[paste(g$s, g$tail, g$tau_k)]] <- lapply(rows, `[[`, "g1")
}
tab <- do.call(rbind, res_rows)
# rank correlation of S8 grant vectors across s (same seed pairing)
rc <- do.call(rbind, lapply(names(tails), function(tn) do.call(rbind, lapply(c(0.3, 1), function(tk) {
  base <- gvecs[[paste(1, tn, tk)]]
  do.call(rbind, lapply(c(0.7, 0.4), function(ss) {
    alt <- gvecs[[paste(ss, tn, tk)]]
    r <- mean(mapply(function(a, b) suppressWarnings(cor(a, b, method = "spearman")), base, alt), na.rm = TRUE)
    data.frame(tail = tn, tau_k = tk, s = ss, rank_corr_vs_s1 = r)
  }))
}))))
write.csv(tab, "sweep_results/E1_heterogeneity/E1b_reduction.csv", row.names = FALSE)
write.csv(rc,  "sweep_results/E1_heterogeneity/E1b_rankcorr.csv", row.names = FALSE)
cat("== E-1b reduction test ==\n"); print(tab, digits = 4, row.names = FALSE)
cat("== E-1b rank corr of S8 grants vs s=1 ==\n"); print(rc, digits = 3, row.names = FALSE)
writeLines(c(paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
             paste("date:", format(Sys.time())),
             "E1a: 100 seeds directions + 30 seeds rule tracking; E1b: 100 seeds, CRN across s"),
           "sweep_results/E1_heterogeneity/RUN_INFO.txt")
cat(sprintf("E-1 DONE in %.0f min\n", as.numeric(Sys.time()-t0, units="mins")))
