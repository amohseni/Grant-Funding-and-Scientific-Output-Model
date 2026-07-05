# ============================================================
# tests/test_T2_reduction.R
# ------------------------------------------------------------
# Regression gate: run_simulation_T(T=2) must reproduce v5
# run_simulation() at identical params. Reports max abs
# difference in per-strategy total_expected, round-1 grant (g1),
# and alpha, across several parameter points and seeds.
# PASS if bit-identical (< 1e-8); WARN if only within MC error.
# ============================================================
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")

src <- readLines("app.R", warn = FALSE)
end <- grep("^shinyApp\\(", src)[1] - 1L
eval(parse(text = paste(src[1:end], collapse = "\n")), envir = globalenv())
source("simulate_T.R")

points <- list(
  list(tag = "default",      n = 30, tau_k = 1.0,  epsilon = 0.3,  k_shape = 2.0, b = 0.5, x_seed = 0.5),
  list(tag = "sharp_signal", n = 30, tau_k = 0.3,  epsilon = 0.3,  k_shape = 2.0, b = 0.5, x_seed = 0.5),
  list(tag = "noisy_signal", n = 30, tau_k = 10.0, epsilon = 0.3,  k_shape = 2.0, b = 0.5, x_seed = 0.5),
  list(tag = "low_growth",   n = 30, tau_k = 1.0,  epsilon = 0.05, k_shape = 2.0, b = 0.5, x_seed = 0.5),
  list(tag = "heavy_tail",   n = 30, tau_k = 1.0,  epsilon = 0.3,  k_shape = 1.3, b = 0.5, x_seed = 0.5),
  list(tag = "small_budget", n = 30, tau_k = 1.0,  epsilon = 0.3,  k_shape = 2.0, b = 0.1, x_seed = 0.25),
  list(tag = "with_prerounds", n = 30, tau_k = 1.0, epsilon = 0.3, k_shape = 2.0, b = 0.5, x_seed = 0.5, n_pre = 3)
)
seeds <- 1:5

max_out <- 0; max_g1 <- 0; max_alpha <- 0
worst <- ""
for (pt in points) {
  npre <- if (is.null(pt$n_pre)) 0 else pt$n_pre
  for (sd in seeds) {
    args <- list(seed = sd, n = pt$n, k_min = 1.0, k_shape = pt$k_shape,
                 r_min = 1.0, r_shape = 2.0, gamma = 1.0, epsilon = pt$epsilon,
                 b = pt$b, n_steps = 50, tau_r = 1.0, tau_k = pt$tau_k,
                 n_pre_rounds = npre, x_seed = pt$x_seed, M = 500, strategies = 1:9)
    v5 <- do.call(run_simulation, args)
    tt <- do.call(run_simulation_T, c(args, list(T_rounds = 2)))
    for (s in 1:9) {
      d_out <- abs(v5$strategies[[s]]$total_expected - tt$strategies[[s]]$total_expected)
      d_g1  <- max(abs(v5$strategies[[s]]$g1 - tt$strategies[[s]]$g_rounds[[1]]))
      a5 <- v5$strategies[[s]]$alpha; at <- tt$strategies[[s]]$alpha
      d_al <- if (is.na(a5) && is.na(at)) 0 else
              if (is.na(a5) != is.na(at)) Inf else abs(a5 - at)
      if (d_out > max_out) { max_out <- d_out; worst <- sprintf("%s seed%d S%d", pt$tag, sd, s) }
      max_g1 <- max(max_g1, d_g1); max_alpha <- max(max_alpha, d_al)
    }
  }
}

cat(sprintf("\n=== T=2 reduction vs v5 (7 points x 5 seeds x 9 strategies) ===\n"))
cat(sprintf("max |Δ total_expected| = %.3e   (worst: %s)\n", max_out, worst))
cat(sprintf("max |Δ g1|             = %.3e\n", max_g1))
cat(sprintf("max |Δ alpha|          = %.3e\n", max_alpha))
tol <- 1e-8
verdict <- if (max(max_out, max_g1, max_alpha) < tol) "PASS (bit-identical)" else
           if (max(max_out, max_g1) < 1e-3) "WARN (within MC error, not bit-identical)" else "FAIL"
cat(sprintf("VERDICT: %s\n", verdict))
