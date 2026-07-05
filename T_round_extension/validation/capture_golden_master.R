# ============================================================
# tests/capture_golden_master.R
# ------------------------------------------------------------
# Snapshot current v5 (app.R) per-strategy behavior as the
# ground-truth regression anchor for the T-round extension.
#
# The new T-round planner, evaluated at T=2, must reproduce
# these numbers bit-for-bit (algorithmic-regression anchor).
# Output is normalizer-independent: run_simulation fixes the
# absolute budget B internally, so the later switch to
# median/b=0.1 normalization does not invalidate this file.
#
# Output: tests/golden_master_v5.rds
#   list of rows, one per (seed x param-point), each holding
#   per-strategy g1, g2, total_expected, alpha.
# ============================================================

setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")

# --- Source app.R without launching Shiny ---
src <- readLines("app.R", warn = FALSE)
end <- grep("^shinyApp\\(", src)[1] - 1L
eval(parse(text = paste(src[1:end], collapse = "\n")), envir = globalenv())

# --- Parameter points to snapshot (small, fast, representative) ---
# Include the CURRENT default plus off-default points that exercise
# each force: signal precision (tau_k), growth (epsilon), tail (k_shape),
# budget (b), seed floor (x_seed).
points <- list(
  list(tag = "default",      n = 30, tau_k = 1.0, epsilon = 0.3, k_shape = 2.0, b = 0.5, x_seed = 0.5),
  list(tag = "sharp_signal", n = 30, tau_k = 0.3, epsilon = 0.3, k_shape = 2.0, b = 0.5, x_seed = 0.5),
  list(tag = "noisy_signal", n = 30, tau_k = 10.0, epsilon = 0.3, k_shape = 2.0, b = 0.5, x_seed = 0.5),
  list(tag = "low_growth",   n = 30, tau_k = 1.0, epsilon = 0.05, k_shape = 2.0, b = 0.5, x_seed = 0.5),
  list(tag = "heavy_tail",   n = 30, tau_k = 1.0, epsilon = 0.3, k_shape = 1.3, b = 0.5, x_seed = 0.5),
  list(tag = "small_budget", n = 30, tau_k = 1.0, epsilon = 0.3, k_shape = 2.0, b = 0.1, x_seed = 0.25)
)
seeds <- 1:5

# --- Capture ---
golden <- list()
for (pt in points) {
  for (sd in seeds) {
    res <- run_simulation(
      seed = sd, n = pt$n,
      k_min = 1.0, k_shape = pt$k_shape,
      r_min = 1.0, r_shape = 2.0,
      gamma = 1.0, epsilon = pt$epsilon,
      b = pt$b, n_steps = 50,
      tau_r = 1.0, tau_k = pt$tau_k,
      x_seed = pt$x_seed, M = 500,
      strategies = 1:9, verbose = FALSE
    )
    row <- list(tag = pt$tag, seed = sd, params = pt)
    for (s in 1:9) {
      r <- res$strategies[[s]]
      row[[paste0("g1_S", s)]]  <- r$g1
      row[[paste0("g2_S", s)]]  <- r$g2
      row[[paste0("out_S", s)]] <- r$total_expected
      row[[paste0("alpha_S", s)]] <- if (is.null(r$alpha)) NA_real_ else r$alpha
    }
    golden[[length(golden) + 1L]] <- row
  }
}

dir.create("tests", showWarnings = FALSE)
saveRDS(golden, "tests/golden_master_v5.rds")

# --- Console summary so we can eyeball the anchor contrasts ---
cat("Captured", length(golden), "golden rows ->  tests/golden_master_v5.rds\n\n")
cat(sprintf("%-14s %5s %10s %10s %10s %10s\n",
            "point", "seed", "signal_fwd", "fwd_v_myoPG", "seed_fwd", "alpha_S8"))
for (row in golden) {
  signal_fwd    <- row$out_S8 - row$out_S7
  fwd_vs_myo_PG <- row$out_S8 - row$out_S5
  seed_fwd      <- row$out_S9 - row$out_S7
  cat(sprintf("%-14s %5d %10.3f %10.3f %10.3f %10.3f\n",
              row$tag, row$seed, signal_fwd, fwd_vs_myo_PG, seed_fwd, row$alpha_S8))
}
