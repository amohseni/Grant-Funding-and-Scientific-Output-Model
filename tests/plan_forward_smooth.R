# ============================================================
# Phase 2 — continuous FORWARD planner (projected-gradient prototype)
# ============================================================
# Drop-in replacement for the greedy inner loop of plan_forward_ce()
# (model.R:1017) in the FORWARD path (S7/S8/S9). Removes the greedy step
# delta entirely: it solves the continuous budget-constrained optimum of the
# receding-horizon certainty-equivalent (CE) objective.
#
# Objective (identical to plan_forward_ce's, summed over researchers):
#   V(G) = sum_i fwd_researcher_value(post0_i, bar1_i(G[i,1]), g_hist_i, G[i,])
# where G is the n x H grant plan, bar1_i = ce_reweight_posterior(post0_i, G[i,1])
# is the info-anticipating round-1 posterior, and each researcher's value
# compounds K along its own path. The objective is SEPARABLE across researchers
# and coupled only through the budget:  sum(G) = B_rem,  G >= L  (L = seed floor).
#
# Method: projected-gradient ascent. Gradient via CENTRAL finite differences of
# fwd_researcher_value (the faithful objective; fwd_marginal is its forward-diff
# gradient, central diff of the same value is strictly more accurate). Projection
# onto {G >= L, sum G = B_rem} is Euclidean simplex projection (Duchi et al.),
# O(nH log nH). Backtracking line search on V; monotone ascent.
#
# STATUS: prototype. NOT yet integrated into model.R. Validation harness at bottom
# (Rscript tests/plan_forward_smooth.R) checks it matches/dominates the fine-greedy
# plan_forward_ce limit and satisfies KKT (funded slots share a common marginal).

# ----- Euclidean projection onto {x >= L, sum x = z} -----
# Projects v onto the shifted simplex. If z <= sum(L) (nothing to give beyond the
# floor), returns L. Standard sort-based threshold on the shifted vector d = v - L.
proj_simplex_lb <- function(v, z, L) {
  d  <- v - L
  zz <- z - sum(L)
  if (zz <= 0) return(L)
  u   <- sort(d, decreasing = TRUE)
  css <- cumsum(u)
  j   <- seq_along(u)
  rho <- max(which(u - (css - zz) / j > 0))
  theta <- (css[rho] - zz) / rho
  L + pmax(d - theta, 0)
}

# ----- Gradient of one researcher's H-round value w.r.t. its plan -----
# Central differences of fwd_researcher_value. Slot s=1 also moves bar1 (the CE
# info term), so it is re-reweighted at g1 +/- dg; slots s>=2 keep bar1 fixed.
# Falls back to a forward difference for a slot pinned at its lower bound (0),
# where a central step would go negative.
researcher_grad <- function(post0, g_hist_i, g_path, dg, gamma, epsilon,
                            lo = rep(0, length(g_path))) {
  H <- length(g_path)
  grad <- numeric(H)
  bar1_fixed <- ce_reweight_posterior(post0, g_path[1], gamma)  # for s >= 2
  for (s in seq_len(H)) {
    central <- (g_path[s] - dg) >= lo[s]
    gp <- g_path; gp[s] <- gp[s] + dg
    gm <- g_path; if (central) gm[s] <- gm[s] - dg
    if (s == 1L) {
      vp <- fwd_researcher_value(post0, ce_reweight_posterior(post0, gp[1], gamma),
                                 g_hist_i, gp, gamma, epsilon)
      vm <- fwd_researcher_value(post0, ce_reweight_posterior(post0, gm[1], gamma),
                                 g_hist_i, gm, gamma, epsilon)
    } else {
      vp <- fwd_researcher_value(post0, bar1_fixed, g_hist_i, gp, gamma, epsilon)
      vm <- fwd_researcher_value(post0, bar1_fixed, g_hist_i, gm, gamma, epsilon)
    }
    grad[s] <- (vp - vm) / (if (central) 2 * dg else dg)
  }
  grad
}

# ----- Continuous forward planner -----
# Signature mirrors plan_forward_ce (minus `delta`). Returns an n x H plan; the
# caller executes column 1 only. H == 1 delegates to the single-round water-fill.
plan_forward_smooth <- function(post0_list, B_rem, gamma, epsilon, H,
                                g_hist = list(), seed_round1 = NULL,
                                dg = 1e-4, max_iter = 400, tol = 1e-10,
                                verbose = FALSE) {
  n   <- length(post0_list)
  ghi <- function(i) vapply(g_hist, function(gh) gh[i], numeric(1))

  if (H == 1L) {
    seed_sum <- if (is.null(seed_round1)) 0 else sum(seed_round1)
    g <- waterfill_round_t(post0_list, B_rem - seed_sum, gamma, epsilon,
                           g_hist = g_hist, init_g = seed_round1)
    return(matrix(g, nrow = n, ncol = 1))
  }

  L <- matrix(0, n, H)
  if (!is.null(seed_round1)) L[, 1] <- seed_round1
  if (B_rem <= sum(L) + tol) return(L)

  recompute_bar1 <- function(plan)
    lapply(seq_len(n), function(i) ce_reweight_posterior(post0_list[[i]], plan[i, 1], gamma))
  value <- function(plan, bar1) {
    v <- 0
    for (i in seq_len(n))
      v <- v + fwd_researcher_value(post0_list[[i]], bar1[[i]], ghi(i), plan[i, ], gamma, epsilon)
    v
  }
  grad <- function(plan) {
    G <- matrix(0, n, H)
    for (i in seq_len(n))
      G[i, ] <- researcher_grad(post0_list[[i]], ghi(i), plan[i, ], dg, gamma, epsilon,
                                lo = L[i, ])
    G
  }

  # init: floor + remaining spread uniformly over all slots
  plan <- L + (B_rem - sum(L)) / (n * H)
  bar1 <- recompute_bar1(plan); V <- value(plan, bar1)

  eta <- NA_real_; prevV <- V
  for (iter in seq_len(max_iter)) {
    G <- grad(plan)
    if (is.na(eta)) eta <- 0.5 * (B_rem / (n * H)) / max(abs(G), 1e-12)
    accepted <- FALSE
    for (bt in 1:50) {
      cand <- matrix(proj_simplex_lb(as.vector(plan + eta * G), B_rem, as.vector(L)), n, H)
      bar1c <- recompute_bar1(cand); Vc <- value(cand, bar1c)
      if (Vc > V + 1e-13 * abs(V)) {
        plan <- cand; bar1 <- bar1c; V <- Vc; accepted <- TRUE; eta <- eta * 1.3
        break
      }
      eta <- eta * 0.5
    }
    if (verbose) cat(sprintf("  iter %3d  V=%.9f  eta=%.3e  %s\n",
                             iter, V, eta, if (accepted) "ok" else "stall"))
    if (!accepted) break
    if (abs(V - prevV) < tol * max(1, abs(V))) break
    prevV <- V
  }
  plan
}

# ------------------------------------------------------------
# VALIDATION HARNESS (run: Rscript tests/plan_forward_smooth.R)
# Confirms: smooth >= fine-greedy plan_forward_ce limit, spends budget, KKT.
# ------------------------------------------------------------
if (sys.nframe() == 0L) {
  root <- "/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model"
  setwd(root); source("model.R"); source("tests/waterfill_round_t.R")

  # objective evaluator (matches plan_forward_ce's objective for any plan)
  plan_obj <- function(posts, plan, gamma, epsilon, g_hist = list()) {
    n <- length(posts); ghi <- function(i) vapply(g_hist, function(gh) gh[i], numeric(1))
    v <- 0
    for (i in seq_len(n)) {
      bar1 <- ce_reweight_posterior(posts[[i]], plan[i, 1], gamma)
      v <- v + fwd_researcher_value(posts[[i]], bar1, ghi(i), plan[i, ], gamma, epsilon)
    }
    v
  }
  # marginal spread across funded slots (KKT: all equal at the optimum)
  kkt_spread <- function(posts, plan, gamma, epsilon, dg = 1e-6, thresh = 1e-6) {
    n <- length(posts); H <- ncol(plan); m <- c()
    for (i in seq_len(n)) {
      bar1 <- ce_reweight_posterior(posts[[i]], plan[i, 1], gamma)
      for (s in seq_len(H)) if (plan[i, s] > thresh)
        m <- c(m, fwd_marginal(posts[[i]], bar1, numeric(0), plan[i, ], s, dg, gamma, epsilon))
    }
    if (length(m) < 2) return(0)
    diff(range(m))
  }

  run_case <- function(seed, n, b, epsilon, gamma = 1, T_rounds = 2) {
    set.seed(seed)
    E_R <- 1 * 2 / (2 - 1); B_rem <- 2 * b * n * E_R      # full 2-round budget
    pop <- draw_initial_population(n, 1, 2, 1, 2, 0)
    p0  <- rpois(n, pmax(lambda_rate(pop$K0, pop$R0, gamma), 1e-12))
    sig <- draw_signals(pop$K0, pop$R0, 1, 1)
    posts <- build_posteriors_hist(p0, sig$sigma_r, sig$sigma_k, 400, 1, 2, 1, 2,
                                   gamma, 1, 1, TRUE, FALSE)
    H <- T_rounds
    t_s <- system.time(psm <- plan_forward_smooth(posts, B_rem, gamma, epsilon, H))[3]
    pgr <- plan_forward_ce(posts, B_rem, B_rem / 2 / 800, gamma, epsilon, H)  # fine greedy ns=800/round
    o_s <- plan_obj(posts, psm, gamma, epsilon)
    o_g <- plan_obj(posts, pgr, gamma, epsilon)
    cat(sprintf("b=%.2f eps=%.1f | smooth obj=%.6f spend=%.5f (t=%.3fs) | greedy800 obj=%.6f | smooth-greedy=%+.2e | ||dg||1/B=%.4f | KKT spread=%.2e\n",
                b, epsilon, o_s, sum(psm), t_s, o_g, o_s - o_g,
                sum(abs(psm - pgr)) / B_rem, kkt_spread(posts, psm, gamma, epsilon)))
    invisible(list(smooth = psm, greedy = pgr))
  }

  cat("=== Phase 2 forward-planner validation (H=2) ===\n")
  for (b in c(0.2, 0.3, 0.5)) run_case(7, 50, b, 0.0)
  for (b in c(0.2, 0.3, 0.5)) run_case(7, 50, b, 0.1)
}
