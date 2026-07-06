# ============================================================
# Phase 1 — exact single-round WATER-FILLING allocator (VALIDATED prototype)
# ============================================================
# Drop-in replacement for greedy_round_t() in the MYOPIC path (S4/S5/S6) and
# for the inner single-round solve. Removes the greedy step delta entirely:
# it computes the EXACT continuous KKT optimum of the round-t objective.
#
# STATUS: mathematically correct + fast. Validated 2026-07-05 (see harness at
# bottom, and docs/SMOOTH_ALLOCATOR_HANDOFF.md). NOT yet integrated into model.R.
#
# Objective (must match greedy_round_t exactly): maximize the POSTERIOR-EXPECTED
# round-t rate, summed over researchers, s.t. sum(g)=budget, g>=0:
#   V(g) = sum_i sum_m w_im * gamma * K_im * (R_im + g_i) / (K_im + R_im + g_i)
# where K_im is atom m's knowledge COMPOUNDED through g_hist. This is an
# expectation over the M-atom posterior, NOT lambda at the posterior mean.
#
# KKT: every FUNDED researcher shares a common marginal (water level) nu:
#   m_i(g_i) = sum_m w_im * gamma * K_im^2 / (K_im + R_im + g_i)^2 = nu
# Solve nu by an outer safeguarded-Newton bisection s.t. sum_i g_i(nu)=budget;
# each g_i(nu) by a warm-started inner Newton. Fully vectorized over the n x M
# atom matrix (base R). ~0.023 s/call at n=50, M=400 (2.5x FASTER than the
# CORRECT ns=800 greedy, and exact).
#
# Signature mirrors greedy_round_t (minus `delta`); init_g = seed floor (S6).

waterfill_round_t <- function(posts, budget, gamma, epsilon, g_hist = list(),
                              init_g = NULL, tol = 1e-9) {
  n <- length(posts); M <- length(posts[[1]]$w)
  g0 <- if (is.null(init_g)) rep(0, n) else init_g
  if (budget <= tol) return(g0)
  # Build n x M atom matrices; compound K through executed history g_hist.
  Amat <- matrix(0, n, M); Rmat <- matrix(0, n, M); Wmat <- matrix(0, n, M)
  for (i in seq_len(n)) {
    K <- posts[[i]]$K0
    for (gh in g_hist) K <- update_knowledge(K, posts[[i]]$R0 + gh[i], epsilon)
    Amat[i, ] <- K; Rmat[i, ] <- posts[[i]]$R0; Wmat[i, ] <- posts[[i]]$w
  }
  base <- Amat + Rmat + g0                             # denominator floor at g=0
  num  <- gamma * Wmat * Amat^2                        # marginal numerator per atom
  m0   <- rowSums(num / base^2)                        # marginal at g=0 (length n)
  gcur <- rep(0, n)                                    # warm-start carried across nu
  # S(nu) = total extra grant demanded at price nu; also returns dS/dnu.
  eval_nu <- function(nu) {
    active <- m0 > nu; g <- gcur; g[!active] <- 0
    for (k in 1:6) {                                   # inner Newton on m_i(g_i)=nu
      d <- base + g; m <- rowSums(num / d^2); mp <- rowSums(-2 * num / d^3)
      step <- (m - nu) / mp; step[!active] <- 0
      g <- pmax(0, g - step)
    }
    g[!active] <- 0; gcur <<- g
    d <- base + g; mp <- rowSums(-2 * num / d^3)
    list(g = g, S = sum(g), Sp = sum(ifelse(active, 1 / mp, 0)))
  }
  lo <- max(m0) * 1e-13; hi <- max(m0)                 # S decreasing in nu; bracket
  nu <- hi / 2
  for (it in 1:40) {                                   # outer safeguarded Newton
    e <- eval_nu(nu); diff <- e$S - budget
    if (abs(diff) < tol * budget) break
    if (diff > 0) lo <- nu else hi <- nu               # keep bracket [lo,hi]
    nu_n <- nu - diff / e$Sp                            # Newton step
    nu <- if (is.finite(nu_n) && nu_n > lo && nu_n < hi) nu_n else (lo + hi) / 2
  }
  g0 + eval_nu(nu)$g
}

# ------------------------------------------------------------
# VALIDATION HARNESS (run: Rscript tests/waterfill_round_t.R)
# Confirms: exact (>= fine-greedy ns=3200, spends budget), fast.
# ------------------------------------------------------------
if (sys.nframe() == 0L) {
  setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
  source("model.R")
  obj <- function(posts, g) sum(vapply(seq_along(posts), function(i)
    post_lambda_round_t(posts[[i]], list(), g[i], 1, 0), numeric(1)))
  set.seed(7); n <- 50; B <- 0.3 * n * 2
  pop <- draw_initial_population(n, 1, 2, 1, 2, 0)
  p0  <- rpois(n, pmax(lambda_rate(pop$K0, pop$R0, 1), 1e-12))
  sig <- draw_signals(pop$K0, pop$R0, 1, 1)
  posts <- build_posteriors_hist(p0, sig$sigma_r, sig$sigma_k, 400, 1, 2, 1, 2, 1, 1, 1, TRUE, FALSE)
  tw <- system.time(for (r in 1:50) gw <- waterfill_round_t(posts, B, 1, 0))[3] / 50
  g3200 <- greedy_round_t(posts, B, B / 3200, 1, 0)
  cat(sprintf("water-fill: obj=%.7f spend=%.7f time=%.5fs/call\n", obj(posts, gw), sum(gw), tw))
  cat(sprintf("vs fine greedy ns=3200: obj=%.7f  wf-greedy=%+.2e (>=0 expected)  ||dg||1/B=%.6f\n",
      obj(posts, g3200), obj(posts, gw) - obj(posts, g3200), sum(abs(gw - g3200)) / B))
  cat(sprintf("water-fill = %.1fx greedy-ns50 (0.00335s), %.1fx greedy-ns800 (0.063s)\n",
      tw / 0.00335, tw / 0.063))
}
