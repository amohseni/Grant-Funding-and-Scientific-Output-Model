# ============================================================
# tests/ce_tax_vs_T.R  —  Is the post-T=3 PG decline a CE artifact?
# ------------------------------------------------------------
# Decomposes the certainty-equivalence (CE) valuation tax along the
# horizon. At each planning round t of the CE-forward policy we compare
#   V_CE_t : CE's point-estimate (p̄) lookahead valuation of the remaining
#            horizon, and
#   V_SP_t : a scenario-integrated (stochastic-programming) valuation of
#            the SAME round-t decision — sampling p_t and re-planning the
#            continuation, i.e. the proper E over realizations.
# gap_t = V_SP_t − V_CE_t is the per-round CE Jensen tax. Summed over
# rounds it estimates the total CE tax at horizon T.
#
# Logic: if Σ gap_t stays tiny as T grows, CE is faithful → the PG decline
# is REAL. If Σ gap_t grows with T and is comparable to the PG deficit,
# the decline is (partly) a CE artifact.
# ============================================================
setwd("/Users/amohseni/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
suppressPackageStartupMessages(library(parallel))
src <- readLines("app.R", warn=FALSE); end <- grep("^shinyApp\\(", src)[1]-1L
eval(parse(text=paste(src[1:end],collapse="\n")), envir=globalenv())
source("simulate_T.R")

# CE self-valuation of a forward plan from posteriors `posts`, given executed
# history g_hist (list of past grant vectors). Mirrors the greedy's objective.
ce_plan_value <- function(posts, B_rem, delta, gamma, epsilon, H, g_hist) {
  plan <- plan_forward_ce(posts, B_rem, delta, gamma, epsilon, H, g_hist = g_hist)
  n <- length(posts)
  val <- sum(vapply(seq_len(n), function(i) {
    ghi <- vapply(g_hist, function(gh) gh[i], numeric(1))
    pb  <- ce_reweight_posterior(posts[[i]], plan[i, 1], gamma)
    fwd_researcher_value(posts[[i]], pb, ghi, plan[i, ], gamma, epsilon)
  }, numeric(1)))
  list(plan = plan, value = val)
}

# CE continuation value for rounds (t+1..T) from a posterior list, given the
# executed history so far. H_cont = T - t; returns 0 when terminal.
ce_continue_value <- function(posts, B_rem, delta, gamma, epsilon, H_cont, g_hist) {
  if (H_cont < 1) return(0)
  ce_plan_value(posts, B_rem, delta, gamma, epsilon, H_cont, g_hist)$value
}

# One trajectory of the CE-forward policy, instrumented to return per-round
# CE tax gap_t and total forward value.
ce_tax_trajectory <- function(seed, T_rounds, n, epsilon, b, tau_k, k_shape,
                              M = 200, n_steps = 50, S = 20, gamma = 1) {
  set.seed(seed)
  E_R <- 2; B <- b * n * E_R; B_total <- 2 * B; delta <- B / n_steps
  pop <- draw_initial_population(n, 1, k_shape, 1, 2, 0); K <- pop$K0; R0 <- pop$R0
  p0  <- rpois(n, pmax(lambda_rate(K, R0, gamma), 1e-12))
  sig <- draw_signals(K, R0, 1, tau_k)

  K_cur <- K; B_rem <- B_total
  p_hist <- list(); g_hist <- list()
  gaps <- numeric(0); Vfwd_self <- NA_real_; realized <- 0

  for (t in seq_len(T_rounds)) {
    posts <- build_posteriors_hist(p0, sig$sigma_r, sig$sigma_k, M, 1, k_shape, 1, 2,
                                   gamma, 1, tau_k, TRUE, TRUE, p_hist = p_hist, g_hist = g_hist)
    H <- T_rounds - t + 1L
    pv <- ce_plan_value(posts, B_rem, delta, gamma, epsilon, H, g_hist)
    g_t <- pv$plan[, 1]
    if (t == 1L) Vfwd_self <- pv$value

    # immediate round-t expected value under posterior
    v_imm <- sum(vapply(seq_len(n), function(i) {
      ghi <- vapply(g_hist, function(gh) gh[i], numeric(1))
      post_lambda_round_t(posts[[i]], ghi, g_t[i], gamma, epsilon)
    }, numeric(1)))
    CEcont_pt <- pv$value - v_imm                       # CE point-estimate continuation

    # scenario-integrated continuation: sample p_t, update posterior, CE-continue
    if (H >= 2) {
      B_rem_c <- B_rem - sum(g_t); gh2 <- c(g_hist, list(g_t))
      SPcont <- mean(vapply(seq_len(S), function(s) {
        p_t <- vapply(seq_len(n), function(i) {
          j <- sample.int(M, 1, prob = posts[[i]]$w)
          Kc <- posts[[i]]$K0[j]
          for (gh in g_hist) Kc <- update_knowledge(Kc, posts[[i]]$R0[j] + gh[i], epsilon)
          rpois(1, pmax(lambda_rate(Kc, posts[[i]]$R0[j] + g_t[i], gamma), 1e-12))
        }, numeric(1))
        postt <- build_posteriors_hist(p0, sig$sigma_r, sig$sigma_k, M, 1, k_shape, 1, 2,
                                       gamma, 1, tau_k, TRUE, TRUE,
                                       p_hist = c(p_hist, list(p_t)), g_hist = gh2)
        ce_continue_value(postt, B_rem_c, delta, gamma, epsilon, T_rounds - t, gh2)
      }, numeric(1)))
      gaps <- c(gaps, SPcont - CEcont_pt)              # gap_t (round t's Jensen tax)
    }

    # advance the ACTUAL policy under truth
    R_t <- R0 + g_t; lam_t <- lambda_rate(K_cur, R_t, gamma)
    realized <- realized + sum(lam_t)
    p_hist[[t]] <- rpois(n, pmax(lam_t, 1e-12)); g_hist[[t]] <- g_t
    B_rem <- B_rem - sum(g_t); K_cur <- update_knowledge(K_cur, R_t, epsilon)
  }
  c(sum_gap = sum(gaps), Vfwd_self = Vfwd_self, realized = realized,
    n_gaps = length(gaps), mean_gap = if (length(gaps)) mean(gaps) else NA_real_)
}

# ---------------- run ----------------
N <- 40; EPS <- 0.30; TAUK <- 1.0; KSH <- 2.0; BB <- 0.5
seeds <- 1:16; S_SCEN <- 20; NCORES <- max(1, detectCores() - 1)

cat(sprintf("CE-tax vs T  (n=%d, eps=%.2f, tau_k=%.1f, k_shape=%.1f, b=%.2f; %d seeds, %d scenarios)\n\n",
            N, EPS, TAUK, KSH, BB, length(seeds), S_SCEN))
cat(sprintf("%3s | %11s %9s | %11s %8s | %s\n",
            "T", "PG(S8-S5)", "PGdefic", "CE_tax_tot", "%ofFwd", "interpretation"))

for (Tr in 2:5) {
  # CE tax along trajectory
  mats <- mclapply(seeds, function(sd) ce_tax_trajectory(sd, Tr, N, EPS, BB, TAUK, KSH, S = S_SCEN),
                   mc.cores = NCORES)
  M <- do.call(rbind, mats)
  tax   <- mean(M[, "sum_gap"]);  tax_se <- sd(M[, "sum_gap"]) / sqrt(nrow(M))
  fwdval<- mean(M[, "realized"]); taxpct <- tax / fwdval * 100

  # realized PG at this exact cell (more seeds, cheap)
  pgv <- vapply(1:120, function(sd) {
    r <- run_simulation_T(seed = sd, T_rounds = Tr, n = N, epsilon = EPS, b = BB,
                          tau_k = TAUK, k_shape = KSH, x_seed = 0.25, M = 200, strategies = c(5, 8))
    r$strategies[[8]]$total_expected - r$strategies[[5]]$total_expected
  }, numeric(1))
  PG <- mean(pgv); PGse <- sd(pgv) / sqrt(length(pgv))

  interp <- if (Tr == 2) "(baseline)" else
            if (abs(tax) > 0.5 * abs(min(0, PG)) && PG < 0) "CE tax >~ deficit -> artifact plausible" else
            if (PG < 0) "CE tax << deficit -> decline looks REAL" else "PG>=0"
  cat(sprintf("%3d | %+6.3f±%.3f  %8.3f | %+8.3f±%.3f  %5.2f%% | %s\n",
              Tr, PG, PGse, min(0, PG), tax, tax_se, taxpct, interp))
}
cat("\nCE_tax_tot = Σ_t (scenario-integrated − CE point-estimate) continuation value.\n")
cat("If CE_tax grows with T and is comparable to |PG deficit|, the high-T decline is partly a CE artifact.\n")
