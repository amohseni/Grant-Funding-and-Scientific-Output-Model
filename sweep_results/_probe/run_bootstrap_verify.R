# Bootstrap verification items (addendum end-matter / notes-bootstrap-integration §4):
# 1. Honest fixed-schedule sweep at the depth corner (no CE planner).
# 2. CE self-consistency at depth (honest mispricing vs optimizer bug).
# 3. 200-seed re-runs of the exploration sweeps (main-text provenance).
suppressPackageStartupMessages(library(parallel))
setwd("/Users/amohseni/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model")
source("model.R"); source("sweep.R"); source("sweep_T.R")
t0 <- Sys.time()
cores <- max(1L, detectCores() - 1L)
outdir <- "sweep_results/bootstrap_verify"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Item 3 first (manifest machinery, checkpointed) ----------
for (nm in c("exploration_corner", "exploration_poverty", "exploration_depth")) {
  sweep_one_T(nm, seeds = 1:200,
              base_params = list(allocator = "smooth", strategies = c(2, 5, 7, 8), M = 200),
              out_dir = "sweep_results/exploration_200", resume = TRUE)
  cat(sprintf("[item3 %s done] %.0fs\n", nm, as.numeric(Sys.time()-t0, units="secs")))
}
writeLines(c(paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
             paste("date:", format(Sys.time())), "seeds: 1:200"),
           "sweep_results/exploration_200/RUN_INFO.txt")

# ---------- Item 1: honest fixed two-block schedules, myopic within-round ----------
# Depth corner: b=3, r_min=0.001, tau_k=100, k_shape=1.3, T=6, budget_ref=K.
# Schedule (x, k): share x of B_total spread evenly over rounds 1..k, remainder
# evenly over k+1..T. Within each round: waterfill on the current posterior
# (honest re-deciding funder, TRUE dynamics, no CE lookahead).
honest_run <- function(seed, sched, Tr = 6, n = 50, b = 3, r_min = 0.001,
                       tau_k = 100, k_shape = 1.3, eps = 1e-4, M = 200, gamma = 1) {
  set.seed(seed)
  E_K <- 1 * k_shape / (k_shape - 1)
  B_total <- 2 * b * n * E_K
  pop <- draw_initial_population(n, 1, k_shape, r_min, 2, 0)
  K <- pop$K0; R0 <- pop$R0
  lam0 <- lambda_rate(K, R0, gamma)
  p_init <- rpois(n, pmax(lam0, LAMBDA_FLOOR))
  sigs <- draw_signals(K, R0, 1, tau_k)
  set.seed(seed * 1000)
  budget_t <- sched * B_total
  p_hist <- list(); g_hist <- list(); tot <- 0
  for (t in seq_len(Tr)) {
    posts <- build_posteriors_hist(p_init, sigs$sigma_r, sigs$sigma_k, M,
                                   1, k_shape, r_min, 2, gamma, 1, tau_k,
                                   TRUE, TRUE, p_hist = p_hist, g_hist = g_hist)
    g_t <- waterfill_round_t(posts, budget_t[t], gamma, eps, g_hist = g_hist)
    lam_t <- lambda_rate(K, R0 + g_t, gamma)
    p_t <- rpois(n, pmax(lam_t, LAMBDA_FLOOR))
    tot <- tot + sum(lam_t)
    p_hist[[t]] <- p_t; g_hist[[t]] <- g_t
    K <- update_knowledge_split(K, R0, g_t, eps)
  }
  tot
}
two_block <- function(x, k, Tr = 6) c(rep(x / k, k), rep((1 - x) / (Tr - k), Tr - k))
scheds <- list(even = rep(1/6, 6))
for (x in c(0.1, 0.3, 0.5, 0.7, 0.9)) for (k in c(1, 2, 3))
  scheds[[sprintf("x%.1f_k%d", x, k)]] <- two_block(x, k)
seeds <- 1:200
for (eps_cell in c(1e-4, 0.3)) {
  tab <- do.call(rbind, lapply(names(scheds), function(nm) {
    v <- unlist(mclapply(seeds, honest_run, sched = scheds[[nm]], eps = eps_cell, mc.cores = cores))
    data.frame(schedule = nm, eps = eps_cell, out = mean(v), se = sd(v)/sqrt(length(v)))
  }))
  tab <- tab[order(-tab$out), ]
  fn <- sprintf("honest_schedules_eps%s.csv", ifelse(eps_cell < 0.01, "0", "03"))
  write.csv(tab, file.path(outdir, fn), row.names = FALSE)
  cat(sprintf("== Item 1 (eps=%g): top schedules ==\n", eps_cell))
  print(head(tab, 6), digits = 5, row.names = FALSE)
  ev <- tab$out[tab$schedule == "even"]
  early <- tab[grepl("k1|k2", tab$schedule) & as.numeric(sub("x([0-9.]+)_.*", "\\1", tab$schedule)) > 0.34, ]
  cat(sprintf("even = %.2f | best early-mass = %s (%.2f) | early beats even by %.3f (z=%.1f)\n",
      ev, early$schedule[which.max(early$out)], max(early$out),
      max(early$out) - ev, (max(early$out) - ev) / sqrt(tab$se[tab$schedule=="even"]^2 + early$se[which.max(early$out)]^2)))
  cat(sprintf("[item1 eps=%g done] %.0fs\n", eps_cell, as.numeric(Sys.time()-t0, units="secs")))
}

# ---------- Item 2: CE self-consistency at depth (20 seeds) ----------
ce_value_plan <- function(posts, plan, gamma, eps) {
  # CE objective value of a full plan (n x H) at t=1, matching the smooth planner
  sum(vapply(seq_along(posts), function(i) {
    p0 <- posts[[i]]
    pb <- ce_reweight_posterior(p0, plan[i, 1], gamma)
    fwd_researcher_value(p0, pb, numeric(0), plan[i, ], gamma, eps)
  }, numeric(1)))
}
cell <- function(seed, Tr = 6, n = 50, b = 3, r_min = 0.001, tau_k = 100,
                 k_shape = 1.3, eps = 1e-4, M = 200, gamma = 1) {
  set.seed(seed)
  E_K <- 1 * k_shape / (k_shape - 1); B_total <- 2 * b * n * E_K
  pop <- draw_initial_population(n, 1, k_shape, r_min, 2, 0)
  K <- pop$K0; R0 <- pop$R0
  p_init <- rpois(n, pmax(lambda_rate(K, R0, gamma), LAMBDA_FLOOR))
  sigs <- draw_signals(K, R0, 1, tau_k)
  set.seed(seed * 1000)
  posts <- build_posteriors_hist(p_init, sigs$sigma_r, sigs$sigma_k, M,
                                 1, k_shape, r_min, 2, gamma, 1, tau_k, TRUE, TRUE)
  plan_S8 <- plan_forward_smooth(posts, B_total, gamma, eps, Tr)
  plan_S5 <- matrix(0, n, Tr)   # even-tranche myopic, CE-anticipated round by round
  gh <- list()
  for (t in seq_len(Tr)) {
    plan_S5[, t] <- waterfill_round_t(posts, B_total / Tr, gamma, eps, g_hist = gh)
    gh[[t]] <- plan_S5[, t]
  }
  c(ce_S8 = ce_value_plan(posts, plan_S8, gamma, eps),
    ce_S5 = ce_value_plan(posts, plan_S5, gamma, eps))
}
cc <- do.call(rbind, mclapply(1:20, cell, mc.cores = cores))
d <- cc[, "ce_S8"] - cc[, "ce_S5"]
cat(sprintf("\n== Item 2: CE self-consistency (depth corner, 20 seeds) ==\nCE(S8 plan) - CE(S5-like plan): mean %.3f (se %.3f) | S8 >= S5 in %d/20 seeds -> %s\n",
    mean(d), sd(d)/sqrt(20), sum(d >= 0),
    ifelse(mean(d) >= 0, "HONEST MISPRICING (planner maximizes its own objective)", "OPTIMIZER BUG")))
write.csv(data.frame(seed = 1:20, ce_S8 = cc[, "ce_S8"], ce_S5 = cc[, "ce_S5"]),
          file.path(outdir, "ce_self_consistency.csv"), row.names = FALSE)
writeLines(c(paste("hash:", system("git rev-parse HEAD", intern = TRUE)),
             paste("date:", format(Sys.time()))), file.path(outdir, "RUN_INFO.txt"))
cat(sprintf("BOOTSTRAP VERIFY DONE in %.0f min\n", as.numeric(Sys.time()-t0, units="mins")))
