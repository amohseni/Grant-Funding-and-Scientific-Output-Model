# C5a: does the records-only funder's allocation get closer to the optimal
# allocation across rounds? T=8, correlation of S4 round-t grants with the
# gap rule computed at the round-t state (true K entering round t, baseline R0).
suppressPackageStartupMessages(library(parallel))
source("model.R")
gap_rule <- function(K, R, budget) {
  f <- function(cc) sum(pmax(cc * K - R, 0)) - budget
  cc <- uniroot(f, c(1e-9, 1e6), tol = 1e-10)$root
  pmax(cc * K - R, 0)
}
Tn <- 8
run1 <- function(ks, sd) {
  rp <- run_simulation_T(seed = sd, T_rounds = Tn, tau_k = 1, k_shape = ks,
                         allocator = "smooth", strategies = c(4), detail = TRUE)
  tranche <- rp$params$B_total / Tn
  s4 <- rp$strategies[[4]]
  sapply(1:Tn, function(t) {
    K_t <- s4$K_rounds[[t]]; R_t <- rp$R0_at_start
    suppressWarnings(cor(s4$g_rounds[[t]], gap_rule(K_t, R_t, tranche)))
  })
}
for (ks in c(1.3, 2)) {
  m <- do.call(rbind, mclapply(1:50, function(sd) run1(ks, sd), mc.cores = max(1L, detectCores()-1L)))
  cat(sprintf("k=%.1f  mean corr by round: %s\n  se: %s\n", ks,
      paste(sprintf("%.3f", colMeans(m, na.rm=TRUE)), collapse=" "),
      paste(sprintf("%.3f", apply(m, 2, sd, na.rm=TRUE)/sqrt(nrow(m))), collapse=" ")))
}
# C5a long horizon: records-only (S4) vs records+review (S5) allocation
# correlation with the round-t gap rule, T=20, tau_k=1.
suppressPackageStartupMessages(library(parallel))
source("model.R")
gap_rule <- function(K, R, budget) {
  f <- function(cc) sum(pmax(cc * K - R, 0)) - budget
  cc <- uniroot(f, c(1e-9, 1e6), tol = 1e-10)$root
  pmax(cc * K - R, 0)
}
Tn <- 20
run1 <- function(ks, sd) {
  rp <- run_simulation_T(seed = sd, T_rounds = Tn, tau_k = 1, k_shape = ks,
                         allocator = "smooth", strategies = c(4, 5), detail = TRUE)
  tranche <- rp$params$B_total / Tn
  out <- c()
  for (S in c(4, 5)) {
    st <- rp$strategies[[S]]
    out <- c(out, sapply(1:Tn, function(t)
      suppressWarnings(cor(st$g_rounds[[t]], gap_rule(st$K_rounds[[t]], rp$R0_at_start, tranche)))))
  }
  out
}
for (ks in c(1.3, 2)) {
  m <- do.call(rbind, mclapply(1:50, function(sd) run1(ks, sd), mc.cores = max(1L, detectCores()-1L)))
  mu <- colMeans(m, na.rm = TRUE)
  cat(sprintf("k=%.1f S4 rounds 1,2,4,8,12,16,20: %s\n", ks,
      paste(sprintf("%.3f", mu[c(1,2,4,8,12,16,20)]), collapse=" ")))
  cat(sprintf("k=%.1f S5 rounds 1,2,4,8,12,16,20: %s\n", ks,
      paste(sprintf("%.3f", mu[Tn + c(1,2,4,8,12,16,20)]), collapse=" ")))
}
# Paired per-seed test: is review's value lower at tau=0.05 than tau=0.3 (heavy tails)?
suppressPackageStartupMessages(library(parallel))
source("model.R")
val <- function(tk, ks, sd) {
  rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks,
                         allocator = "smooth", strategies = c(4, 5), detail = TRUE)
  rp$strategies[[5]]$total_expected - rp$strategies[[4]]$total_expected
}
for (ks in c(1.3, 2)) {
  d <- simplify2array(mclapply(1:50, function(sd) val(0.05, ks, sd) - val(0.3, ks, sd),
                               mc.cores = max(1L, detectCores() - 1L)))
  cat(sprintf("k=%.1f: mean diff (tau .05 - .3) = %.3f, se = %.3f, z = %.2f\n",
              ks, mean(d), sd(d)/sqrt(50), mean(d)/(sd(d)/sqrt(50))))
}
