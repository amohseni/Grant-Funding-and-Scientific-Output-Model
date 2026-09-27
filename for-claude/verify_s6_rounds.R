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
# Sharp-end reversal probe: decompose S5's value at tau 0.05 vs 0.3 (k=1.3, T=2)
# H1 exploration/dynamics: tau=0.05 wins round 1, loses round 2.
# H2 plug-in/CE artifact: tau=0.05 already loses round 1.
suppressPackageStartupMessages(library(parallel))
source("model.R")
probe <- function(tk, sd, ks = 1.3) {
  rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks,
                         allocator = "smooth", strategies = c(4, 5), detail = TRUE)
  s5 <- rp$strategies[[5]]; s4 <- rp$strategies[[4]]
  r1_5 <- sum(s5$lam_full[[1]]); r2_5 <- sum(s5$lam_full[[2]])
  r1_4 <- sum(s4$lam_full[[1]]); r2_4 <- sum(s4$lam_full[[2]])
  gini <- function(x) { x <- sort(x); n <- length(x); if (sum(x)==0) return(NA); sum((2*seq_len(n)-n-1)*x)/(n*sum(x)) }
  c(v1 = r1_5 - r1_4, v2 = r2_5 - r2_4, tot = s5$total_expected - s4$total_expected,
    gini1_5 = gini(s5$g_rounds[[1]]))
}
cores <- max(1L, detectCores() - 1L)
a <- do.call(rbind, mclapply(1:50, function(sd) probe(0.05, sd), mc.cores = cores))
b <- do.call(rbind, mclapply(1:50, function(sd) probe(0.3, sd), mc.cores = cores))
d <- a - b
for (col in c("v1","v2","tot")) {
  x <- d[,col]
  cat(sprintf("%s: mean diff (.05-.3) = %.3f  se = %.3f  z = %.2f\n", col, mean(x), sd(x)/sqrt(50), mean(x)/(sd(x)/sqrt(50))))
}
cat(sprintf("round-1 grant gini: tau=.05 %.4f  tau=.3 %.4f  paired z = %.2f\n",
    mean(a[,"gini1_5"]), mean(b[,"gini1_5"]),
    mean(a[,"gini1_5"]-b[,"gini1_5"])/(sd(a[,"gini1_5"]-b[,"gini1_5"])/sqrt(50))))
suppressPackageStartupMessages(library(parallel))
source("model.R")
probe <- function(tk, sd, ks = 1.3) {
  rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, k_shape = ks,
                         allocator = "smooth", strategies = c(4, 5), detail = TRUE)
  out <- c()
  for (S in c(4,5)) {
    st <- rp$strategies[[S]]
    r <- sapply(1:2, function(t) sum(lambda_rate(st$K_rounds[[t]], rp$R0_at_start + st$g_rounds[[t]], 1)))
    out <- c(out, r, sum(r) - st$total_expected)  # identity check
  }
  names(out) <- c("r1_4","r2_4","chk4","r1_5","r2_5","chk5")
  out
}
cores <- max(1L, detectCores() - 1L)
a <- do.call(rbind, mclapply(1:50, function(sd) probe(0.05, sd), mc.cores = cores))
b <- do.call(rbind, mclapply(1:50, function(sd) probe(0.3, sd), mc.cores = cores))
cat(sprintf("identity check max |sum(rounds)-total|: %.2e\n", max(abs(c(a[,c("chk4","chk5")], b[,c("chk4","chk5")])))))
v1 <- (a[,"r1_5"]-a[,"r1_4"]) - (b[,"r1_5"]-b[,"r1_4"])
v2 <- (a[,"r2_5"]-a[,"r2_4"]) - (b[,"r2_5"]-b[,"r2_4"])
cat(sprintf("round-1 value diff (.05-.3): %.3f se %.3f z %.2f\n", mean(v1), sd(v1)/sqrt(50), mean(v1)/(sd(v1)/sqrt(50))))
cat(sprintf("round-2 value diff (.05-.3): %.3f se %.3f z %.2f\n", mean(v2), sd(v2)/sqrt(50), mean(v2)/(sd(v2)/sqrt(50))))
suppressPackageStartupMessages(library(parallel))
source("model.R")
val <- function(tk, sd, tr, ks = 1.3) {
  rp <- run_simulation_T(seed = sd, T_rounds = 2, tau_k = tk, tau_r = tr, k_shape = ks,
                         allocator = "smooth", strategies = c(4, 5))
  rp$strategies[[5]]$total_expected - rp$strategies[[4]]$total_expected
}
cores <- max(1L, detectCores() - 1L)
for (tr in c(1, 0.01)) {
  d <- simplify2array(mclapply(1:50, function(sd) val(0.05, sd, tr) - val(0.3, sd, tr), mc.cores = cores))
  cat(sprintf("tau_r=%.2f: value diff (.05-.3) = %.3f  se = %.3f  z = %.2f\n",
              tr, mean(d), sd(d)/sqrt(50), mean(d)/(sd(d)/sqrt(50))))
}
