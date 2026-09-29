# Referee point 1: how far is the two-round complete-information optimum from the
# gap rule applied each round with the budget split evenly? Deterministic growth
# K' = K + eps*lambda, so the T = 2 problem is a smooth deterministic program.
set.seed(1)
lam <- function(K, R) K * R / (K + R)
gap_rule <- function(K, R, B) {
  if (B <= 1e-12) return(rep(0, length(K)))
  f <- function(cc) sum(pmax(cc * K - R, 0)) - B
  cc <- uniroot(f, c(1e-9, 1e7), tol = 1e-12)$root
  pmax(cc * K - R, 0)
}
rpareto <- function(n, xm, a) xm / runif(n)^(1 / a)
one <- function(seed, n = 50, aK = 2, kmin = 1, eps = 0.1, b = 1) {
  set.seed(seed)
  K <- rpareto(n, kmin, aK); R <- rpareto(n, 1, 2)
  B <- b * n * 2                       # b n E[R0], E[R0] = 2 at alpha_R = 2, scale 1
  total <- function(g1, B2) {
    K2 <- K + eps * lam(K, R + g1)
    g2 <- gap_rule(K2, R, B2)
    sum(lam(K, R + g1)) + sum(lam(K2, R + g2))
  }
  # benchmark: gap rule each round, half the budget each round (the paper's oracle)
  g1b <- gap_rule(K, R, B / 2); bench <- total(g1b, B / 2)
  none <- 2 * sum(lam(K, R)) + eps * sum(lam(K, R)) * 0  # no funding: K grows too
  K2n <- K + eps * lam(K, R); none <- sum(lam(K, R)) + sum(lam(K2n, R))
  # dynamic optimum: g1 = B * s * softmax(theta), s in (0,1)
  obj <- function(par) {
    s <- plogis(par[1]); w <- exp(par[-1] - max(par[-1])); w <- w / sum(w)
    g1 <- B * s * w
    -total(g1, B * (1 - s))
  }
  init <- c(0, log(g1b + 1e-3))
  o <- optim(init, obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-9))
  s <- plogis(o$par[1]); w <- exp(o$par[-1] - max(o$par[-1])); w <- w / sum(w); g1d <- B * s * w
  dyn <- -o$value
  c(seed = seed, aK = aK, eps = eps, b = b,
    gain_bench = bench - none, gain_dyn = dyn - none,
    loss_share = (dyn - bench) / (dyn - none),     # share of funding's gain lost by the static rule
    share_r1 = s,                                  # dynamic optimum's round-1 spending share
    cor_g1 = cor(g1d, g1b), n_funded_bench = sum(g1b > 1e-6), n_funded_dyn = sum(g1d > 1e-3),
    top_bench = max(g1b) / sum(g1b), top_dyn = max(g1d) / sum(g1d))
}
grid <- expand.grid(aK = c(1.3, 3.5), eps = c(0.1, 0.85), b = c(0.2, 1))
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  kmin <- 2 * (g$aK - 1) / g$aK
  rows <- t(sapply(1:6, function(s) one(s, aK = g$aK, kmin = kmin, eps = g$eps, b = g$b)))
  colMeans(rows)
}))
res <- as.data.frame(res)
print(round(res[, c("aK", "eps", "b", "loss_share", "share_r1", "cor_g1", "n_funded_bench", "n_funded_dyn", "top_bench", "top_dyn")], 3))
write.csv(res, "dynopt_results.csv", row.names = FALSE)
