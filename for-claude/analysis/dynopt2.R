# Package A: (1) two-round complete-information dynamic optimum vs the gap rule applied each
# round with an even split, full grid, 20 populations; (2) T = 5 check restricted to the
# spending shares (gap-rule recipients each round). Fixed mean capability 2, n = 50.
lam <- function(K, R) K * R / (K + R)
gap_rule <- function(K, R, B) {
  if (B <= 1e-12) return(rep(0, length(K)))
  f <- function(cc) sum(pmax(cc * K - R, 0)) - B
  cc <- uniroot(f, c(1e-9, 1e7), tol = 1e-12)$root
  pmax(cc * K - R, 0)
}
rpareto <- function(n, xm, a) xm / runif(n)^(1 / a)
draw <- function(seed, n, aK) { set.seed(seed); list(K = rpareto(n, 2 * (aK - 1) / aK, aK), R = rpareto(n, 1, 2)) }
two_round <- function(seed, aK, eps, b, n = 50) {
  p <- draw(seed, n, aK); K <- p$K; R <- p$R; B <- b * n * 2
  total <- function(g1, B2) { K2 <- K + eps * lam(K, R + g1); g2 <- gap_rule(K2, R, B2); sum(lam(K, R + g1)) + sum(lam(K2, R + g2)) }
  g1b <- gap_rule(K, R, B / 2); bench <- total(g1b, B / 2)
  K2n <- K + eps * lam(K, R); none <- sum(lam(K, R)) + sum(lam(K2n, R))
  obj <- function(par) { s <- plogis(par[1]); w <- exp(par[-1] - max(par[-1])); w <- w / sum(w); -total(B * s * w, B * (1 - s)) }
  o <- optim(c(0, log(g1b + 1e-3)), obj, method = "BFGS", control = list(maxit = 2000, reltol = 1e-9))
  s <- plogis(o$par[1]); w <- exp(o$par[-1] - max(o$par[-1])); w <- w / sum(w); g1d <- B * s * w; dyn <- -o$value
  c(aK = aK, eps = eps, b = b, loss_share = (dyn - bench) / (dyn - none), share_r1 = s, cor_g1 = cor(g1d, g1b))
}
five_round <- function(seed, aK, eps, b, n = 50, T = 5) {
  p <- draw(seed, n, aK); K0 <- p$K; R <- p$R; B <- b * n * 2
  run <- function(shares) { K <- K0; tot <- 0; for (t in 1:T) { g <- gap_rule(K, R, B * shares[t]); tot <- tot + sum(lam(K, R + g)); K <- K + eps * lam(K, R + g) }; tot }
  none <- run(rep(0, T)); even <- run(rep(1 / T, T))
  obj <- function(th) { w <- exp(c(0, th) - max(c(0, th))); -run(w / sum(w)) }
  o <- optim(rep(0, T - 1), obj, method = "BFGS", control = list(reltol = 1e-10))
  w <- exp(c(0, o$par) - max(c(0, o$par))); w <- w / sum(w)
  c(aK = aK, eps = eps, b = b, loss_share_even = (-o$value - even) / (-o$value - none), com = sum(seq_len(T) * w) / (T + 1), share_last = w[T])
}
grid <- expand.grid(aK = c(1.3, 2, 3.5), eps = c(0.1, 0.5, 0.85), b = c(0.2, 1, 2))
t0 <- Sys.time()
r5 <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) { g <- grid[i, ]; colMeans(t(sapply(1:20, function(s) five_round(s, g$aK, g$eps, g$b)))) }))
write.csv(as.data.frame(r5), "dynopt5_results.csv", row.names = FALSE)
cat("T=5 done", round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n")
r2 <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; rows <- t(sapply(1:20, function(s) two_round(s, g$aK, g$eps, g$b)))
  out <- colMeans(rows); out["loss_share_max"] <- max(rows[, "loss_share"]); out["cor_g1_min"] <- min(rows[, "cor_g1"])
  cat(i, round(as.numeric(Sys.time() - t0, units = "mins"), 1), "min\n"); out
}))
write.csv(as.data.frame(r2), "dynopt2_results.csv", row.names = FALSE)
cat("DYNOPT DONE\n")
