suppressPackageStartupMessages(source("model.R"))
set.seed(2)
res <- list()
for (aK in c(1.3, 2, 3.5)) for (tk in c(0.05, 0.3, 1, 3)) for (M in c(200, 2000)) {
  kmin <- 2*(aK-1)/aK
  ess <- c(); relerr_top <- c(); relerr_all <- c()
  for (rep in 1:20) {
    K <- rpareto(50, kmin, aK); R <- rpareto(50, 1, 2)
    p <- rpois(50, lambda_rate(K, R, 1)); sr <- R + rnorm(50, 0, 1); sk <- K + rnorm(50, 0, tk)
    for (i in 1:50) {
      po <- posterior_samples_single(p[i], sr[i], sk[i], M, kmin, aK, 1, 2, 1, 1, tk, TRUE, TRUE)
      ess <- c(ess, 1/sum(po$w^2))
      pm <- sum(po$w * po$K0); e <- (pm - K[i])/K[i]
      relerr_all <- c(relerr_all, e); if (K[i] > quantile(K, 0.9)) relerr_top <- c(relerr_top, e)
    }
  }
  res[[length(res)+1]] <- data.frame(aK=aK, tau_k=tk, M=M, ess_median=median(ess), ess_q10=quantile(ess,0.1), ess_min=min(ess),
     bias_top=mean(relerr_top), mad_top=mean(abs(relerr_top)), bias_all=mean(relerr_all))
}
print(do.call(rbind, res), digits=3, row.names=FALSE)
