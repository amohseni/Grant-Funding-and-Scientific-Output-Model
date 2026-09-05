"""Generates fig2-curve.dat: value of targeting vs budget scale, exact from Proposition 1.
Population: n=75, K,R0 ~ Pareto(alpha=2, min 1), rho=0, A=1, numpy default_rng(659)."""
import numpy as np
rng = np.random.default_rng(659)
n, A = 75, 1.0
K  = rng.pareto(2.0, n) + 1.0
R0 = rng.pareto(2.0, n) + 1.0
def phi(K, R): return 2*A*K*R/(K+R)
def gap_alloc(B):
    lo, hi = 0.0, 1.0
    G = lambda c: np.maximum(c*K - R0, 0).sum()
    while G(hi) < B: hi *= 2
    for _ in range(200):
        mid = (lo+hi)/2
        lo, hi = (mid, hi) if G(mid) < B else (lo, mid)
    return np.maximum(((lo+hi)/2)*K - R0, 0)
f0 = phi(K, R0).sum(); S = R0.sum()
with open("fig2-curve.dat", "w") as f:
    f.write("b v\n")
    for b in np.concatenate([np.arange(0.01, 0.1, 0.01), np.arange(0.1, 2.001, 0.02)]):
        B = b*S
        fu = phi(K, R0 + B/n).sum()
        fo = phi(K, R0 + gap_alloc(B)).sum()
        f.write(f"{b:.3f} {(fo-fu)/(fu-f0)*100:.3f}\n")
