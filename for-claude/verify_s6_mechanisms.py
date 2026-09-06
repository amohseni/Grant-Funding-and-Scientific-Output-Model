"""Mechanism tests for Section 6 claims (see section-6-claims.md).
C5b: records-only identification decays (K compounds, R fixed; Fisher info about K).
C6: share of optimal-funding output gain from the top 10% by capability.
C7: top-decile recovery vs middle-pair ordering under a noisy signal s = K + N(0, tau).
Run: python3 verify_s6_mechanisms.py
"""
import numpy as np

A, eps = 1.0, 0.1
print("C5b: K compounds with fixed R; lambda -> 2AR; per-round Fisher info about K decays")
K, R = 2.0, 1.5
for t in range(1, 201):
    lam = 2*A*K*R/(K+R); dldK = 2*A*R*R/(K+R)**2
    if t in (1, 10, 50, 100, 200):
        print(f"  t={t:3d}: K={K:7.1f} lam={lam:.3f} (2AR={2*A*R}) Fisher={dldK**2/lam:.2e}")
    K += eps*2*A*K*R/(K+R)
K, tot = 2.0, 0.0
for t in range(20000):
    lam = 2*K*R/(K+R); tot += (2*R*R/(K+R)**2)**2/lam; K += eps*2*K*R/(K+R)
print(f"  cumulative Fisher info, 20000 rounds: {tot:.4f} (finite)")

def solve_c(K, R0, B):
    def G(c): return np.maximum(c*K-R0,0).sum()
    lo, hi = 0.0, 1.0
    while G(hi) < B: hi *= 2
    for _ in range(200):
        m=(lo+hi)/2; lo,hi=(m,hi) if G(m)<B else (lo,m)
    return (lo+hi)/2
phi = lambda K,R: 2*K*R/(K+R)

print("C6: share of optimal-gain from top 10% by K (b=0.1, n=75, 400 draws)")
rng = np.random.default_rng(11)
for alpha in (1.3, 2.0, 3.5):
    shares = []
    for _ in range(400):
        K  = rng.pareto(alpha, 75)+1; R0 = rng.pareto(alpha, 75)+1
        c = solve_c(K, R0, 0.1*R0.sum()); g = np.maximum(c*K-R0, 0)
        gain = phi(K, R0+g) - phi(K, R0)
        shares.append(gain[np.argsort(-K)[:8]].sum()/gain.sum())
    print(f"  alpha={alpha}: {np.mean(shares):.2f}")

print("C7: top-10% recovery vs middle-pair ordering, s = K + N(0,tau) (2000 draws/cell)")
rng = np.random.default_rng(12)
for alpha in (1.3, 2.0, 3.5):
    for tau in (1.0, 3.0, 10.0):
        hit, mid = [], []
        for _ in range(2000):
            K = rng.pareto(alpha, 75)+1; s = K + rng.normal(0, tau, 75)
            hit.append(len(set(np.argsort(-K)[:8]) & set(np.argsort(-s)[:8]))/8)
            m = np.argsort(-K)[18:56]; i, j = rng.choice(m, 2, replace=False)
            hi, lo = (i, j) if K[i] > K[j] else (j, i)
            mid.append(1.0 if s[hi] > s[lo] else 0.0)
        print(f"  alpha={alpha} tau={tau:4.1f}: top recovery {np.mean(hit):.2f}, middle-pair order {np.mean(mid):.2f}")
