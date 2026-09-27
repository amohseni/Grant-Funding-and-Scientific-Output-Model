# S8 mechanism computation: single-round expected-output prices of egalitarian
# schemes, exact given the allocation (no output noise enters expectations).
# Population: n=50, K,R0 ~ Pareto(alpha_K / alpha_R=2, min 1), A=1. Budget B = b*n*E[R],
# E[R]=2. Signal s = K + N(0, tau). Schemes: gap rule (complete info; benchmark),
# screened equal split (top-m by s, B/m each, m=q*n), screened lottery (pool top-2m
# by s, m funded at random, B/m each; exact expectation over the draw), uniform (B/n),
# full lottery (m winners of B/m, exact expectation). Metric: (scheme - no-funding) /
# (gap rule - no-funding), averaged over 2000 populations.
import numpy as np
rng = np.random.default_rng(20260908)
n, q, npop = 50, 0.2, 2000
m = int(q*n)
def phi(K, R):  # expected output, A=1
    return 2*K*R/(K+R)
def gap_alloc(K, R0, B):
    lo, hi = 1e-9, 1e9
    for _ in range(200):
        c = 0.5*(lo+hi)
        s = np.maximum(c*K - R0, 0).sum()
        if s > B: hi = c
        else: lo = c
    return np.maximum(0.5*(lo+hi)*K - R0, 0)
def run(alpha, b, tau):
    out = np.zeros(5)
    base = 0.0; opt = 0.0
    for _ in range(npop):
        K  = (1 + rng.pareto(alpha, n))
        R0 = (1 + rng.pareto(2.0, n))
        B  = b*n*2.0
        f0 = phi(K, R0).sum()
        fstar = phi(K, R0 + gap_alloc(K, R0, B)).sum()
        s = K + rng.normal(0, tau, n)
        order = np.argsort(-s)
        # screened equal split: top-m get B/m
        g = np.zeros(n); g[order[:m]] = B/m
        f_ses = phi(K, R0+g).sum()
        # screened lottery: pool top-2m, each funded w.p. 1/2 with B/m
        pool = order[:2*m]
        f_slot = phi(K, R0).sum() + 0.5*(phi(K[pool], R0[pool]+B/m) - phi(K[pool], R0[pool])).sum()
        # uniform
        f_u = phi(K, R0+B/n).sum()
        # full lottery: each researcher funded w.p. m/n with B/m
        f_flot = f0 + (m/n)*(phi(K, R0+B/m) - phi(K, R0)).sum()
        denom = fstar - f0
        out += np.array([f_ses-f0, f_slot-f0, f_u-f0, f_flot-f0, fstar-f0])/denom
    return out/npop
print(f"{'alpha':>5} {'b':>4} {'tau':>4} | {'scr-split':>9} {'scr-lott':>9} {'uniform':>8} {'full-lott':>9}")
for alpha in (1.3, 2.0, 3.5):
    for b in (0.1, 0.5):
        for tau in (0.3, 1, 3, 10):
            r = run(alpha, b, tau)
            print(f"{alpha:5.1f} {b:4.1f} {tau:4.0f} | {r[0]:9.3f} {r[1]:9.3f} {r[2]:8.3f} {r[3]:9.3f}" if tau>=1 else
                  f"{alpha:5.1f} {b:4.1f} {tau:4.1f} | {r[0]:9.3f} {r[1]:9.3f} {r[2]:8.3f} {r[3]:9.3f}")

# --- Extension (same session): wide split (pool 2m, B/2m each) added to isolate the
# within-pool randomization price; MC SEs reported; cells restricted to the displayed
# (alpha, b) pairs. See verify_s8_OUTPUT.txt addendum for the run and its output.
