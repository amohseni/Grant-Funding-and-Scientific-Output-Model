"""Numerical verification of the gap rule (Proposition 1).

For random populations (K, R0 Pareto-drawn), compare:
  (a) the gap-rule allocation g_i = max(c K_i - R_i0, 0), c found by bisection on
      the budget identity, against
  (b) a direct numerical solution of the constrained maximization (scipy SLSQP).

Checks: allocations agree; objective values agree; equal-marginal-value condition
holds at the formula allocation.
"""
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(20260814)


def phi(g, K, R0, A):
    return 2 * A * K * (R0 + g) / (K + R0 + g)


def gap_rule(K, R0, B):
    """Bisection on c for sum(max(c K - R0, 0)) = B."""
    def G(c):
        return np.maximum(c * K - R0, 0.0).sum()
    lo, hi = 0.0, 1.0
    while G(hi) < B:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if G(mid) < B:
            lo = mid
        else:
            hi = mid
    c = 0.5 * (lo + hi)
    return np.maximum(c * K - R0, 0.0), c


def direct_opt(K, R0, B, A):
    n = len(K)
    x0 = np.full(n, B / n)
    cons = [{"type": "eq", "fun": lambda g: g.sum() - B}]
    bounds = [(0.0, B)] * n
    res = minimize(lambda g: -phi(g, K, R0, A).sum(), x0, bounds=bounds,
                   constraints=cons, method="SLSQP",
                   options={"maxiter": 2000, "ftol": 1e-14})
    return res.x, -res.fun


def run_trial(n, alpha_K, alpha_R, b, A, trial):
    K = (rng.pareto(alpha_K, n) + 1.0)          # Pareto with x_m = 1
    R0 = (rng.pareto(alpha_R, n) + 1.0)
    if trial % 3 == 0:                          # exercise the R0 = 0 boundary
        R0[rng.integers(0, n)] = 0.0
    B = b * R0.sum() if R0.sum() > 0 else b * n
    g_formula, c = gap_rule(K, R0, B)
    g_direct, f_direct = direct_opt(K, R0, B, A)
    f_formula = phi(g_formula, K, R0, A).sum()

    # equal-marginal check at the formula allocation
    mv = 2 * A * K**2 / (K + R0 + g_formula) ** 2
    funded = g_formula > 1e-9
    nu = 2 * A / (1 + c) ** 2
    mv_funded_dev = np.abs(mv[funded] - nu).max() if funded.any() else 0.0
    mv_unfunded_ok = (mv[~funded] <= nu + 1e-9).all()

    return {
        "obj_gap": f_direct - f_formula,          # >0 would mean formula suboptimal
        "alloc_dev": np.abs(g_formula - g_direct).max() / max(B, 1.0),
        "budget_dev": abs(g_formula.sum() - B) / B,
        "mv_funded_dev": mv_funded_dev,
        "mv_unfunded_ok": mv_unfunded_ok,
        "n_funded": int(funded.sum()),
        "c": c,
    }


results = []
configs = []
for alpha in [1.3, 2.0, 3.5]:
    for b in [0.1, 0.5, 1.0]:
        for A in [0.5, 1.0, 2.0]:
            configs.append((alpha, b, A))

for t, (alpha, b, A) in enumerate(configs):
    for rep in range(3):
        r = run_trial(50, alpha, alpha, b, A, t * 3 + rep)
        r.update(alpha=alpha, b=b, A=A)
        results.append(r)

worst_obj = max(r["obj_gap"] for r in results)
worst_alloc = max(r["alloc_dev"] for r in results)
worst_budget = max(r["budget_dev"] for r in results)
worst_mv = max(r["mv_funded_dev"] for r in results)
all_unfunded_ok = all(r["mv_unfunded_ok"] for r in results)

print(f"trials: {len(results)} (n=50 each; alpha in 1.3/2.0/3.5, b in 0.1/0.5/1.0, A in 0.5/1/2)")
print(f"max objective shortfall of formula vs direct optimizer: {worst_obj:.3e}")
print(f"  (negative means the formula BEAT the numerical optimizer)")
print(f"max allocation deviation (relative to B): {worst_alloc:.3e}")
print(f"max budget identity violation (relative): {worst_budget:.3e}")
print(f"max |marginal value - nu| among funded: {worst_mv:.3e}")
print(f"unfunded marginal values all <= nu: {all_unfunded_ok}")
fr = [r["n_funded"] for r in results]
print(f"funded counts ranged {min(fr)}..{max(fr)} of 50 (frontier active in small-b trials)")

# Corollary 1(ii): c strictly increasing in B on a fixed population
K = rng.pareto(2.0, 50) + 1.0
R0 = rng.pareto(2.0, 50) + 1.0
cs = [gap_rule(K, R0, b * R0.sum())[1] for b in [0.05, 0.1, 0.25, 0.5, 1.0, 2.0]]
print(f"c(B) on fixed population, b = .05,.1,.25,.5,1,2: "
      + ", ".join(f"{c:.4f}" for c in cs)
      + f"  strictly increasing: {all(x < y for x, y in zip(cs, cs[1:]))}")
