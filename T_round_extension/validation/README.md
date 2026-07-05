# validation/ — proof the code is correct

Evidence that the T-round model is faithful and the results are trustworthy. Run any `.R` script
from the **project root** (`Rscript validation/<script>.R`).

| File | What it proves / does |
|------|-----------------------|
| `test_T2_reduction.R` | **The anchor.** T-round model at T=2 == v5 bit-identically (max |Δ| ≈ 4e-14) across 7 param points × 5 seeds × 9 strategies. The core correctness guarantee. |
| `golden_master_v5.rds` | Snapshot of v5 per-strategy output at several param points; the reference the anchor test compares against. Re-captured after the fix. |
| `capture_golden_master.R` | Regenerates `golden_master_v5.rds` from the current `app.R`. |
| `assertions_T.R` | The §6 validity battery (10 checks): T=1 triviality, ε→0 limit, τ_K limits, heavy-tail finiteness + ESS, budget conservation, M-convergence, CE-vs-scenario gap, sign-path. |
| `assertions_log.txt` | Output of the above — **all pass**. Read this for the one-line status of each check. |
| `benchmark_T.R` | Measures per-run cost vs T and projects full-manifest wall-clock. |
| `benchmark_report.txt` | Its output: O(T²) scaling; full manifest ~30 min at M=400 on 7 cores. |
| `launch_manifest_fixed.R` | **The exact production driver** used to generate `../data/` (200 seeds, M=400, fixed code). |

**Headline validation facts:**
- T=2 reproduces the (fixed) v5 model exactly.
- All §6 assertions pass; CE-vs-scenario valuation gap at T=2 is −0.26% (the CE approximation is
  near-exact there).
- Forward advantage is stable as greedy granularity is refined (post-fix) — see `../diagnostics/`.
