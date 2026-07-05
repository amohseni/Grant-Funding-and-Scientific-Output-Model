# Compute Environment (for reproducibility)

The results were produced with the following environment. Reproductions should use matching
major versions; the model is pure base R plus a few well-established packages, so exact patch
versions are not critical.

| Component | Version |
|-----------|---------|
| R | 4.4.2 (2024-10-31) |
| Platform | aarch64-apple-darwin20 (Apple Silicon, macOS) |
| shiny | 1.10.0 |
| ggplot2 | 3.5.1 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| purrr | 1.0.4 |
| scales | 1.3.0 |
| parallel | 4.4.2 (base) |

**Dependencies by component:**
- `model.R` — **base R only** (no packages). This is the model engine.
- `sweep.R` / `sweep_T.R` — `dplyr`, `tidyr`, `ggplot2`, `scales`, `parallel`.
- `app.R` — the above plus `shiny`.

**Reproducibility notes:**
- Every run seeds its own RNG from the `seed` argument (`set.seed(seed)` inside `run_simulation_T`);
  sweeps use `seed = trial index`, shared across cells (common random numbers).
- Parallelism uses `parallel::mclapply` (fork); results are independent of core count because
  each task is seeded explicitly. Cores default to `detectCores() - 1`.
- Install packages with: `install.packages(c("shiny","ggplot2","dplyr","tidyr","purrr","scales"))`.

To regenerate this file: `Rscript -e 'for (p in c("shiny","ggplot2","dplyr","tidyr","purrr","scales")) cat(p, as.character(packageVersion(p)), "\n")'`.
