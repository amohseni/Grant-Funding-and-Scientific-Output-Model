# archive/ — older model versions & scripts

Superseded code kept for history. None of this is used by the current model; the live code is at
the repo root (`model.R`, `app.R`, `sweep.R`, `sweep_T.R`).

## models/
- **`simulation_v3.R`** — the **v3** Shiny app (2-round, older forward planner). Superseded by
  `app.R` (v5), which is itself the T=2 special case of the current default T-round model (`model.R`). Kept to trace the model's evolution.

## sweep_scripts/
Older parameter-sweep scripts, predating the current `sweep.R` / `sweep_T.R` infrastructure:
- `parameter_sweep.R` — an early micro-sweep script.
- `sweep_11-05-2026.R`, `sweep_12-04-2026.R` — dated iterations of the sweep infrastructure
  (v5-era rewrites). Filenames de-colon'd from the originals for macOS friendliness.

## Model lineage (newest → oldest)
`model.R` (T-round, **default**) → `app.R` (v5, 2-round = T=2 special case) →
`archive/models/simulation_v3.R` (v3, 2-round).
