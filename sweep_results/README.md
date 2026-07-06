# sweep_results/ — raw parameter-sweep outputs

Raw `.rds` / `.png` outputs written by the sweep runner. **For analysis, prefer the cleaned,
documented copy in [`../T_round_extension/data/`](../T_round_extension/data/)** (same files, with
a data dictionary).

## Current (canonical)
- **`T_run_fixed/`** — the **corrected T-round manifest**: 16 sweeps × {summary, raw, rawlong}
  `.rds` + plots, 200 trials/cell, M=400, from the **fixed** planner. This is the real dataset.
  (`_run_launch.log` inside it is that run's console log.) Column meanings:
  `../T_round_extension/DATA_DICTIONARY.md`.

## legacy/ — superseded, kept for provenance (do not analyze)
| Folder / file | What it is | Why superseded |
|---------------|-----------|----------------|
| `T_round_buggy_M200/` | The **first** T-round manifest run (M=200) | Ran on the pre-fix planner with the `ce_reweight` SIR bug — high-T results are artifacts. See `../T_round_extension/STATE_OF_PLAY.md` §5. |
| `sweep_2round_11-05-2026/` | An older **2-round (v5)** parameter sweep | Predates the T-round extension |
| `micro_run_12-04-2026/` | A small early micro-sweep iteration | Exploratory |
| `plots_DHA_inner/`, `plots_DHA_toplevel/` | Old hypothesis plots (hyp_A … hyp_H) | Earlier analysis pass |
| `old_csv/` | CSV outputs (`forward_vs_myopic_main_spotcheck`, `optimal_alpha_v1/v2`) | The CSVs referenced by `../docs/CHAT_HANDOFF_PROMPT.md`; superseded by the `.rds` manifest |
| `pg_focused_summary.rds` | A focused PG(T) probe | Pre-fix; superseded by `horizon_growth` in `T_run_fixed/` |
| `logs/` | Console logs from superseded runs | — |

> New sweeps: the runner writes wherever you set `out_dir` (e.g. the proposed supplements in
> `NEXT_SIMULATIONS.md` target `sweep_results/T_run_supplement/`).
