# docs/ — process & metadata documents

| File | What it is |
|------|-----------|
| `PROGRESS.md` | The **running development log** for the T-round extension — chronological, detailed, includes the full bug story. The clean, organized version of this status lives in `../T_round_extension/STATE_OF_PLAY.md`. |
| `CHAT_HANDOFF_PROMPT.md` | **Historical.** A prompt written for an earlier analysis handoff ("Sweep B / optimal alpha") that referenced CSV outputs now archived under `../sweep_results/legacy/old_csv/`. Superseded by the current T-round dataset and its dictionary; kept for reference. |
| `assertions_log.txt` | Output of the §6 validity assertion suite on the **current** (fixed) model — all checks pass. Regenerate with `../tests/assertions_T.R`. |
| `benchmark_report.txt` | Runtime benchmark: per-run cost vs horizon T and the full-manifest wall-clock projection. |

For the model, findings, and how to use the data, start at the repo root `README.md` and
`../T_round_extension/`.
