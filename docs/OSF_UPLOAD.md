# OSF Upload Plan (split: Materials + Data)

A concrete checklist for archiving this study on the Open Science Framework as **one project with
two components** (per the chosen "split data / materials" structure). Most of this is mechanical;
the paste-ready text is below.

## Before you upload — two must-dos

1. **⚠ Materialize the iCloud files.** This repo is on iCloud Drive, so many files are 0-byte
   placeholders until downloaded. From the project root, run:
   `find . -type f -not -path './.claude/*' -exec cat {} + > /dev/null`
   (or Finder → right-click the folder → **Download Now**). Verify with
   `find . -name '*.rds' -size 0` returning nothing.
2. **Push a clean commit to GitHub** (the remote already exists). `.gitignore` keeps `.claude/`
   (agent worktrees/plans) and R session cruft out. OSF can then link the GitHub repo directly.

## OSF project structure

**Project title:** *Grant Funding and Scientific Output: A T-Round Bayesian Model of Research
Allocation* (match the paper).

**Component 1 — "Materials" (code, app, docs)** — license: MIT.
Contents:
- `model.R`, `app.R`, `sweep.R`, `sweep_T.R`  (the model engine, Shiny app, sweep infrastructure)
- `tests/`  (validation + diagnostic scripts)
- `figures/make_figures.R`  (regenerates the figures)
- `reproduce.R`  (one-command end-to-end regeneration)
- `README.md`, `ROADMAP.md`, `RESULTS.md`, `docs/` (PROGRESS, STATE_OF_PLAY equivalent, ENVIRONMENT,
  assertions_log, benchmark_report), `CITATION.cff`, `LICENSE`
- `T_round_extension/{code,validation,diagnostics,report}/`  (the self-contained analysis bundle,
  minus its `data/` which goes to the Data component — or keep the bundle whole and cross-link)

**Component 2 — "Data" (sweep results)** — license: CC-BY-4.0.
Contents:
- `T_round_extension/data/`  — the canonical corrected dataset (16 sweeps × {summary, raw, rawlong}
  `.rds`, plots) — **this is the primary data deposit**
- `T_round_extension/DATA_DICTIONARY.md`  — the schema/codebook (upload alongside the data)
- (optional) `figures/*.pdf` — the rendered publication figures
- **Exclude** `sweep_results/legacy/` (superseded runs incl. the pre-bugfix data) and
  `sweep_results/T_run/` — these are provenance-only, not for the archive.

> The simplest mapping: the `T_round_extension/` package already *is* a self-contained
> materials+data bundle. If you prefer one component, upload it whole and note the dual license in
> the wiki. The split above is the more conventional OSF layout.

## Steps

1. Create the OSF project; add collaborators.
2. Add two components ("Materials", "Data") with the licenses above.
3. **Link the GitHub repo** to the Materials component (OSF → Add-ons → GitHub), so code stays in
   sync. Upload the `.rds` data to the Data component directly (GitHub isn't ideal for large binaries).
4. Fill the project wiki with the paste-ready description below; add the data dictionary link.
5. Add the citation (from `CITATION.cff`) and, once you have co-authors/ORCIDs, update it.
6. When the results are frozen, **register** the project (creates a timestamped, immutable snapshot
   with a DOI) — do this at or after submission.
7. Put the OSF DOI back into `CITATION.cff` and the paper.

## Paste-ready project description (wiki / OSF summary)

> A simulation study of how a science funder should allocate grants across researchers and over T
> funding rounds to maximize scientific output. Researchers differ in latent knowledge (K) and
> resources (R); per-round output is λ = γ·KR/(K+R) (bottlenecked by the scarcer input), and
> knowledge compounds between rounds. We compare nine allocation strategies — naive, myopic-Bayesian,
> and forward-looking (certainty-equivalent) planners — that infer researchers' hidden types from
> noisy publication counts and two signals. This OSF project archives the model code, an interactive
> Shiny explorer, the full parameter-sweep dataset (16 sweeps, 200 trials per cell), and a
> validation/reproducibility suite. **Materials** (code) are MIT-licensed; **Data** are CC-BY-4.0.
> Reproduce everything with `Rscript reproduce.R` (see `docs/ENVIRONMENT.md` for versions).

## Data-component readme (paste into the Data component)

> Corrected T-round parameter-sweep dataset. Each of the 16 sweeps has three files:
> `<sweep>_summary.rds` (one row per parameter cell, with `_mean`/`_se` for every metric — the main
> file), `<sweep>_raw.rds` (one row per trial), `<sweep>_rawlong.rds` (per-round funding schedule).
> Load with R `readRDS()`. Full schema, all sweep grids, and parameter definitions are in
> `DATA_DICTIONARY.md`. Provenance and the (superseded, excluded) pre-bugfix run are described in the
> repository's `STATE_OF_PLAY.md`.
