# Roadmap to Publication

Concrete, ordered steps to take this project from its current state (core results in, model &
app validated) to a finished, shareable, OSF-archived study ready to write up. Each step is
tagged **[you]** or **[assistant]** (things Claude can execute) and marked with rough effort.

> ## ✅ Status (all assistant-executable steps complete)
> - **Stage 1 — DONE.** 3 supplementary sweeps run (`tail_map`, `info_value`, `horizon_long`).
>   `horizon_long` was re-run at converged granularity after a diligence check caught that
>   `n_steps=50` inflates the high-T advantage ~4× (T=1–5 unaffected). CE robustness handled via the
>   lightweight check: CE gap < 0.4% across T=2–5, so the SP-lite planner was **not needed**.
> - **Stage 2 — DONE.** 16/16 sweeps consolidated into `T_round_extension/data/`; `DATA_DICTIONARY.md`
>   updated; publication figures (`figures/`, 8× PNG+PDF) generated; `RESULTS.md` written.
> - **Stage 3 — DONE (assistant parts).** `reproduce.R`, `LICENSE` + `LICENSE-DATA.md`, `CITATION.cff`
>   (needs your ORCID/co-authors), `.gitignore`, `docs/ENVIRONMENT.md`, `docs/OSF_UPLOAD.md`.
> - **Stages 4–5 — YOURS:** commit/push, redeploy the app, create the OSF project, write the paper.
>
> The remaining action items that only you can do are listed under "Decisions/actions" per stage
> below (ORCID + co-authors in `CITATION.cff`, GitHub push, Shiny redeploy, OSF creation).

**Where we are now:** the T-round model is built, validated (bit-identical to the 2-round v5 at
T=2), and a subtle planner bug was found and fixed. The full 13-sweep manifest has been run
(corrected, 200 trials/cell) and packaged in `T_round_extension/`. The Shiny app is refactored to
the T-round model. What remains is (1) a few supplementary sweeps to close analysis gaps, (2)
paper-ready figures and stats, and (3) reproducibility + OSF packaging.

---

## Stage 1 — Complete the results *(gets "all data in")*

**1.1 Run the supplementary sweeps.** [assistant] · ~15–25 min compute
Defined in `T_round_extension/NEXT_SIMULATIONS.md`. Recommended minimum (≈31 cells):
- `tail_map` — resource-side inequality (the one real gap; may qualify the "inequality" framing)
- `info_value` — isolates the information channel at ε→0 (backs the compounding-vs-info decomposition)
- `horizon_long` — does the forward advantage saturate past T=5?

Optional (robustness, add if reviewers will ask): `signal_complementarity`, `horizon_noise_hieps` /
`horizon_scale_hieps`, `pop_horizon`. → writes to `sweep_results/T_run_supplement/`, same schema.

**1.2 Decide the CE-approximation robustness scope.** [you decide] · then [assistant] if yes
The deferred **SP-lite planner** check (`STATE_OF_PLAY.md` §6e) would bound the certainty-equivalence
error at high T with a scenario-based reference. The CE gap is already ~0.24% at T=2 and results are
granularity-stable, so this is *optional rigor* for a methods appendix. **Decision: build it, or cite
the existing CE-vs-scenario assertion (#9) and scope it out.**

---

## Stage 2 — Analysis & figures *(paper-ready)*

**2.1 Consolidate all sweeps.** [assistant] · ~10 min
Fold `T_run_supplement/` into the analysis package (`T_round_extension/data/`), refresh
`DATA_DICTIONARY.md` with the new sweeps.

**2.2 Produce publication figures.** [assistant] · ~30–60 min
A dedicated `figures/` set (PNG + PDF, consistent theme) for the paper's core claims:
- The (T×ε) phase diagram (forward vs myopic) — the headline
- The signal-value law heatmap (inequality × precision) + a fitted curve
- Correlation robustness; seeding-never-helps; the new `tail_map`, `info_value` decomposition,
  `horizon_long` saturation
Deliver as standalone script(s) that regenerate every figure from the `.rds` data.

**2.3 Statistical results write-up.** [assistant drafts → you verify] · ~30–60 min
A `RESULTS.md` capturing, per claim: effect sizes with SE/z, regime-boundary locations, and simple
parametric fits (e.g. `signal_fwd` as a function of α_K and τ_K; the ε-slope vs ε→0 intercept of the
forward gain). This is the skeleton you write the paper's Results section from.

---

## Stage 3 — Reproducibility & OSF packaging *(shareable)*

**3.1 One-command reproduce script.** [assistant] · ~15 min
A top-level `reproduce.R` that regenerates the entire dataset from seeds (`main_sweep_T` + the
supplements) and re-runs the validation suite, so a collaborator can verify everything end-to-end.
Environment is recorded in `docs/ENVIRONMENT.md` (done).

**3.2 License & citation.** [you choose → assistant writes] · ~10 min
- **Pick a license** (common for research code/data: MIT for code + CC-BY-4.0 for data, or a single
  permissive choice). → I'll add `LICENSE`.
- Provide author list / ORCID / intended citation. → I'll add `CITATION.cff` (GitHub/OSF read it).

**3.3 Repo hygiene.** [assistant] · ~10 min — partly done
- `.gitignore` excludes `.claude/` (agent worktrees, plans), `.DS_Store`, R session files (done).
- **Decide on the legacy data** (`sweep_results/legacy/`): keep it in-repo (labeled, small) or move it
  out of the OSF upload. Recommendation: keep in GitHub for provenance, exclude from the OSF data
  component (OSF should hold only the canonical corrected data).

**3.4 ⚠ Materialize iCloud files before uploading.** [you] · a few min
This project lives on **iCloud Drive**, so many files are "dataless" placeholders (they show 0 bytes
until downloaded). **Before uploading to OSF or zipping**, force-download everything: in Finder,
right-click the project folder → **"Download Now"**, or run
`find . -type f -exec cat {} + > /dev/null` from the project root. Otherwise the archive will contain
empty files.

**3.5 Build the OSF structure.** [you decide scope → assistant assembles]
Two clean options:
- **(a) Whole cleaned repo** to OSF (mirrors GitHub) — simplest, keeps code+data+docs together.
- **(b) Split** — OSF *data component* = the corrected sweep data + dictionary; OSF *materials* = code
  + docs; link the GitHub repo. More conventional for OSF.
Recommendation: **(a)** now (fast), optionally reorganize to (b) later. The `T_round_extension/`
package is already a self-contained, documented bundle that maps cleanly to an OSF component.

---

## Stage 4 — Share & deploy

**4.1 Commit & push to GitHub.** [you, or assistant on request] · the remote already exists
(`github.com/amohseni/Grant-Funding-and-Scientific-Output-Model`). A single clean commit of the whole
T-round work + reorganization. *(No commits have been made during this work — the tree is staged for
your review.)*

**4.2 Redeploy the Shiny app.** [you] · ~5 min
`app.R` now sources `model.R`; the deploy bundle must include both. `rsconnect::deployApp()` (or the
RStudio "Publish" button). Also fixes the earlier `ce_reweight` bug in the live app.

**4.3 Create the OSF project.** [you] · link GitHub, upload the data component, add the wiki/README,
and (when ready) register the project to timestamp it.

---

## Stage 5 — Write the paper *(you)*

Draft from `RESULTS.md` (2.3) and the `figures/` (2.2). The narrative arc the results support:
forward planning beats myopic and by more as knowledge compounds and the horizon lengthens
(front-loading mechanism); the peer-review signal's value is governed by inequality × precision and
survives realistic K–R correlation; uniform seeding never helps. The methods appendix covers the
CE planner, the T=2 = v5 validation, and (if 1.2) the CE-vs-SP bound.

---

## Decisions I need from you to proceed

1. **Supplementary sweeps:** run the recommended 3, or all 6? (Stage 1.1)
2. **SP-lite robustness check:** build it, or scope it out with the existing assertion? (Stage 1.2)
3. **License:** which one for code/data? (Stage 3.2)
4. **OSF scope:** whole-repo (a) or split data/materials (b)? (Stage 3.5)

Give me those and I'll execute Stages 1–3 (sweeps → figures → reproduce script → license/citation →
packaging) in sequence; Stages 4–5 are yours (push/deploy/upload/write), and I can prep each.
