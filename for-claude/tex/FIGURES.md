# Figure registry (grant-funding paper)

Per the paper-figures skill (step 0). Document: amsart, canonical preamble, letter,
4 cm side margins; full slot = \textwidth = 5.35 in (13.59 cm). Figure text: \small
tick and direct labels, captions \small. All current figures are native pgfplots
tikzpictures \input into main.tex, so they inherit the document font at true size
(no scaling). Tool choice settled by Aydin 2026-09-05: pgfplots throughout this
paper, with styles matched across all figures.

| slug | source | slot | family | status |
|---|---|---|---|---|
| frontier | fig1-frontier.tex (+ fig1-funded.dat, fig1-unfunded.dat) | full (0.82\textwidth axis) | analytic, seed-659 population | placed, Fig 1, \S4 |
| targeting-value | fig2-targeting-value.tex (+ fig2-curve.dat, make-fig2-data.py) | full (0.82\textwidth axis) | analytic, seed-659 population | placed, Fig 2, \S4 |
| thin-grants | fig3-thin-grants.tex (+ fig3-curves.dat) | full (0.82\textwidth axis) | analytic, schematic axes (Fig 1 convention) | placed, Fig 3, \S5 |
| review-value | fig4-review-value.tex (+ fig4-value.dat) | full (0.82\textwidth axis) | simulation, fine D-4 grid (17 tau x 3 alpha x 50 seeds; run_D4_fine.R, model.R hash 21b0d9a, bit-identity validated) | built, Fig 4, \S6, awaiting section verdict |

Family invariants: one population (n = 75, K, R0 ~ Pareto(alpha = 2, min 1), rho = 0,
A = 1, numpy default_rng(659); seed chosen by search for Aydin's requested features:
~1/5 funded at b = 0.1 with slope 1/c = 1.20, no high-K high-R outlier, several
over-resourced low-capability researchers, two big-gap funded researchers); black data ink, gray secondary, no grid, axis lines left
and bottom only, direct labels, no legends. Data regenerable: make-fig2-data.py; the
frontier data snippet is recorded in state.md (2026-09-05), b = 0.1, c = 0.833.

Simulation-figure conventions (set by Fig 4): log x-axis where the sweep is
logarithmic, plain numeric tick labels (no powers), lines without marks through the
measured grid (grid density and seeds stated in the caption), direct labels at the
right edge, black ink throughout. Seed counts follow the generating run (Fig 4: 50
seeds per point, the D-4 convention).

Owed: simulation figures for \S7-\S9 in this style (data exported from the R runs as
.dat files).
