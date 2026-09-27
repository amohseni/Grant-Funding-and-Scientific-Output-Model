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
| records-rounds | fig5-records-rounds.tex (+ fig5-rounds.dat) | full (0.82\textwidth axis) | simulation, T=20 trajectory (verify_s6_rounds.R, 50 seeds) | placed, renders as Fig 4, \S6 |
| review-value | fig4-review-value.tex (+ fig4-value.dat) | full (0.82\textwidth axis) | simulation, fine D-4 grid (17 tau x 3 alpha x 50 seeds; run_D4_fine.R, model.R hash 21b0d9a, bit-identity validated) | placed, renders as Fig 5, \S6 |
| overtrust | fig6-overtrust.tex (+ fig6-overtrust.dat) | full (0.82\textwidth axis) | simulation, D-2 grid (200 seeds; error bars 1 SE, the family's first, justified by the within-noise caveat) | placed, Fig 6, \S6 |
| depth-schedule | fig7-depth-schedule.tex (+ fig7-schedule.dat) | full (0.82\textwidth axis) | simulation, exploration_depth cell re-run in-session (b=3: 24 seeds, SE<=0.006, round-1 share matches canonical 200-seed run) | placed, Fig 8, \S7 |
| timing-boundary | fig8-timing-boundary.tex (+ fig8-boundary.dat) | full (0.82\textwidth axis) | simulation, resource_regime map (200 seeds) | placed, Fig 7, \S7 |
| schedule-vs-signal | fig9-schedule-vs-signal.tex (data inline) | full (0.82\textwidth axis) | simulation, horizon_growth + horizon_long (64-200 seeds) | placed, Fig 9, \S7 |
| floor-cost | fig10-floor-cost.tex (+ fig10-floor-cost.dat) | full (0.82\textwidth axis) | simulation, D3_seed_signal (200 seeds) | placed, Fig 10, \S8 |
| lottery-prices | fig11-lottery-prices.tex (+ fig11-heavy.dat, fig11-even.dat) | two panels, 0.46\textwidth each (the family's first panel pair; marks identified in caption) | single-round pricing computation (price_lotteries.py, 2000 pops/point, MC SE <= 0.002) | placed, Fig 11, \S8 |
| concentration | fig13-concentration.tex (+ fig13-heavy01.dat, fig13-def1.dat, fig13-even1.dat) | 0.82\textwidth; change in Gini from each curve's Cobb-Douglas level; Leontief as detached marks at x=-14; "our form" dotted at gamma=-1; no zero rule (direct label needs the space) | sigma_tierB_concentration + Leontief endpoint + 400-seed refinement, single round, 200 seeds (verify_s9.R) | placed, RENDERS AS Fig 12, \S8 (slug historical) |
| signal-robustness | fig12-signal-robustness.tex (+ fig12-left.dat, fig12-right.dat) | two panels, 0.46\textwidth each; marks identified in caption (filled = Cobb-Douglas, open = gamma -3, triangles = Leontief) | sigma_tierA_gc0/_gcm3/_leontief signal_value summaries, two rounds, 100 seeds (50 Leontief) | placed, RENDERS AS Fig 13, Appendix C |

NOTE on numbering: file slugs are historical; rendered figure numbers follow placement
order (\S6 places records-rounds before review-value, so fig5-*.tex renders as Figure 4
and fig4-*.tex as Figure 5; \S7 will place timing-boundary before depth-schedule).
Renaming files to match was skipped deliberately: the device bridge cannot delete, so
renames would leave stale duplicates in the repo. Trust \ref, not slugs.

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

Table 2 (Appendix C): scheduling value and schedule center of mass by technology, T=5
(sigma_tierA_gc0/_gcm3 horizon_growth summaries).

No figures owed for \S\S4-8 or Appendices A-C. Label fig:lottery-prices renamed
fig:lottery-schemes (file slug kept).
