# Funding the Gap: an interactive guide

A Shiny app that walks through the core results of "A model of optimal science funding:
targeting the capability-resource gap" (Mohseni, DeDeo, Zollman) in seven steps, one
interactive figure per step: the production function, the gap rule and the two heuristics,
what records reveal, what review is worth, seed grants and lotteries, a regime map on which
a funder can place its own program, and takeaways with the model's limits.

The analytic parts (harmonic-mean output, the gap rule, the heuristics, the certain-grant
versus lottery comparison) are computed live in base R. Every simulation result is read from
the paper's own figure data in `data/`, copied from `for-claude/tex/` and `for-claude/analysis/`
(fig4-value, fig5-rounds, fig10-floor-cost, fig11-heavy/even, fig16-grid, fig16-targ/rev).
If those change, copy them again; nothing here recomputes them.

Run locally: `shiny::runApp("explainer")` (needs shiny, bslib, ggplot2).
Deploy: `Rscript explainer/deploy.R` from the repository root, with an rsconnect account
configured (`rsconnect::setAccountInfo`). The app name is `funding-the-gap`.

Column note for `data/fig16-grid.dat`: `rev` is review's value as a share of the research
output that review-informed funding adds (the paper's panel B); `revgain` is the share of the
no-review shortfall recovered, which is unstable at very small budgets and is not used.

Figure sizing: each figure renders at a fixed width and aspect ratio chosen for its content (single
plots 640 px wide at 4:3 or 3:2, two-panel plots 720 px), so plot text stays in proportion to the page
text. The container scales a figure down on narrow screens and never up. Below 992 px the figure moves
above its text and controls.
