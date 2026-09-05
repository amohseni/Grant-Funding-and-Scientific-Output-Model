# Context transfer for the drafting session (2026-08-06)

Copy of the handoff block given in chat; the new session should read this file plus state.md first.

OBJECTIVE. Write the grant-funding paper (Mohseni + Simon + Kevin Zollman; Templeton
deliverable) step-by-step to conclusion: venue final check, outline, findings
consolidation, intro, model spec, objections/implications, conclusion, appendices.
Invoke the paper-writing skill; read drafts/grant-funding/state.md FIRST. It is the
paper's tracker and carries the full findings register, resolved decisions, and
do-not-claim list.

FIRST TASK (Aydin's instruction): build complete command of all simulations before
any prose. Read, in the model repo Grant-Funding-and-Scientific-Output-Model
(branch smooth-allocator): docs/PAPER_INTEGRATION_HANDOFF_2026-08-06.md (CANONICAL:
20/20 preregistered verdicts, file map, 9 drafting flags, integration order), then
docs/DIAGNOSTICS_RESULTS_2026-08.md and T_round_extension/RESOURCE_REGIME_RESULTS.md.
Any number re-derivable: Rscript sweep_results/_probe/verify_all_claims.R (needs R;
spot-run owed). Root RESULTS.md is stale; ignore where they conflict.

THE PAPER. Story 1 locked: fund the talent-resource gap. Closed form g* = cK - R
(water-filling; cite as precedent, claim only first statement for science funding).
Thesis: the value of peer review is the cost of not knowing the gap; now measured
(D-4: corr 0.34 -> 0.97, oracle gap 17.2 -> 1.7% of S1; the thesis figure). Type D
primary + B engine + E outer layer. Architecture: puzzle/two-rules-refuted -> rule +
frontier -> obstacle (K unobserved; thin grants talent-uninformative) -> peer review
(inequality x precision; AUC calibration tau_K > 20; overtrust caveat) -> timing
(demoted; back-loading, bootstrap subsection with seed-and-harvest figure
exploration_depth_schedules.png; scope to gamma_ces < 0) -> floors (price law: cost =
targeting value x floored share; engage lotteries where price lowest) -> robustness
(budget-conditional concentration law + Cobb-Douglas boundary; these are results) ->
discussion (free dominates paid; B regenerates C, D regenerates E; one quantity,
dispersion of marginal returns, prices everything). Base sources: Aydin's report +
model write-up PDFs ("Optimal_Funding_Strategies_for_Scientific_Output" x2; ask him
to re-attach; LaTeX sources not in either repo). Apply
drafts/grant-funding/spec-errata.md (12 items: notation A/gamma_ces/sigma, Eq 10
budget fix, defaults table, TBD-to-Kevin flags).

CONSTRAINTS. No em dashes ever. One concept one term. Findings-forward intro;
conditional theses; opposition at full strength. Venue realistic range: Research
Policy down through philosophy of science (general science requires calibration
that does not exist). Robustness-not-run claims carry TBD notes to Kevin where
flagged. Do-not-claim list in state.md governs every sentence.

OPEN DECISIONS (Aydin's): venue final; K terminology (talent recommended vs
knowledge). Word-level choices always his; critique before revision; ranked
candidates, never silent substitutions; gate each section on his approval.

State-folder contract: update state.md at every milestone; sessions push directly
(queue in memory /areas/for-claude-writeback.md only on push failure).
