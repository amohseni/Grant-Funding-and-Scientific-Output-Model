# Section 7 claims ledger (v1, 2026-09-08)

All rows re-derived 2026-09-08 in-container from the canonical sweep files (staged from
the Mac; verify_s7.R + verify_s7_OUTPUT.txt in this folder), matching
docs/PAPER_INTEGRATION_HANDOFF_2026-08-06.md. Strategy names: S2 = uniform seed,
S5 = myopic with review, S7/S8 = forward without/with review; PG = S8 - S5;
b_idx = budget-schedule center of mass (0.5 = even).

| # | Claim (as formulated in the text) | Modal force | Ground | Status |
|---|---|---|---|---|
| T1 | Up to now the funder has spent equal parts per round | descriptive of \S4-\S6 runs (tranche = B/T) | run_D4.R and defaults | OK |
| T2 | With no resources at all, all-in-first-round and all-in-last-round produce identical expected output | analytic, limiting case | R0 = 0: unfunded output is 0, so no observation or compounding between rounds under either lump; both allocate on the prior; same K, same allocation, same expected lump output | OK (derivation, note 7) |
| T3 | Best above-even-early two-block schedule loses to even: by 1.6% (eps ~ 0) and 3.5% (eps = 0.3) of the even schedule's output | measured, planner-free design; stated in relative terms per the 2026-09-08 rule | bootstrap_verify/honest_schedules_*.csv (raw: 7.8 of 498.2; 20.8 of 602.6) | VERIFIED |
| T4 | Every early-leaning cell has eps <= 0.03; boundary eps* ~ 0.02 across baseline resources; poverty raises the early share (0.17 vs 0.06) and mutes back-loading (b_idx 0.52 vs 0.65) without reversing it | measured over the r_min x eps map, 200 seeds | resource_regime_summary.rds: 10 front-load cells, max eps 0.03; eps=0.85 row | VERIFIED |
| T5 | Min b_idx across the 32 exploration cells = 0.497 (do-not-claim guard: no strict "never front-loads") | precision honesty | exploration_200 summaries | VERIFIED (guard held) |
| T6 | Thin grants: informed funder gains 0.3% of uniform funding's output at standard depth; 9.3% at b=3, 10.0% at b=6 (sharp free signal at standard depth: 9-11%); first-round share 0.110 at b=3 (even 0.167); full schedule in the depth-schedule figure (b=3 re-run in-session, 24 seeds, SEs <= 0.006, round-1 0.104 vs canonical 0.110) | measured, relative terms; poverty + heavy tail + no review signal + T=6 + budget_ref=K context stated | exploration_depth + corner summaries; fig7-depth-schedule | VERIFIED, figure-backed |
| T7 | Attribution: pure grant-fed channel leans mildly early (b_idx 0.481-0.489, SE <= 0.0023); pure free channel leans late (0.529-0.680); free dominates at ~1/8-1/3 of the grant-fed rate | measured, 200 seeds, T=5 | bload_decouple_summary.rds + transects (crossings: eps_free in (0.05,0.1) at eps_paid=0.3; (0.1,0.15) at 0.85) | VERIFIED |
| T8 | Scheduling's gain < the signal's value in ALL 32 tested (T, eps) cells: ratio 0.000-0.075 at eps <= 0.3 (all horizons), max 0.66 at the joint extreme (T=10, eps=0.85); as % of no-funding output, gain 0-0.7% vs signal 0.8-7.5% | sweep-ranging statement per the 2026-09-08 rule (the earlier single-number 1/30 formulation retired) | horizon_growth + horizon_long; fig9-schedule-vs-signal | VERIFIED, figure-backed |
| T9 | PG grows ~quadratically for T <= 5 (exponents 2.3, 2.1), saturates T=5-10 (1.4, 0.8) | fitted exponents, scoped to T<=5 | horizon_growth + horizon_long | VERIFIED |
| T10 | Cobb-Douglas: PG(T=5) <= 0.014, b_idx <= 0.512; gamma=-3: PG 3.77, b_idx 0.660; body claim "exists only where capability and resources are complements" | scoped claim; CD = boundary of the family | sigma_tierA_gc0 / _gcm3 | VERIFIED |
| T11 | CE self-consistency: the deliberate planner's schedule beats even under its own objective 20/20 seeds (+3.33) | verification held in reserve (S8-S5<0 cells excluded from body) | ce_self_consistency.csv | VERIFIED (not cited in body) |
| T12 | "For a funder thinking about how to maximize their impact, choosing whom to fund outweighs choosing when" (his wording) | comparative, licensed by T8 (ratio < 1 in all 32 cells) | T8; fig9 | OK given T8 |

Notes.
- b gloss RESOLVED (draft note 8): the depth cells use budget_ref = "K" (B = b n E[K]),
  so b scales the budget against capability; the footnote now says so.
- Result 5 of the campaign memo (deliberate planner loses to myopic re-deciding at
  depth, S8-S5 = -16 to -45) is deliberately absent from the body; if it enters, the
  honest-CE-mispricing caveat (Boot-2, T11) must attach.
- Poverty non-monotonicity of information value (peaks at r_min ~ 1) is available but
  unused; candidate for \S10 discussion.
