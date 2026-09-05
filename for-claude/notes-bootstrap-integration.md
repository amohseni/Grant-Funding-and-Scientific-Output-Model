# Integrating the resource-regime / exploration suite into the paper

2026-08-05. Source: Aydin's 64-seed suite (resource_regime, exploration_corner, exploration_poverty,
exploration_depth; write-up at T_round_extension/RESOURCE_REGIME_RESULTS.md; data
sweep_results/T_run_smooth_supplement/). The main diagnostic suite (Packages A/B/C) is still
running; this memo integrates what exists now and marks what waits.

## 0. The one-sentence integration

The suite answers the strongest natural objection to back-loading (the bootstrap objection) and, in
doing so, extends the paper's thesis rather than qualifying it: **in poor fields even output stops
revealing the talent–resource gap, because thin grants pin output to the grant; the first dollar's
job is to make the gap observable, and money then follows the information it bought.** Aydin's
gloss ("fund a little up front, just enough for an informative signal, then watch") is exactly what
the optimal schedule does at depth, and the paper should say it in nearly those words.

## 1. Claims established, typed and calibrated

- **C-R1 (no reversal).** Across all 54 cells of the suite, b_idx never falls below 0.5 beyond
  noise. The b_idx = 0.5 boundary is vertical at ε* ≈ 0.02 (0.005–0.05 across r_min); below it the
  schedule is FLAT (b_idx ≈ 0.4985–0.502, within noise of even), and timing is worthless there
  (S8−S5 ≈ 0). Poverty mutes back-loading (ε=0.85: b_idx 0.651 at r_min=3 → 0.521 at r_min=0.001;
  round-1 share 0.059 → 0.174) but never flips it. CALIBRATION RULE: the paper claims "there is no
  front-loading regime," not "front-loading exists below ε*." A b_idx of 0.4985 is even-split, and
  claiming a front-load regime there hands a referee a free correction.
- **C-R2 (identification requires depth; new mechanism).** With R0 ≈ 0, λ = γKg/(K+g), so g ≪ K
  pins λ ≈ γg: output reveals the grant, not the talent (g=1/3: K=1 gives λ=0.25, K=100 gives
  0.33). Discrimination value S8−S2 in the pure-exploration corner: +0.9 at b=0.5 vs +17.5 with a
  sharp free signal. Threshold: output becomes talent-limited when g ≳ K; with the decoupled purse
  this is budget scale b ≳ T/2 (per-round even tranche per researcher 2bE[K]/T ≳ E[K]).
- **C-R3 (seed-and-harvest).** At b=3 discrimination value jumps 0.9 → 39 (b=6: 56). Optimal
  round-1 share 0.107 < 1/6 even; mass deploys from round 2 once posteriors have separated.
  Learning front-loads the observation and back-loads the money. Strictly positive round-1 spend is
  guaranteed by complementarity (the lump-equivalence correction: early money enables, later money
  exploits; at R0=0 the two lump schedules are payoff-identical, so the conjecture licenses positive
  early spend, not early mass).
- **C-R4 (paid forces resupply free forces; QUALITATIVE ONLY).** Two halves of one symmetry: the
  funder's own early grants raise K and thereby regenerate the free-growth condition that rewards
  spending late (why back-loading is over-determined in poverty); and re-deciding each round on
  posteriors updated by one's own funding captures paid information automatically (why deliberate
  deferral loses to even-tranches-plus-re-deciding). B regenerates C; D regenerates E. The free
  forces dominate even when their exogenous sources are switched off.

## 2. What does NOT enter the main text

- **The S8−S5 = −16 / −45 magnitudes.** These are not findings about planning; they are the CE
  planner mispricing information at depth. The paper's own history makes this radioactive: a
  forward planner scoring below myopic was previously the signature of a bug (STATE_OF_PLAY §5),
  and the CE validation (±0.4%, F8) was run at base parameters only. Appendix with the caveat, or
  nothing. The defensible sentence: "no force favors early mass: statics favor even, information
  favors late," supported by S5 beating S8 and by the schedule sweep in item 4.1 below once run.
- **F8 scope correction (must happen).** The methods section's CE claim narrows to: near-exact at
  base parameters (T=2–5, ±0.4%); the CE information term misprices at exploration depth, where
  conclusions rest on the myopic re-deciding benchmark, not on the CE schedule being optimal.
- **A "front-loading regime."** See C-R1 calibration rule.
- **ε* as a sharp constant.** 0.02 spans 0.005–0.05 across r_min; report the band.

## 3. Placement map (Story 1 architecture)

| Piece | Where | Why |
|---|---|---|
| budget_ref decoupling | Setup, one sentence + computational appendix | Already in the errata's Eq (10) replacement |
| C-R2 thin-grants mechanism | The obstacle section (K unobserved), beside the pubs-confound | It is the second way output fails as a signal; sharpens "cost of not knowing the gap" and feeds the peer-review section: in poor fields, buy review or buy depth |
| b ≳ T/2 threshold | Peer-review/policy discussion, one line | The quotable policy quantity for nascent fields |
| Bootstrap subsection (conjecture → lump-equivalence → C-R1 → C-R3) | Extends the timing section (old Part V) | The objection at full strength, then corrected; keeps the timing section honest without re-inflating it |
| Seed-and-harvest schedule figure | The subsection's single figure | exploration_depth_schedules.png; the shape IS the argument |
| C-R4 symmetry (B→C, D→E) | Discussion (the free-dominates-paid close) | Strongest version of Story 2; one paragraph, qualitative, with the CE caveat in a footnote |
| resource_regime heatmap, suite tables, R4 magnitudes + caveat | Appendix | Robustness/completeness; not load-bearing |
| Templeton note | Cover memo to Templeton, not the paper | Field-building is their program; the nascent-fields result is the deliverable's most actionable piece |

Proportion discipline: the timing section stays demoted (its effects are still ~1/30 of the signal
story). The bootstrap subsection earns its space as the objection to it, not as a second headline.
Target: 1.5–2 pages + one figure in the main text; everything else appendix.

## 4. Verification before the integration is final

1. **Honest schedule sweep (the missing certification).** The seed-and-harvest claim currently
   rests on an uncertified CE schedule. Cheap fix, independent of the CE planner: grid over fixed
   two-block schedules (share x of budget in rounds 1..k, remainder after; x × k grid), each
   executed with the myopic within-round allocator under true dynamics, in the depth corner (b=3)
   and one poverty-ε cell. If no early-mass schedule beats even/late-mass, C-R3 is certified
   planner-free. Add to the running suite as a follow-up item; ~minutes.
2. **CE self-consistency check at depth:** does S8 maximize its own CE objective there (honest
   mispricing) or score below S5 on it (bug)? One cell, one diagnostic, distinguishes "model
   limitation to state" from "code to fix."
3. **Seeds 64 → 200** for any figure or number that enters the main text (suite machinery exists;
   the rest of the paper is 200-seed and a referee will ask).
4. **ε* boundary phrasing** per §2.

## 5. Findings register updates

- F9 RESOLVED: no front-loading regime; boundary vertical in ε at ε* ≈ 0.02 (band 0.005–0.05);
  poverty mutes back-loading, never flips it.
- F8 SCOPE NARROWED: CE near-exact at base; misprices information at exploration depth.
- NEW F10: thin grants are talent-uninformative (λ pinned to g when g ≪ K); identification
  threshold g ≳ K, i.e. b ≳ T/2 at the decoupled purse.
- NEW F11: at depth, paid information works and the optimal schedule is seed-and-harvest
  (discrimination value 0.9 → 39 → 56 across b = 0.5, 3, 6; round-1 share 0.107 vs 0.167 even);
  across the suite b_idx never < 0.5 beyond noise.
- NEW F12 (qualitative): paid forces resupply free forces (B regenerates C, D regenerates E); the
  free-force dominance survives switching off its exogenous sources. Carries the CE caveat.

## 6. What waits for the main suite

Package C (ε_free/ε_paid decoupling) now matters MORE, not less: C-R1's "the paid bootstrap
resupplies the free-growth condition" is exactly the B-regenerates-C claim, and the decoupling
sweep is its direct test. If Package C confirms (pure-paid front-loads, pure-free back-loads, and
the coupled model back-loads because paid growth feeds the free channel), F12 upgrades from
qualitative to demonstrated and the discussion's symmetry paragraph gets numbers. Seed-floor
section and Story 5 unchanged, still gated on Packages A and B.
