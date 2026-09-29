# Grant-funding paper: state

Created 2026-08-05. Maintained by Claude. Working title: TBD.
Model repo: `~/Documents/GitHub/Grant-Funding-and-Scientific-Output-Model`.
Coauthors: Simon, Kevin. Templeton report obligation: yes.

## Source of record

Two documents supplied by Aydin 2026-08-05, and they supersede the repo docs:
- `Optimal_Funding_Strategies_for_Scientific_Output.pdf` (16pp) - the Templeton-style report, his
  favourite framing so far. Uses the EXACT allocator and a common normalization (percent of
  no-funding output S1). Repo `RESULTS.md` is stale relative to it.
- `Optimal_Funding_Strategies_for_Scientific_Output_1.pdf` (14pp) - the model write-up. Confirms the
  CES generalization exists in the SPEC (not in `model.R`, which still hardcodes harmonic).

Claude's critique delivered 2026-08-05 as `framing-critique.md`.

## Corrections to the 2026-08-05 session-1 assessment

- The exact-allocator re-read is DONE. F3/F4/F5 stand; the "RE-READ" flags are cleared.
- Back-loading is NOT a headline candidate. On the common normalization, planning is worth at most
  0.7% of no-funding output against 22.3% for the peer-review signal, roughly 30x. Demote to a section.
- The empirical anchoring is stronger than assessed: Fang et al.'s AUC 0.54 is used as an estimate of
  a very noisy tau_K, locating real NIH review above the precision elbow. That is a calibration hook.

## The central analytic result (was not in the repo, not known session 1)

Closed-form optimal grant: `g*_i = c K_i - R_i`, `c = sqrt(gamma/nu) - 1`, nu the Lagrange multiplier;
c set by the budget. Fund the talent-resource gap. Funding frontier `K = R/c`. Verified in simulation
at r = 0.95 to 1.00. Derived under full information; the funder allocates on posterior-expected
marginals, and the gap between the two is exactly where peer-review value comes from.

## Candidate stories (full statements in framing-critique.md)

| # | Story | Verdict |
|---|---|---|
| 1 | Fund the gap. Both intuitive rules wrong for the same reason; value of peer review = cost of not knowing the gap | RECOMMENDED SPINE. Type D primary, B engine, E outer layer. ~70% written. |
| 2 | Free forces dominate paid forces | Best closing move. Belongs as Story 1's discussion, not as a rival. |
| 3 | Back-loading (F1+F2) | Demoted to a section. Smallest effect in the paper. |
| 4 | The precision elbow; review capacity should be allocated by field | Story 1's section 4. Best policy hook. |
| 5 | The concentration dispute is a disagreement about sigma | Highest ceiling, currently a conjecture with two endpoints and one interior probe. Gate on the sigma sweep. |

Proposed thesis sentence: **the value of peer review is the cost of not knowing the talent-resource gap.**

## Open contradictions: RESOLVED 2026-08-30 (items 1-2), one check owed (item 3)

1. Back-loading attribution: RESOLVED. E11 withdrew the textual contradiction (both source passages
   say free C over paid B; the "paid over free" reading was Claude's misreading), and Package C
   settled it empirically (pure free back-loads, b_idx -> 0.680; pure paid front-loads, 0.481).
   [said 2026-08-30] Aydin confirms: paid-over-free is an error.
2. Production-family (gamma) robustness: RESOLVED, the sweep WAS run (Package B, session 3).
   Verified 2026-08-30 against the model repo (branch smooth-allocator): verdict table in
   docs/PAPER_INTEGRATION_HANDOFF_2026-08-06.md section 2.4; data on disk in
   sweep_results/sigma_tierA_gc0/ (20 files), _gcm3/ (20), _leontief/ (6; 50 seeds, n_steps=800,
   kink caveat), sigma_tierB_concentration/ (15; 200 seeds + 400-seed refinement). Standing verdicts:
   P-B1 MIXED (information story survives the whole CES family including Leontief; Cobb-Douglas
   kills planning/back-loading, so timing claims scope to gamma_ces < 0); P-B2/B3 REFUTED
   (concentration budget-conditional, interior max near gamma = -6 at light tails). The write-up's
   "replicated at gamma in {0,-3,-inf}" claim is no longer aspirational, but it is true only with
   the Cobb-Douglas exception named.
3. Novelty of `g* = cK - R`: skeptical-friend pass done (notes-gap-rule-novelty.md: math is
   standard water-filling, cite as precedent; the statement for science funding appears novel on a
   moderate search). [said 2026-08-30] Aydin believes it novel; a deeper literature check is owed
   before the paper leans on the novelty, on his list.

## Decisions made 2026-08-05 (Aydin)

- STORY LOCKED: Story 1 spine, Story 2 as discussion, Story 4 as section 4, Story 3 demoted to a
  section, Story 5 gated on the sigma sweep.
- Dynamics/production robustness checks NOT done yet; frame with an explicit TBD note to Kevin
  (co-author/editor). Resolves W9: the model write-up's "replicated at gamma in {0,-3,-inf}" claim
  is aspirational, not actual; must not survive into the draft as fact.
- Seed-floor framing DISPUTED: Aydin rejects "floors nearly free"; resolution + diagnostics in
  notes-seed-floor.md. Section held until D3/D4 run.
- Gap-rule novelty: skeptical-friend pass done, notes-gap-rule-novelty.md. Math = standard
  water-filling (cite as precedent); statement-for-science-funding = novel as far as a moderate
  search shows.
- Literature engagement (lotteries, concentration, peer-review reliability), notation fix, claim
  modulation: all agreed.

## Sweep handoff issued 2026-08-05

`docs/SWEEP_HANDOFF_2026-08-05.md` in the MODEL repo: Package A (seed D3/D4/D1, new strategies
S10/S11 = S5/S8 + seed, persistent-floor option), Package B (CES switch with ces_gamma naming and
the /2 bit-compatibility convention; Tier A robustness at gamma_ces in {0,-3} + Leontief boundary;
Tier B concentration-vs-sigma gate at T=1 with Gini/coverage metrics and a monotonicity refinement
rule), Package C (decoupled eps_free/eps_paid knowledge growth to attribute back-loading;
preregistered: pure-free back-loads, pure-paid front-loads). All predictions preregistered; results
memo lands as docs/DIAGNOSTICS_RESULTS_2026-08.md in the model repo.

New flags found while writing it: model write-up Eq (10) says budget is not transferable across
rounds, but the implemented forward strategies allocate the whole remaining horizon budget (report
p.8 schedule shares); and the write-up's default table (T=10, n=100, B=50) disagrees with the code
base params (T=2, n=50, b=0.5). Both must be fixed in the model spec before section 5 is drafted.

## Errata issued 2026-08-05

`spec-errata.md` (this folder): 12 items with exact replacement text for the report + model
write-up sources (not in either repo; apply at the LaTeX source). Highlights: canonical notation
(A = productivity, gamma = CES exponent, rho = K-R correlation, alpha_K/alpha_R tails); Eq (10)
rewritten to match the implemented budget (fixed total purse, forward transfers across rounds);
defaults table replaced with the sweep base; TBD-to-Kevin flags on the unrun robustness sweep, the
unimplemented lookahead-h interpolation and baseline options, and the quadratic-in-T claim.
E11: the W4 "textual contradiction" (V.2 vs p.4 margin) was Claude's misreading, WITHDRAWN; both
passages say free C over paid B. The open question is evidential only, tested by Package C.
K-terminology (talent vs knowledge) is Aydin's word-level call; recommendation: talent.

## CAMPAIGN COMPLETE 2026-08-06: canonical results of record

`docs/PAPER_INTEGRATION_HANDOFF_2026-08-06.md` in the MODEL repo (branch smooth-allocator) is
canonical over all earlier memos. 20/20 preregistered verdicts in; every number re-derivable via
`Rscript sweep_results/_probe/verify_all_claims.R` (not runnable from cloud/VM sessions: no R;
spot-run owed from a Mac session). Verdict summary:

CONFIRMED: P-A1 (floor on strong optimizer costs 6.11% of S1 at heavy+sharp), P-A2 (exact),
P-A3 (convex), P-C1/C2/C3 (free back-loads 0.680 / paid front-loads 0.481 / diagonal identical),
P-D1 (implied tau_K > 20 at AUC 0.54), P-D3 (resource signal redundant <=1.1%), P-D4 (gap-rule
convergence 0.34 -> 0.97, oracle gap 17.2 -> 1.7% of S1), P-E1a (scaled rule r 0.87-0.98),
P-E2 (n/b robust), Boot-1/2/3 (seed-and-harvest planner-free; negative PG = honest CE mispricing
20/20; 64-seed stable at 200).
MIXED: P-B1 (info story survives whole CES family; Cobb-Douglas KILLS planning: scope timing
claims to gamma_ces < 0), P-D2 (overtrust negative only at moderate tails; heavy tails attenuate).
REFUTED: P-B2/B3 (concentration NOT monotone in sigma: budget-conditional + interior max at
gamma ~ -6; Story 5 dead as stated, replaced by the budget-conditional law), P-E1b (product
reduction fails materially: 1-D talent is a substantive commitment; A and K identified only by
dose-response).
SCOPED: A5 (exponent 2.1-2.3 for T<=5; 0.8-1.4 over T=5-10: decelerating, say neither
"quadratic" unqualified nor "no saturation").

Held decisions RESOLVED: (1) seed section = the price law, "floors cost what targeting is worth
times the floored share," lotteries engaged where the price is lowest; (2) Story 5 dead, F13
budget-conditional concentration law + Cobb-Douglas boundary become the robustness section's two
real results; (3) back-loading mechanism text evidenced (free channel dominates above
eps_free/eps_paid ~ 1/8-1/3), and the bootstrap conjecture's kernel gets its precise home: pure
front-loading is real exactly when knowledge grows only through funded work (the
fields-that-don't-exist-yet corner).

Claude's two additional flags beyond the handoff's nine: (i) long-T terminology: "decelerates,"
retire both "no saturation" and unqualified "quadratic"; (ii) state "no front-loading (coupled
model)" and "pure-paid front-loads (decomposition)" together with the F2 threshold, else they
read as a contradiction.

Unified discussion spine now available: one quantity, the dispersion of posterior-expected
marginal returns (= the value of knowing the gap), prices peer review, targeting, floors, and,
via complementarity, scopes the timing story.

NEXT: the outline (workflow step 2), unblocked, on Aydin's word. Then intro (step 4).

## Bootstrap/exploration suite integrated 2026-08-05

Aydin ran his own 64-seed suite (resource_regime, exploration_corner/poverty/depth; write-up
T_round_extension/RESOURCE_REGIME_RESULTS.md). Integration memo: notes-bootstrap-integration.md
(this folder). Register updates: F9 RESOLVED (no front-loading regime anywhere; boundary vertical
at eps* ~ 0.02, band 0.005-0.05; poverty mutes back-loading, never flips it); F8 SCOPE NARROWED
(CE near-exact at base only; misprices information at exploration depth); NEW F10 (thin grants
talent-uninformative, threshold g >~ K i.e. b >~ T/2); NEW F11 (seed-and-harvest at depth,
discrimination 0.9 -> 39 -> 56, round-1 share 0.107 vs 0.167 even); NEW F12 qualitative (paid
forces resupply free forces, B->C and D->E; CE caveat attached).
DO NOT put in main text: the S8-S5 = -16/-45 magnitudes (CE mispricing, appendix with caveat);
any "front-loading regime" claim (b_idx 0.4985 is even-split); eps* as a sharp constant.
Verification owed before final: honest fixed-schedule sweep at depth (certifies seed-and-harvest
planner-free); CE self-consistency check (mispricing vs bug); 64 -> 200 seeds for main-text
figures. Package C of the main suite now doubles as F12's direct test (B-regenerates-C).
Main diagnostic suite (Packages A/B/C) running; results memo pending.

## Addendum handoff issued 2026-08-05

`docs/SWEEP_HANDOFF_ADDENDUM_2026-08-05.md` in the MODEL repo. Package D (information integrity,
all promised by the current text): D-1 model-implied AUC calibration inverting Fang et al.'s 0.54
to a tau_k (makes footnote 13 quantitative); D-2 misspecified trust (tau_k_true vs tau_k_belief;
the peer-review caveat curve; predicted asymmetry: overtrust can go negative, undertrust only
attenuates); D-3 resource-signal ablation (report Appendix A's own flagged unrun sweep; config
only); D-4 allocation convergence to the gap rule (corr(g, cK-R) vs tau_k plus an oracle gap; the
thesis figure). Package E (insurance): E-1 heterogeneous productivity, closed form survives as
g* = c_i K_i - R_i with c_i = sqrt(2A_i/nu)-1, sweep confirms; E-2 headline surface at n=200 and
b x T at high epsilon. Package F (GATED, Aydin's go required): grant persistence phi, the one
structural perturbation under which back-loading could genuinely reverse; registered as uncertain;
decline-path = scope the claim in setup + limitations.

## Decisions made 2026-08-05 (evening)

- Package F CANCELLED by Aydin, with rationale recorded in the addendum: grants do not buy durable
  capacity outside endowments; the durable residue of a grant (skills, publications, reputation)
  IS K, so the paid-growth channel already carries it. Consumable-grant assumption goes in the
  setup as a substantive claim; limitations gets one line (endowments; capital-intensive/Leontief
  fields, folding into the substitutability discussion).
- E-1 redesigned around Aydin's product-reduction hypothesis (variation in c folds into the
  distribution of cK): not exact algebraically (ceiling 2AK vs saturation scale K), so E-1b tests
  whether it holds approximately: same log-variance of T = AK split between latent factors,
  s in {1, 0.7, 0.4}, readouts output / grant Gini / signal value / who-gets-funded. Approximate
  reduction keeps 1-D talent in the paper; failure goes to the appendix. E-1a (scaled rule
  g* = c_i K_i - R_i, A observable) unchanged.
- Empirical note banked for the discussion: A and K are separately identified only by dose-response
  (output vs funding depth), never by cross-sectional output.

## Session 3 (ran 2026-08-14; files named 2026-08-06 after the context-transfer date): sweeps landed, command built

All diagnostics complete: 20/20 preregistered verdicts in
`docs/PAPER_INTEGRATION_HANDOFF_2026-08-06.md` (model repo, branch smooth-allocator; CANONICAL,
with `DIAGNOSTICS_RESULTS_2026-08.md`). Read this session together with both PDFs and all notes
files. Register consequences:

- D-4 CONFIRMED: the thesis figure exists (corr 0.34->0.97 heavy tail; oracle gap 17.2->1.7% of S1).
- D-1 CONFIRMED: implied tau_K at AUC 0.54 is >20 (heavy/moderate tails), elbow ~2.5. Footnote 13
  survives with a wide margin.
- D-2 MIXED: overtrust asymmetry confirmed; the negative regime is tail-dependent (negative only at
  alpha_K=2; at 1.3 it only attenuates, 13.99 vs 17.20).
- D-3 CONFIRMED: resource signal redundant (<=1.1% of S1 everywhere). Closes report Appendix A.
- Package A: floors are NOT nearly free on a strong optimizer: 6.11% of S1 at the focal cell; price
  law cost ~ kappa x targeting value x floored share, kappa in [0.10, 0.86]; convex in x_seed.
  Seed-floor section UNFROZEN (E12 resolved); write as price schedule.
- Package B: P-B1 MIXED: information story survives whole CES family incl. Leontief; Cobb-Douglas
  KILLS planning/back-loading, scope Part V to gamma_ces<0. P-B2/B3 REFUTED: concentration is
  budget-conditional (tight concentrates, ample equalizes toward Leontief; interior max at light
  tails). Story 5 dies as stated; replacement law is a headline robustness result.
- Package C: attribution resolved. Pure free back-loads (b_idx->0.680); pure paid front-loads
  (0.481, the only strict-and-real front-loading); free dominates once eps_free/eps_paid > ~1/8-1/3.
  F12 upgraded from qualitative to demonstrated (B regenerates C confirmed).
- A5: T-exponent 2.05-2.33 for T<=5; 0.78-1.44 over T=5-10 (saturation). E10 resolved: cite the
  exponent, scope quadratic to T<=5.
- E-1a CONFIRMED (scaled rule survives, r 0.87-0.98); E-1b REFUTED: product reduction fails
  (signal value -69 to -83%); 1-D talent is a substantive commitment, say so in limitations.
- E-2: n=50 insured (corr 0.996/0.999 vs n=200).
- Bootstrap verifications done: seed-and-harvest planner-free (honest schedules); negative PG at
  depth is honest CE mispricing (20/20), caveat wherever S8-S5<0 appears; 200-seed provenance for
  all main-text exploration figures (use exploration_200/ and resource_regime/, not the 64-seed
  originals).

Do-not-claim additions from the handoff (flag 9): "no front-loading" means min b_idx 0.497
(<=0.003 below 0.5, worthless there); Leontief cells carry the kink-degeneracy caveat.

Session 3 decisions (Aydin, 2026-08-06): VENUE = Research Policy. K TERMINOLOGY = "epistemic
capability" ("capability" after first use; the capability-resource gap; symbol stays K). Aydin
supplied the term himself (rejected talent/knowledge; chose from his candidates research capacity
vs epistemic capability on Claude's ranking: capacity conflates with capacity-building usage that
includes R). Errata E3 applies with "epistemic capability" substituted.

Outline v1 delivered: drafts/grant-funding/outline-2026-08-06.md. Eleven sections, Type D spine
with the lottery misprescription as the proposed section-8 closing move (Aydin has not yet
approved the outline or the misprescription move). E9 superseded: Package B ran, so the
robustness sentence becomes a true claim with the Cobb-Douglas exception named.

Outline v1 APPROVED (Aydin, 2026-08-06): adversarial lottery engagement in (he holds the
pro-lottery prescription wrong and potentially harmful where it matters most); E9 supersession
in; no intro vignette; TITLE = "Funding the capability-resource gap: optimal science funding
across multiple rounds".

ABSTRACT LOCKED (2026-08-06, verbatim in intro-draft-v1.md): 100 words, no numbers, no
two-rules hook; thesis sentence retained; closes "Whom to fund matters far more than when to
fund." Abstract rules learned: no numbers, no puzzle, <100 words target, qualitative results
stated plainly.

Intro draft v1 delivered (intro-draft-v1.md, ~1,060 words, 10 paragraphs): awaiting Aydin's
critique. Known open items in it: lottery-funder citations are placeholders [SNSF, NZ HRC,
Volkswagen; check]; the "part of the case for lotteries leans on review's predictive record"
warrant needs a citation check; the 17-to-2-percent and one-thirtieth numbers carry
at-our-settings scope.

Intro v2 delivered (intro-draft-v2.md) after Aydin's style pass on v1. Style rules learned
(recurring; candidates for the paper-writing skill): no epigrams/riddle-compression; no
philosopher's-register tags ("on our account", "hereafter", "locates", "present themselves");
enact, never announce; he supplied three verbatim sentences (opener "Science funding is large,
critical, and its optimal allocation is still not well understood."; "the decision is this:";
"We model this decision situation as follows."). NEW ERRATUM E13 (report): the Part I funding
sentence mixes current-dollar total ($993B, 2024) with constant-2017-dollar funder amounts
($592B/$148B); replace with current-dollar throughout ($937B 2023; $993B est. 2024; business
75%/$743B; federal obligations $194B FY2024; NSF 26-314 + NSB S&E Indicators 2026).

Intro v3 delivered (intro-draft-v3.md) after Aydin's line-edit pass on v2; his lessons
extracted to style-lessons-2026-08-06.md (13 rules + consistency flags; skill-gaps candidates).
Key new rules: motivation experiential (reader in funder's shoes, second person); jargon only
as working terms of art; spell out compressed references; method labels wait for results
sections; NEVER claim novelty explicitly; no result numbers in intro; results as method-neutral
claims ("We show...", never "in simulation" before methods); never foreclose the paper's arc;
capacity over identity claims ("review can provide"); clarity beats compression. Evidence
recalibration: Fang = "barely better than chance", funded-grants scope noted, positive-validity
studies cited alongside (Li & Agha 2015; Park, Lee, Kim 2015; Gallo 2014; verified via
secondary survey, originals to check). PENDING AYDIN: terminology "uniform seed grants" (his
edit) conflicts with locked abstract's "uniform funding floors"; proposed abstract emendation
in v3 notes. "Options" adopted for the two intuitive allocations.

Intro v3 timing/roadmap emendations applied (Aydin's parallel timing sentence; roadmap "We
proceed as follows." + \S). PAPER-WRITING SKILL UPDATED: core-style.md addendum (evidence
calibrated against own machinery; capacity over identity; clarity over compression; parallel
contrast sentences; anti-patterns 19-21: no epigrams, no philosopher's register, no novelty
claims) + section-patterns.md 03a addendum (motivation register rules; roadmap standing form);
repackaged skills/paper-writing.skill in repo and delivered; Aydin must save it.

TEX STARTED: tex/main.tex (canonical amsart preamble; abstract + intro; compiles clean,
5pp with bibliography) + tex/references.bib (21 entries, CHECK notes on unverified fields)
+ main-preview.pdf. Container lacks chicago.bst so preview used plainnat; source keeps
\bibliographystyle{chicago} (present in MacTeX). TODOs in tex: Simon's surname/affiliation;
abstract emendation "uniform seed grants" pending Aydin sign-off; \S3-\S11 hard-coded until
sections exist; lottery-funder cites.

Section 2 draft v1 delivered (section-2-draft-v1.md): four-literature funnel (concentration,
review predictivity, lotteries, formal models), show-not-say positioning, water-filling +
poverty-gap precedents placed here per Aydin's intro edict, Azoulay incentives paragraph
optional. All [check] citations need original-source verification before lock.

FRAMING DECISION (Aydin, 2026-08-14): the paper is NOT framed as settling the peer-review
dispute. Paragraph 1 reframed to the question family: disputes (concentration, review
predictivity, lotteries) unified as one problem (fixed budget, imperfect knowledge of
researchers); contribution = working out the logic of optimal funding (whom, when, what
information is worth, what egalitarian alternatives cost). Thesis sentence stays at full
strength in the review paragraph as the central finding, not the advertised purpose.
Concentration fact returned to P1 as a member of the question family (still also in \S2).
Funding sentence in one accounting frame ("roughly a fifth of it funded by the federal
government"); "obligated"/"fiscal" banned as jargon; $194B obligations figure dropped
(frame-mixing). Author block complete: Simon DeDeo, Social and Decision Sciences, CMU.
Running header short title: "Funding the capability-resource gap".

Section 2 revised per Aydin's feedback (2026-08-14): opening = "Our contribution sits at the
intersection of four literatures."; formal models characterized precisely (Gross-Bergstrom
contest model verified, PLoS Biology 17(1):e3000065; Avin BJPS 70(3):629-656 landscape
simulation verified, now Avin 2019b beside Mavericks 2019a; Park-Lee-Kim Research Policy
44(6):1145-1159 verified); "makes both latent" replaced with the bottleneck formulation;
"central objects" sentence replaced with the three-declarative chain (rule, records fail,
review's value). Bib updated accordingly.

Intro COMPRESSED (2026-08-14, calibrated against Aydin's three published intros: HARKing,
Causation, Social Sciences Hard): ~1,050 -> ~650 words; five result paragraphs to two;
mechanisms relocated (thin grants -> S5; price-law structure -> S8; scoping detail -> S9);
model paragraph to five sentences; commitments to three. Lesson 15 added to style-lessons
(findings at sentence strength; ~600-700 word target; not yet in the .skill package, batch
at next repackage). Azoulay/incentives paragraph moved from S2 to NEW section-10-notes.md
(discussion ledger, includes previously banked discussion items).

LINE-EDIT PASS APPLIED (2026-08-14, lesson 16 in style-lessons): "reads" banned (use
"accounts for"); assumptions not commitments; "our model" not "the model"; NO formulas in the
intro (gap rule now described in words; g* = cK - R debuts in \S4); "technology" cut
(never introduced); "Peer review supplies an additional signal of capability"; his verbatim
lines applied (P1 "In spite of this...", "Whether peer review is effective is disputed",
"Each option can be intuitively compelling", "captures most of this value", "The same
reasoning elucidates", "We show how strategic timing...").

MAJOR DECISION (Aydin, 2026-08-14): the capability-is-one-dimensional thread is DROPPED from
the entire paper (assumptions paragraph now two assumptions; E-1b heterogeneity appendix OUT
of the paper plan; outline's limitation item superseded). Claude's flag, recorded once:
E-1b is a real fragility (signal value falls 69-83% when latent productivity takes talent
variance); with the thread dropped the paper is silent on it; recommended fallback if Kevin
objects: one limitation sentence in \S10, no appendix. Aydin's rationale: distracts from
the clarity of the results.

TIMING SENTENCE CALIBRATION (2026-08-14): Aydin's proposed "spending later or earlier can be
more effective as a function of one's budget relative to researchers' total resources" was
adjusted to "how much of the budget to spend early rather than late shifts with the funder's
budget relative to researchers' total resources": the do-not-claim list bars any
front-loading-regime reading (no strict-and-valuable front-loading exists in the model); the
dependence is in the schedule's shape and strength, not its direction.

SECTION 2 v2 (2026-08-14, reflection pass approved and applied, section-2-draft-v2.md; v1
renamed -superseded): lotteries paragraph rebuilt around the cost-side/value-side complement
(Gross-Bergstrom price what a lottery saves, effort; we price what it forgoes, targeting;
moved G&B there from the models paragraph); timing-coverage sentence added (none treats
timing of spending as the funder's choice; [check] Avin dynamics); Sorzano characterized
precisely; warrant softenings ("orders of magnitude" out; "most often pressed" out);
precedents moved to NEW section-4-notes.md (water-filling, poverty-gap, plus banked \S4
items: frontier term, formula debut, dissociation refutation, figures); closing chain
replaced with synthesis close (the four literatures meet in the gap); "our model defines the
quantity at issue"; "long-running dispute" aligned in intro and \S2; "a single quantity"
replaces "one-dimensional" for rival models.

SECOND \S2 PASS (2026-08-14): "efficacy of peer review" is now the term for the disputed
property (funnel + topic sentence; study descriptions keep factual verbs). The
concentration-positioning sentence ("Our model treats the two questions as one") FAILED
Aydin's read: no two questions had been posed, and the unification described nothing the
model does. Root cause named: form-first writing (rhetorical schema chosen, content
back-filled, driven by the per-paragraph pivot template). Rewritten around the true claim
(concentration answer depends on budget relative to researchers' total resources; tight
concentrates, ample spreads). Lesson 17 (forward-only reading test; form never precedes
content) added to style-lessons AND to the .skill package (core-style addendum; repackaged,
third revision today, Aydin must re-save). Forward-only sweep of intro + \S2 completed.

ORDER OF VIRTUES codified (2026-08-14): Aydin's clarity-above-all statement is now hot rule
0 of the paper-writing skill and a core-style section, with two enforcement procedures:
plain-first drafting (generation) and the paraphrase test (audit); the skill's own
punch-rewarding rules (bold closes, verdict resets) explicitly subordinated. His statement
recorded verbatim in worldview.md section 5 with [said] tag. Skill repackaged (fourth
revision today; the newest supersedes all). Lesson 18 in style-lessons.

Session 3 next: Aydin's verdict on \S2 v2 second pass, then model section (\S3), to be
drafted plain-first. Context-transfer copy added this session.

GAP-RULE PROOF (2026-08-14, session 3 continued): Aydin's five \S4 line edits applied
(sweep sentence in main.tex now "supported by proofs / established across"; "We suppose
first"; the "Within a round" footnote with his exact wording; "who can produce more value
with it"; objective display confirmed on its own line). Then, per his directive to derive
rather than assume the gap rule: NEW gap-rule-proof-v1.md, self-contained and
appendix-ready. Structure: meticulous setup (A1: primitives; A2: assumptions, with the
note that nothing beyond A>0, K_i>0, R_{i0}>=0, B>0 is used; A3: problem (P)); Lemma 1
(return to a grant: phi_i' formula, strict concavity, saturation at 2AK_i); Lemma 2
(existence, uniqueness by strict concavity + compactness, budget exhaustion); Lemma 3
(equal marginal value, iff: necessity by transfer argument = the \S4 sentence made
formal, sufficiency by concavity gradient inequality; NO KKT machinery, flagged as note
1); Proposition 1 in four steps (nu in (0,2A) so c>0, NOT automatic, proved from
existence of a funded researcher; funded inversion; unfunded complementary slackness;
c pinned uniquely by the budget identity via G(c) strictly increasing where positive);
Corollary 1 (frontier K_i > R_{i0}/c; c(B) continuous strictly increasing; nu(B) =
2A/(1+c)^2 decreasing); Remarks (water-filling with capability-scaled fill level; FGT as
the K_i = K special case; scope across rounds). Boundary detail surfaced: R_{i0} = 0
researchers are funded at every budget. NUMERICALLY VERIFIED (verify_gap_rule.py): 81
trials (n=50; tails 1.3/2/3.5; b 0.1/0.5/1; A 0.5/1/2; R0=0 boundary exercised) against
scipy SLSQP: formula never below the optimizer, allocations agree to 6.5e-7 of budget,
marginal-value conditions hold to machine precision, funded counts 2..50 (frontier
active), c(B) strictly increasing spot-checked. Awaiting Aydin's verdict on the proof
document; then LaTeX conversion (amsart theorem environments) and \S4 finalization.

\S4 NUMBERS-OUT REVISION (2026-08-14, Aydin's direction): \S4 is now qualitative only;
every claim in it is derivable from the gap rule. His topic sentence applied ("The gap
rule explains why each of the intuitive funding options discussed earlier fails";
options pluralized, flagged). Simulation numbers (23.5/26.1, 30.1/29.8, 12->5 percent)
and report Fig 2 banked in section-4-notes.md for the quantitative sections. Both
refutations now analytic: proof doc gains Example 1 (track-record-proportional loses to
uniform: K=(8,4), R0=(8,1), B=2; outputs 10.75 < 11.14 < 11.43, exact rationals
570/53 < 568/51 < 80/7) and Example 2 (all-to-scarcest loses to uniform: K=(10,1/2),
R0=(5,0), B=2; 112/15 < 49/6 < 42/5), both verified in exact arithmetic. Budget static
now analytic: NEW Corollary 2 (optimal-minus-uniform output difference vanishes as
B -> infinity, absolutely and relative to uniform's gain; proof via the capability
ceiling 2A*sum(K)). ERROR CAUGHT in the process: old \S4 sentence "the optimal
allocation converges toward uniform funding" is false (grant shares converge to
capability-proportional); replaced by the true outputs-converge claim. Closing replaced
with Aydin's sentence, one flagged alignment (funds -> resources, order matched to
capability-resource gap); three second-sentence candidates ranked in draft notes, (a)
roadmap version currently in text. Figure 1 relabeled analytic illustration (computed
from Proposition 1, no seed provenance needed). Awaiting his verdict on both documents.

\S4 CLARITY PASS (2026-09-05, Aydin's line edits): "Consider the marginal value of a
dollar" (guiding over bombastic; skill anti-pattern 23). Proposition 1 restated to
preempt the every-researcher misreading: "fills the shortfall, if any, between the
resources they hold and a target proportional to their capability"; interpretation
paragraph rebuilt around the target mechanism (c sets targets cK_i; all positive
shortfalls funded in full simultaneously; budget disciplines targets; frontier defined
plainly). Water-filling + FGT precedents and both appendix-example pointers moved to
footnotes (his body/footnote rule; skill core-style block + lesson 21). TERMINOLOGY
CHANGE: "schemes" replaces "options" for the two intuitive allocations, propagated to
the intro in main.tex (three spots). His capacity sentence applied ("who is at capacity",
typo fixed from "whose"); flagged: capacity/capability collision, alternative offered.
His closing "determined by their difference" NOT applied as written, flagged for truth
(marginal value is determined by the ratio R/K, not the difference; drafted as "how a
researcher's resources compare with their capability"). Budget paragraph rewritten
(ceiling first, converge under ample budget, reverse under tight). "Dollars allocated"
per his edit. Preview rebuilt. Awaiting verdicts on the flags (draft note 4).

\S4 SECOND CLARITY PASS (2026-09-05): his four edits applied ("First, consider the simple
case where..." replacing "We suppose first," because complete information is an idealized
case considered, not an assumption made; "Now, consider the marginal value of a dollar
for such a funder."; "exactly" cut from the budget clause; "We can understand the rule as
follows."). His largest-gap-first confusion resolved by rebuilding the interpretation
paragraph around the rising-c picture: targets are endogenous to the budget (c stops
rising when grants sum to the budget), so the budget never falls short of the targets;
the first dollars do go to the researchers with the least resources relative to their
capability, which is the correct form of his intuition (priority is by the ratio R/K,
not the absolute difference). Water-filling footnote unchanged; frontier sentence kept.

\S4 FIGURATIVE-LANGUAGE PURGE (2026-09-05, after Aydin's ceiling correction): the run of
failures he caught ("disciplines," "picture c rising," "stops c early," "ceiling,"
"saturates," and earlier "prices both at once") named as one pattern, figurative language
imported into the description of formal objects, now skill anti-pattern 24 + a core-style
block with a substitution test (every sentence about a formal quantity must translate a
statement of the formalism, no surplus) + lesson 22. THE TRUTH POINT: expected output
under harmonic Lambda increases with every dollar (diminishing returns) and is bounded
above by 2AK_i, approached but never attained; "ceiling"/"saturates" wrongly suggested an
attained plateau. Fixed in \S4 (budget paragraph now states the bound precisely), in the
proof doc (Lemma 1 commentary, Corollary 2 commentary, "pinned" and "loses to" removed);
Corollary 2 itself unaffected (its proof uses only the bound and the limit, both true).
Interpretation paragraph rewritten statically: gap defined as target minus resources
(signed), grant = positive gap, c defined as the value at which grants sum to the budget,
"the budget always covers the gaps: a question of which gaps to fund first does not
arise"; small/large budget comparative static in increases/decreases vocabulary.
Proposition 1 now says "fills the gap, if any" (shortfall dropped; one concept one term).
His within-a-round footnote moved into the body as his sentence ("Further, we initially
restrict our attention to a single round of funding..."). Equal-marginal sentence
completed (unfunded condition added). "Increases/decreases" standardized: abstract
("signal's value increases with"; NB abstract was locked, changed under his "standardly
throughout"), \S2 tex+md ("increases or decreases as grants concentrate"). "At capacity"
replaced with "whose resources approach or exceed their target" (his ceiling correction
applies to it; flagged for veto). "Returns twice" -> "recurs"; "drives" -> "determines."
Remaining "grows" flagged for his call: \S3 "Capability grows through research" (his
approved text). Skill repackaged (sixth revision; supersedes all).

\S4 INTO MAIN.TEX (2026-09-05): Aydin approved figure option iii (relative value of
targeting). Caption cleaned per his edits (Proposition-1 clause cut; Corollary 2 kept
parenthetical; no redundant final restatement). \S4 typeset and inserted as
section 4 (\label{sec:optimal}, prop environment with footnoted c and Appendix [X]
placeholders matching \S3's convention; objective written with A*Lambda notation;
water-filling/FGT and appendix-example footnotes; frontier-map figure as TODO comment;
Figure~\ref{fig:targeting-value} cited in the budget paragraph and \input after it).
pgfplots added to preamble (compat 1.17). Wording applied earlier this day: "Capabilities
increase through research" (main.tex + \S3 draft), "including grants" (\S4). Compiles
clean: 12pp, 0 errors, no undefined references. Next: his verdict on the typeset \S4,
then \S5 (track-record estimation; thin-grants ledger staged); appendix conversion of
gap-rule-proof-v1.md still owed (will resolve the Appendix [X] and Corollary 2
references).

FIGURE 1 (FRONTIER) + TABLE 1 FIX (2026-09-05): Table 1 now footnotesize on a tabular*
spanning exactly \textwidth (extracolsep fill; no line breaks, no overflow). NEW
fig1-frontier.tex + fig1-funded.dat + fig1-unfunded.dat: the gap rule in one population,
(R, K) plane, reduced from the GUI trio to the minimal Proposition-1 content (one round,
complete information, optimal grants; no bottleneck coloring, no correlation
annotations): same seed-7 population as Figure 2, b = 0.25 (25/50 funded, the clean
split; b stated in caption per the figure-defaults convention), c = 1.099; frontier
R = cK drawn and labeled, funded filled with horizontal quiver arrows ending exactly on
the frontier (arrow length = gap), unfunded open. Cited at the frontier sentence.
Compiles 13pp clean; figure pages inspected (frontier line trimmed to plot area,
"carries" -> "moves" in caption). Numbering: frontier = Fig 1, targeting value = Fig 2.

PAPER-FIGURES SKILL ARRIVED mid-task (Aydin saved it; synced after both figures were
built). Conformance check run: spec/takeaway-first captions, slot sizing (native pgfplots
inherits document font at true size), Tufte rules, direct labels, lint (figlint on both
figure pages: no label collisions in the figures; sub-7pt glyphs are LaTeX math
scriptscript sizes in footnotes/captions, standard, not figure text), page inspection
done. DEVIATIONS FLAGGED for Aydin: (1) tool choice: skill prescribes ggplot2 via
tikzDevice for data plots, pgfplots only for TikZ-annotated plots; Fig 1 qualifies
(quiver annotations), Fig 2 is a pure data plot in pgfplots. R is unavailable in cloud
sessions, and the two figures match each other; decision owed: keep analytic figures in
pgfplots and match the ggplot2 style when simulation figures arrive, or regenerate Fig 2
in ggplot2 from a Mac session. (2) Figures are \input tikzpictures, not standalone-
compiled PDFs per the skill's kit; works, but the kit (figures/figstyle.tex, FIGURES.md)
should be adopted before the simulation figures. Minimal FIGURES.md started in tex/.

FIGURE POPULATION REVISED (2026-09-05, Aydin's notes on Fig 1): tool choice settled,
pgfplots throughout with matched styles [said]. New population found by SEED SEARCH
rather than hand-editing points, so the drawn-from-the-defaults caption stays true:
n = 75 (his "a few more points"), seed 659, b = 0.1, c = 0.833 (frontier slope 1.20 vs
prior 0.91; a steeper slope is forced by the 1/5 coverage under iid Paretos, this is the
gentlest that clears it), 16/75 funded (~1/5, his target), no high-K high-R point (his
outlier objection; search excluded K>4 & R>4.5), three over-resourced low-capability
researchers, two funded researchers with gaps > 2 (long arrows: the few-big-gaps image
he wants carried to \S8). Fig 2 regenerated on the SAME population (family invariant
kept; curve verified monotone, 163% at b=0.01 to 13.5% at b=2; y-axis rescaled to 175).
Captions and data script updated (n=75, b=0.1); FIGURES.md updated. Compiles clean;
both figure pages inspected.

FIG 1 SCHEMATIC + \S5 DRAFTED (2026-09-05): Fig 1 per Aydin's notes: leader line to the
frontier removed, axis numbers removed (values do not matter), conventional arrow tips on
both axes, same-population cross-reference cut from the caption; compiled and inspected.
NEW section-5-draft-v1.md (plain-first): P1 cross-sectional underdetermination made exact
(iso-output pairs; the record bounds capability below at lambda/2A but cannot separate;
the two researchers it cannot separate are the ones the gap rule treats most
differently); P2 the measured cost (records-only allocation correlation stays 0.13-0.18,
output 17.2% of S1 short of complete information; the span framing per the do-not-claim
list); P3 thin grants (analytic 2AR limit, K=1 vs K=100 example, depth requirement;
exploration-corner numbers kept qualitative with \S7 pointer, option flagged); P4 bridge
to \S6 (review as the missing capability signal, lesson-10 formulation). Eight notes
incl. the strategy-mapping verification flag (pubs-only vs records+resource-signal).
Awaiting his verdict.

\S5 REVISION + THIN-GRANTS FIGURE (2026-09-06, per Aydin): P2 qualitative in body
(weakly correlated / falls well short), numbers footnoted; telegraph sentence added (the
records-only baseline lives in \S6's convergence figure; \S5 does not duplicate it,
per Claude's recommendation Aydin accepted by directing the qualitative form); bridge
now "under the right conditions. \S6 measures what such a signal is worth, and states
the conditions." NEW PROPOSED fig3-thin-grants (analytic, schematic axes): output vs
resources for K = 3, 10, 30, all slope 2A at origin, coinciding when thin, separating at
depth; vindicates P3 analytically without simulation. FIGURES.md row added (PROPOSED).
Working home now the model repo's for-claude/ (relocated 2026-09-05); this session
edited via bridge staging after a container reset.

\S5 TYPESET + \S6 DRAFTED (2026-09-06): Aydin approved \S5 and the thin-grants figure
("looks great"). \S5 inserted into main.tex as section 5 (\label{sec:records};
Figure 3 = fig3-thin-grants \input after the thin-grants paragraph; the P2 numbers
footnote carries two TODO comments: verify_all_claims re-derivation and the
strategy-mapping check; "increasing precision" -> "increasing informativeness" for the
locked term). Compiles 15pp clean, \S5 pages inspected. FIGURES.md: fig3 status ->
placed. NEW section-6-draft-v1.md (function analysis + plain-first draft): definition
paragraph with the licensed value formulation + output measure; convergence result with
footnoted D-4 numbers (0.22->0.95 defaults, 0.34->0.97 heavy, 17.2->1.7 span) and
FIGURE 4 placeholder (design sketch in note 4, x = informativeness, y = allocation
correlation vs output-shortfall alternative, records-only flat baseline; needs .rds ->
.dat export from a Mac session); determinants paragraph (dispersion, budget inherited
from \S4, informativeness flattening with elbow footnote tau ~ 2.5, sweep 0.05-20);
efficacy-dispute payoff (AUC 0.54 -> tau > 20 footnote; no single number settles it);
overtrust caveat (negative at alpha 2, reduced at 1.3; z-honesty in footnote TODO);
bridge to \S7 (whom -> when). Seven notes incl. sharp/coarse flag and the
gathered-conditions option. Awaiting his verdict on \S6 and the Figure 4 design.

STANDING RULE FROM AYDIN (2026-09-06, skill-gaps observation to promote on recurrence or
his word): re-introduce the paper's central objects of reference at section boundaries,
within reason; readers skim or enter partway. Applied: \S5 opening "evidence about the
capability-resource gap" (was "the gap"); \S6 definition sentence "uncertainty regarding
the capability-resource gap." Named rules (the gap rule) keep their names.

\S6 v2 (2026-09-06, Aydin's restructure): convergence DEMOTED (his verdict: near-trivial
direction, worth mentioning only to explain the logic; "the central result is
convergence" was also unclear). MAIN RESULT recentered on the field delineation, led by
the heavy-tailed case with his two mechanisms (most addable output sits with the few
most capable; they are easiest to identify, standing far apart, so even a noisy signal
separates them) and his corollary at sentence strength (fine distinctions among
borderline applications add comparatively little; effort spent there may be
misallocated; honesty flag: uniform-tau model licenses this via the flattening +
whale mechanism, rank-local discernment not separately simulated). SECOND STANDING RULE
from Aydin: minimize cross-references by number; restate borrowed facts so sections are
self-contained (readers do not have section numbers memorized); applied throughout \S6
v2 (back-refs removed, "The next section" for the adjacent forward pointer); typeset
sections' forward pointers stand pending his word. "Long-tailed" mapped to locked
"heavy-tailed." Both candidate rules (object re-introduction; self-containedness)
logged for skill promotion on his word.

CROSS-REFERENCE RULE REFINED (2026-09-06, Aydin): forward references are often
acceptable (they promise "it will come in \S y"); backward references by number are not
(they demand memory of a number); for the immediately next section, prefer "Next we..."
over "the next section" or "\S n" (concise, unambiguous, no redundancy). Applied
outside the intro: \S3's close ("Next we solve the case..."), \S4's close ("Next we
turn to how such estimation can work: first from the track record alone, then with peer
review added"), \S5's telegraph and close ("Next we ask..." / "Next we measure..."),
\S6 v2 bridge ("Next we turn to when..."). Non-adjacent forward pointers ((\S8),
(\S9), (\S7), and \S2's closing trio) kept per the rule. No backward-by-number
references existed in \S2-\S5. Intro roadmap untouched per his exclusion. 15pp clean.

\S6 v3 + CLAIMS LEDGER (2026-09-06): his line edits applied ("We can quantify this
value in terms of output"; "Two further effects are less obvious" with enumeration,
replacing the spatial "beneath"; "Review can inform whom to fund"; his bridge "The next
section turns to when to fund. In particular, we examine..."). HIS SANITY CHECK on
never-catches-up ANSWERED AND CLAIM CORRECTED: his intuition (resources grow, output ->
2AK, K identified) fails against the spec because resources do NOT accumulate while
capability compounds; output converges to the resource-limited 2AR and carries
vanishing information about K (Fisher per round 7.9e-2 -> 6.1e-7 over 200 rounds;
cumulative info finite at 0.44 over 20k rounds); "never" replaced with the scoped
"does not catch up" plus the structural cause in the body. NEW ANTI-PATTERN 25 in the
skill (SEVENTH revision, packaged; Aydin must save): inflated relational verbs
(tells/determines/settles) reserved for deterministic relations; informs/bears
on/shifts otherwise; his diagnosis (wanting to sound strong; accuracy categorically
first) recorded as lesson 23. NEW section-6-claims.md: all 14 claims with modal force,
ground, status; mechanism tests run THIS SESSION (verify_s6_mechanisms.py): C5b
verified; C6 verified (top-10%-by-K share of optimal gain 0.92/0.86/0.72 at alpha
1.3/2.0/3.5); C7 verified (top-decile recovery at tau=3: 0.68 heavy vs 0.21 light;
middle-pair ordering 0.50-0.54, near chance, the corollary's direct mechanism); C3/C4/
C5a/C8/C12/C13/C14 magnitudes owed to the R re-derivation before lock.

1. K terminology: talent (recommended) or knowledge. Word-level, Aydin's.
2. Approve the bootstrap verification items (schedule sweep, CE check, 200-seed re-runs) if not
   already queued in the running suite.
2. Venue, downstream of the locked story: Research Policy vs philosophy of science.
3. After sweeps land: seed-floor section claim (D3/D4), Story 5 promotion (Tier B), back-loading
   mechanism text (Package C).

## Do not claim

- Any forward advantage in tail_map, regime_map, info_value, horizon_noise, funder_scale, or low-eps
  horizon_scale. Discretization artifacts.
- The T=3 information peak.
- "Quadratic in T" without a fitted exponent; the range tested supports "superlinear".
- "Uniform seed floors reduce output" as a prohibition. The measured cost is 0.03% to 0.40%. State it
  as a price.
- (Added 2026-08-14, Aydin's calibration) "The value of peer review is the cost of allocating
  without knowing the gap." REJECTED as misleading: the funder has partial knowledge of the gap
  before review and imperfect knowledge after, so review's value is NOT the full cost of
  ignorance. The licensed formulation, everywhere (updated 2026-08-14): review's value is the expected
  value of its reduction of the funder's uncertainty regarding the gap. The oracle gap (D-4's 17.2% to 1.7%) may be
  described as the span from records-only to complete information, never as "what review is
  worth."

## 2026-09-08 (verification pass + Figure 4 + section 6 v4)

Aydin's directives executed: tau clause removed from S6's first sentence; "only
moderate" -> "even a signal of moderate informativeness captures most of its value"
(the only/even ambiguity swept). Then the full claims verification: container R
reproduced canonical D_gap_convergence cells BIT-IDENTICALLY (1e-14), licensing
in-session re-derivation. New runs (scripts + outputs in this folder):
run_D4_fine.R -> gap_convergence_fine.csv (17 tau x 3 alpha x 50 seeds, adds
alpha_K=3.5); verify_s6_rounds.R (+_OUTPUT.txt): T=8/T=20 trajectories, paired
sharp-end test.

Verification outcomes now in section-6-draft-v1.md (v4) and section-6-claims.md:
- Corrected: corr at sharpest 0.94/0.96 (not 0.95/0.97); 17.2->1.7% is HEAVY-tailed,
  defaults are 10.4->1.6%; "nearly all" -> "most" (86/93% recovery).
- FALSIFIED and rewritten: "further rounds of records leave the funder no closer"
  (corr rises slowly, 0.08->0.41 over 20 rounds vs review-informed 0.81 in round 1;
  old 0.13-0.18 range was across alpha at round 1, not across rounds).
- New measured wrinkle (footnoted, his verdict pending): value dips ~2% at the very
  sharpest signal under heavy tails (paired z=4.5); flat at defaults.
- Elbow rescoped: tau~2.5 was defaults' half-value point. New footnote: value at 3x
  default noise 82/44/10% (alpha 1.3/2/3.5); half-value noise ~ tau 10/2.5/0.75.
- C10 closed with alpha=3.5 magnitudes (shortfall 3.5% of no-funding output; review's
  max value 2.5% vs 23% heavy).
- Efficacy paragraph rescoped: accuracy-to-noise mapping is field-dependent (AUC 0.54
  = tau~20 at defaults only); closing sentence now the licensed pair (same
  informativeness -> several times the value; same accuracy -> different
  informativeness).
- Overtrust verified (C13/C14) with the 2-SE caveat kept; asymmetry grounded
  (tenfold overtrust loses more than full calibrated value; tenfold undertrust
  retains over half).

FIGURE 4 BUILT: tex/fig4-review-value.tex + fig4-value.dat (share of records-only
shortfall recovered vs noise, log x, three alphas, family style; no AUC marker, see
draft note 8). FIGURES.md updated with simulation-figure conventions. Preview
compiled clean. S6 typeset into main.tex awaits Aydin's verdict on v4.

Open on Aydin: sharp-end reversal footnote (keep/drop/investigate); half-value line
promotion to body; C5c mechanism-isolation run; S6 verdict then typeset.

## 2026-09-08 (second pass: reversal mechanism, S6 typeset, S7 prepared)

Aydin's verdicts executed: reversal INVESTIGATED (kept footnoted, now with mechanism);
half-value line KEPT in footnote; S6 TYPESET.
- Reversal mechanism isolated (reversal_probe 1-3, appended to verify_s6_rounds.R):
  round decomposition refutes exploration (deficit in round 1 itself, z=-3.7);
  intervention on resource noise flips the sign (tau_r=0.01: +1.43, z=+5.8). Cause:
  interaction with the noisy resource estimate. Footnote in S6 updated with the
  mechanism; ledger row C3b added.
- S6 typeset into main.tex (section 6, Figure 4 placed via \input, footnotes 9-15;
  bridge recast to "Next we" per his stated preference, flagged in draft note).
  17pp clean, pages 13-16 inspected. FIGURES.md: fig4 -> placed.
- S5's TODO-flagged footnote CORRECTED in main.tex: old text had the falsified
  0.13-0.18-across-rounds claim and misattributed 17.2% (it belongs to the S5-at-tau-20
  cell, not records-only at defaults). New: corr 0.08 -> 0.41 over twenty rounds;
  records-only shortfall 11.4% of no-funding output (24.6% heavy). Strategy-mapping
  TODO resolved: Table 1's "Myopic, without review" = S4 = the D-4 corr_S4 series.
- S7 PREPARED: all headline numbers re-derived in-container from staged canonical
  sweep files (verify_s7.R + OUTPUT, all match the integration handoff);
  section-7-draft-v1.md (function analysis + plain-first draft + 9 notes incl. title
  candidates, naming decision on "seed-and-harvest", figure candidates ranked) and
  section-7-claims.md (T1-T12, all verified or analytic). Timing-vs-signal ratio
  stated eps-dependently (1/30 at eps=0.3, 1/4 at 0.85).
Open on Aydin: S7 v1 verdict; S7 figure choice (my lean: round-share schedule at
increasing depth); naming and title calls (draft notes 2, 3, 6, 8); repo push.

## 2026-09-08 (third pass: figure-or-proof rule, relative terms, appendix, skill rev 8)

Three new rules from Aydin, all promoted to the skill (EIGHTH revision, delivered; he
must save it): (1) figure-or-proof backing, footnotes alone don't cut it
(anti-pattern 26); (2) results as percent gains over the relevant comparison, never
raw output counts (27); (3) no single-number summaries of sweep-ranging quantities;
state ranges over tested settings (28). Style lessons 24-26; results-section block in
section-patterns; three checklist items.

Actions (audit artifact: figure-or-proof-audit-2026-09-08.md):
- APPENDIX A built (tex/appendix-gap-rule.tex): gap-rule derivation converted from
  gap-rule-proof-v1.md; Lemmas 1-3, Prop 1 restated (starred env), Corollaries 1-2,
  Examples 1-2, Remarks 1-3; all four \S4 Appendix [X] pointers resolved by \ref;
  fig2 caption now refs Corollary~\ref{cor:targeting-vanishes}. Preamble: lemma/
  corollary/example/remark counters made independent (none were used elsewhere).
  Remaining [X]s (2) = owed simulation-specification appendix, TODO-tagged.
- FOUR NEW FIGURES, family style, all data re-derived or re-run in-session:
  fig5-records-rounds (renders Fig 4; T=20 trajectory, S4 vs S5, replaces the C5a
  footnote), fig6-overtrust (Fig 6; D-2 curves as % of no-funding output, 1-SE bars,
  replaces the C13 footnote), fig7-depth-schedule (for \S7; b=3 schedule re-run,
  24 seeds), fig8-timing-boundary (for \S7; resource_regime map),
  fig9-schedule-vs-signal (for \S7; PG/signal ratio, parity line). FIGURES.md notes
  slug-vs-rendered-number mapping (renames skipped: bridge cannot delete).
- \S7 draft v2: his two sentences verbatim; all results in percent terms; the
  one-thirtieth single number retired for the sweep-ranging statement (ratio < 1 in
  all 32 cells; order of magnitude+ at default rate; max 0.66 at T=10, eps=0.85 -- my
  earlier "at most a quarter" was wrong, caught in re-derivation); b gloss fixed
  (depth cells are budget_ref=K). T3/T6/T8/T12 ledger rows restated.
- \S6/\S5 tex: overtrust footnote in relative terms; records footnote removed
  (figure); \S5 display promise recast. Paper builds 23pp clean, appendix inspected.
- Deferred to Aydin: T7 attribution figure (would be \S7's fourth); C5b cumulative-
  Fisher appendix proposition; renames of figure slugs (recommend against).

## 2026-09-08 (fourth pass: S7 typeset, S8 prepared)

S7 v2 APPROVED ("no more figures needed" = T7 figure declined) and TYPESET into
main.tex (section 7, Figures 7-9 placed; fig9 label cluster resolved by caption
identification; boundary footnote's center-of-mass wording corrected at typeset,
flagged in the draft's supersession note). 26pp clean, pages 16-20 inspected.

S8 PREPARED as a function-analysis-first deliverable (section-8-prep.md), NOT a
draft: the lotteries half has design decisions that are Aydin's (and plausibly
Kevin's/Simon's). Verified inputs: Package A re-derived (focal floor cost 6.1% of S1;
kappa 0.10-0.86; convexity; P-A2 exact; D3's heavy-cells-are-sharp-review pairing
caveat found and recorded); D3 relative targeting value extracted ((S5-S2)/(S2-S1):
1.17->0.72 heavy+sharp, 0.74->0.25 base, decreasing in b, matching S2's promise;
absolute %S1 rises over b<=1, objects must be picked explicitly). NEW single-round
scheme-pricing computation (price_lotteries.py, exact expected output, 2000 pops/
cell): full lottery worst everywhere (Jensen); screened schemes capture most of
targeting's value at heavy tails, degrade gracefully in screen noise; within-pool
randomization costs 0.05-0.19 of targeting's value; at even spread + b=0.5 UNIFORM
captures 0.83 and beats every concentrated scheme (the screen is the mistake there).
Jensen proposition drafted for Appendix A (his call). Five decisions listed in the
prep file (computation-as-ground; screened-lottery design; Jensen placement;
abstract-sentence flag; funded-share scope).

## 2026-09-08 (fifth pass: S8 drafted)

Aydin's four S8 decisions applied (single-round computation as ground; screened-
lottery design with wide-split displayed; Jensen in appendix, plain explanation in
body; funded share stated not optimized). Extended computation adds the wide split
and MC SEs (0.002): decomposition at heavy tau=1: the screened lottery's 0.18
shortfall vs the ranked split = 0.14 discarded ranking + 0.04 gamble. DRAFTED
section-8-draft-v1.md (5 paragraphs per approved structure) + section-8-claims.md
(E1-E12, all verified or proven) + fig10-floor-cost + fig11-lottery-prices (two
panels, marks identified in caption) + tex/appendix-lottery.tex (Proposition 2 +
proof + scope remark; renders as Appendix B at typeset). Key modal guard: the
ABSOLUTE "lotteries cheap at even spread" is false in the computation (concentration
is the mistake there; uniform 0.83); the drafted defensibility claim is the
comparative one (chance vs ranking), and the abstract's "cheap exactly where review
is worth little" is flagged for the lock-time abstract pass. Open on Aydin: S8 v1
verdict; ranked/wide "division" vs "split" terminology; then typeset + appendix B
input.

## 2026-09-08 (sixth pass: terminology aligned, S8 typeset, S9 drafted)

"Equal division" adopted (Aydin's call): body/figure/ledger aligned, instances
distinguished by pool. S8 TYPESET into main.tex (Figures 10-11 placed, Appendix B
input after A; Proposition 2 renders on the prop counter). 30pp clean, pages 20-23
inspected. FIGURES.md fully updated (figs 7-11 placed).

S9 DRAFTED (function analysis + draft v1 + ledger G1-G10 + figures 12-13), all
Package B numbers re-derived in-session (verify_s9.R): signal-value structure
survives CD/gamma=-3/Leontief and scales with complementarity (6.8/29.3/33.7 %S1 at
heavy+sharp); CD kills scheduling (reused verify_s7 numbers); concentration is
budget-conditional with the interior max at even spread (peak ~gamma=-6, z=3.5 at
400-seed refinement). Figure 13 redesigned mid-build: absolute Gini hid the shapes,
now displays change-from-Cobb-Douglas with absolute levels in the caption (flagged,
draft note 2). Object discipline: "the informed funder's grants," never "the optimal
allocation," for the tierB object. Held out: family seed-floor costs (strategy-pair
column check owed) and correlation robustness (not re-derived); ledger G7-G8.
Open on Aydin: S9 v1 verdict (notes 1, 2, 4); then typeset; then S10 (discussion)
function analysis.

## 2026-09-08 (seventh pass: S8 rewritten for plainness; skill rev 9)

Aydin's critique of the typeset S8 (unmotivated opening; "payline" and the whole
paying/pricing frame; undescriptive "this/that"; inverted lists; neologisms like
"floored share"; run-on sentences). S8 REWRITTEN in main.tex: motivated opening
(targeting vs the two alternatives, each with its rationale, then what the model can
say); every sentence names its subject; plain terms ("fraction of the budget given
out as seed grants," "selection by review," "partial lottery," "lottery over all
researchers"); losses stated in output units relative to named comparisons; short
sentences. The convexity claim made concrete and verified across all eight D3 curves
(half the budget as seed grants loses more than twice what a quarter loses). Fig 10
and Fig 11 labels/captions and the Appendix B remark aligned; the S7 bridge recast.
31pp clean; no "screen"/"payline" remains in body or captions.

Skill NINTH revision (anti-patterns 29-32: undescriptive anaphora; coined labels;
metaphoric framing vocabulary incl. inverted lists; run-on sentences; four checklist
items); style lessons 27-29. Delivered; he must save.

FLAGGED for Aydin (approved earlier text, same objection): the pricing frame in S1's
roadmap ("S8 prices seed grants and lotteries"), S2 (four instances: "prices what a
lottery saves," "at low paylines," "the payline," "What remains unpriced," "price the
choice field by field"), and S4 ("it sets the price of overriding the optimal
allocation"). Plain replacements proposed in chat; not edited without his word.

## 2026-09-08/09 (eighth pass: price-language sweep; S9 dissolved into S8 + Appendix C)

Paper-wide sweep replacing price/pay/buy/cheap/ledger/payline/floor language with
output-loss language (abstract, S1-S8, fig6/fig10/fig11 captions and comments; label
fig:lottery-prices -> fig:lottery-schemes; literal "review costs"/"costly effort"
kept). Commit 2253df2.

Aydin's decision on S9: no section. Robustness content -> NEW tex/appendix-technology.tex
(Appendix C: CES family defined as the power mean A*M_gamma(K,R); Fig 13 signal
robustness; Table 2 scheduling value + schedule center of mass by technology at T=5,
replacing the S7 footnote that cross-quoted greedy-allocator magnitudes beside the
body's). Concentration result -> closing movement of S8, retitled "Spreading funds:
seed grants, lotteries, and concentration" (Fig 12; non-monotone clause kept at his
word; bottleneck-top-up mechanism hedged as "a pattern consistent with"). Pointers
updated: S1 roadmap (S9 discussion, S10 conclusion), S2 x2, S3, S4, S6 (new pointer
sentence), S7. fig13 label collision fixed (zero rule removed). Preview 34pp clean.
Float note: Fig 12 currently lands at the top of Appendix A's first page because no
S9/S10 body text follows S8 yet; will settle when the discussion lands.

Device bridge disconnected mid-turn; write-back queued in memory
(/areas/for-claude-writeback-2.md, AH cont. 12) and changed files delivered in chat.
Next: S9 (discussion: implications + limitations) function analysis.
S9 Discussion v2 TYPESET into main.tex (\section{Discussion}, label sec:discussion;
roadmap pointer converted to \S\ref{sec:discussion}); 37pp preview clean. Not made:
the one-clause S6 consistency edit (awaiting Aydin's word). Next: S10 conclusion.
Overleaf bundle delivered (main.tex + references.bib at root, figures/ for the rest;
sources unwrapped to one line per paragraph, text verified identical). Full style and
clarity audit written: style-audit-2026-09-27.md (6 READ, 3 VERDICT, 14 global terms,
~95 local items with proposed rewrites); no changes made, awaiting Aydin's approvals.
Audit APPLIED (2026-09-27) with Aydin's exceptions (kept: "how to read the evidence on
review's efficacy" once; "underdetermine" throughout; "Consider the challenge of the
grant funder"; "sits at the intersection of four literatures"; "a dynamic that formal
models show"; "effectiveness requires concentrating"; "epistemic landscape"; his own
wording for the two-assumptions sentence; "What applies to a real field is the pair of
features that organizes the results." without the second clause; "each perform better
in some fields than in others"). All other items applied to main.tex, Appendix A,
Appendix C, and captions (tranche->budget, purse->budget, seeds->simulated populations,
depth->grant size, production form->production function, harder->stronger
complementarity, efficacy->how well review predicts research outcomes, uncertainty
definition->output definition, etc.). Preview 37pp clean; Overleaf bundle regenerated.
Sources are now stored unwrapped (one line per paragraph). Skill TENTH revision:
anti-patterns 33-35 (straightforwardness audit with his named exceptions; simulation
reporting language; one operational definition) + three checklist items + SKILL.md
rule 13. Delivered; he must save.
Compression plan written (compression-plan-2026-09-27.md): layout accounts for 7 of 37 pages; content plan cuts ~30% (captions to a third, footnotes by 60%, S9 and S3 and S8 trimmed, Appendix A compressed, new Appendix D absorbs generating context). Awaiting Aydin's decisions.

## 2026-09-27 (tenth pass: compression, layout kept)

Aydin's instructions: keep the layout (1.5 spacing, 4 cm margins) and figure sizes;
results may stay stated in abstract, intro, and conclusion; prune in the conclusion;
merge S5 into S6; keep a separate S10 conclusion; two-sources-of-growth to one line;
sharp-end reversal cut entirely; Heesen sentence cut. Executed: S2/S3/S4/S5/S6/S7/S8
trimmed per the plan; S5 merged into S6 as "What the funder can learn: track records
and peer review" (sec:review, with sec:records as a starred subsection label);
hard-coded section numbers converted to \ref; S9 restatement paragraphs cut to one;
NEW S10 Conclusion (200 words, big picture); all 13 captions rewritten (1311 -> 667
words) with settings moved to NEW Appendix D (appendix-simulation.tex: budget
normalization 2bnE[R0], A = 1/2 in the paper's notation, defaults, scoring, per-figure
settings table, code pointer); Appendix A compressed (2179 -> 1737: Lemma 2 proof
shortened, Corollary 1 stated without proof, remarks deleted, example commentary cut);
Appendix C opening trimmed. Body 9514 -> 8356 words, footnotes 1377 -> 767, captions
1311 -> 667, appendices 3227 -> 3405 (incl. new D). Pages 37 -> 34 at the fixed
layout (body 27 -> 24). Under 20 pages is not reachable at this layout without
removing content that carries results; options reported to Aydin.

## 2026-09-29 (intro refined; PNAS plan)

Intro refined per Aydin's draft + my memo (intro-refinement-2026-09-29.md): stakes hook
with the number in the first line; new second paragraph (why a model; four questions;
theoretical contribution); funder-facing sentence in the results paragraph; his wording
"how seed grants and lotteries perform relative to targeted funding" (intro + roadmap)
and "concentrated among a few researchers" (global). Commit 4aff978. 35pp.
PNAS plan written (pnas-plan-2026-09-29.md): inventory of 30 results with main/SI
disposition; four composite figures; 4,000-word budget (intro 480, Results 2,300 in
four subsections, Discussion 750, Methods 470); SI structure; decisions (a)-(d).
Appendix A proofs verified by hand and examples recomputed exactly (commit 611a5d6).
PNAS DRAFT v1 written (tex/pnas/pnas-main.tex + fig-p1..p4 composites from existing
panels, stacked vertically with A/B/C labels; article class for review, port to
pnas-new.cls at submission). Word counts: intro 583, Results 1,777, Discussion ~880,
Methods 474 = ~3,700 main text; Significance 120; abstract 145. Fig 1 = frontier only
(Aydin's call); Fig 2 (targeting in budget) -> SI as Fig S1; overtrust S2; grant-size
schedule S3; concentration S4. "Two heuristics a funder might use" per Aydin. SI Appendix
not yet assembled. Bridge dropped before write-back; delivered in chat.
Same day, later: title of both versions changed to "A model of optimal science funding:
targeting the capability-resource gap" (Aydin). Coined-term sweep across main.tex,
appendices, captions, and the PNAS draft (term-replacements-2026-09-29.md: sharp ->
reliable, coarse -> noisy, thin -> small, resource-poor -> few resources, lean ->
shift toward later rounds, exceptional few -> few researchers of unusually high
capability, fundable -> above a review threshold, money is tight -> budget small
relative to researchers' resources, stand-in -> proxy; "Less than information" ->
"less than the signal provided by grant peer review"; "pair of features" ->
"qualitative features of the dynamics"; predictions paragraph now names the tests).
Kept for his call: overtrust/undertrust, tight/ample budget, no-funding output,
records-only. PNAS draft v2 moved onto the amsart template (his layout), figures
inline at their sections, 14pp; Figs 7/10/11 and PNAS Fig 4 narrowed so their
outside labels stay inside the margin. Both files compile with no overfull boxes.
Same day, second pass (Aydin's six requests): "using a Bayesian approach"; "research
output" at first mention per paragraph in framing sections and captions; NORMALIZATION
CHANGED from percent of no-funding output to percent of the research output that
funding adds (funder named at each use): seed grants 11/17/3.5% (was 6/1), schedule
0-4.5% vs signal 7-27% over 30 settings (was 0-0.7 vs 0.8-7.5, 32), records-only
shortfall 27/40% of complete-information gain, review max 9 vs 39%, App C 16/47/51%,
Table 2 to 6.52; figs 6, 10, 12 re-plotted (data regenerated from staged sweeps;
review cells rerun, review_cells.csv). Heuristics sentence in his words; track record
is a signal (review "a further signal"); "In practically all cases"; Fig 1 redesigned
with two frontiers (slopes 3/2, 2/3; black/gray/open) in both versions; PNAS panel
letters raised; b renormalized so b = 1 is the field's one-round baseline resources
(budget = b n E[R0]; all b doubled; code's b' = b/2 noted in App D); SI Appendix built
into pnas-main.tex (si-*.tex, figS1-S5, claim-to-section table; refs \ref-based;
Proposition 2 counter fix). Appendix C retitled "production function". Long version
bibliography regenerated (six cited keys had been missing from the stale bbl; chicago
bst not installed here, bbl built with plainnat, source keeps chicago). Both compile
clean: main 35pp, pnas 34pp (14 main + SI). Details: term-replacements-2026-09-29.md.
Claims ledgers E1-E12, T-rows, G-rows still carry the old normalization; superseded.
Same day, third pass: PNAS subsection "When to fund" -> "Strategic timing of funding"
(and Fig 3 caption lead); long version keeps "The timing of funding" (his instruction
named the subsection head). Coined-term rule promoted: saved to his Claude preferences
(applies on every surface) and proposed as standing rule 3 + hot rule 14 (comparative
metrics) of the paper-writing skill (rev 11 card; he must save). Overleaf bundles
rebuilt: overleaf-main-2026-09-29.zip (main.tex + references.bib at root, 17 tex + 18
dat in figures/, chicago bst) and overleaf-pnas-2026-09-29.zip (pnas-main.tex + bib at
root, 16 tex + 18 dat in figures/); both test-compiled clean with paths rewritten.
Same day, fourth pass: the factor 2 dropped from the production function. lambda_i =
A K_i R_i/(K_i + R_i), "proportional to the harmonic mean"; every 2A -> A (about 50
places: marginal value A K^2/(K+R+g)^2, bound A K_i, lower bound K > lambda/A,
small-grant approximation A R, c = sqrt(A/nu) - 1, nu(B) = A/(1+c)^2); simulations set
A = 1 (the code's formula, no conversion note); Appendix A examples recomputed at A = 1
(285/53, 284/51, 40/7; 56/15, 49/12, 21/5; c and g* unchanged; verified exactly). PNAS
Results defines A at first use. Fig 3 axes are schematic, data unchanged. Both compile
clean (main 36pp, pnas 34pp). Overleaf bundles rebuilt.
Same day, fifth pass: NEW appendix figure (fig14-heuristics.tex; PNAS figS-heuristics,
SI Text 3): the two intuitive schemes on the Fig 1 population, the smaller-frontier
budget B1 = 8.52 split equally among nine researchers; left = nine highest expected
output, right = nine fewest baseline resources; gains over no funding on this
population (A = 1): gap rule 4.01, track record 1.98, under-resourced 2.09, uniform
1.99 (computed in-session). Parenthetical references added in S4 (long) and Results
(PNAS). Data: fig1h-all/track/under.dat. Both compile clean; Overleaf bundles rebuilt.
