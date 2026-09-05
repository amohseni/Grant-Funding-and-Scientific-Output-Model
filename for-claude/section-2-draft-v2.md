# Section 2 (Related literature), draft v2 (2026-08-14)

Revision implementing the reflection pass (all items approved by Aydin). Changes: lotteries
paragraph rebuilt around the cost-side/value-side complement (Gross and Bergstrom moved here);
timing-coverage sentence added to the models paragraph; Sorzano characterized precisely;
warrant items softened; water-filling and poverty-gap precedents moved out (now in
section-4-notes.md); closing chain replaced with a synthesis close; "our model defines the
quantity at issue"; "long-running dispute"; "a single quantity" for rival models. [check]
marks citations not yet read in the original.

---

## 2. Related literature

Our contribution sits at the intersection of four literatures: the empirical debate over
concentrating or dispersing funds, the evidence on the efficacy of peer review, the case for
funding lotteries, and formal models of science funding.

Whether funders should concentrate resources on a few researchers or spread them across many
is a long-running dispute, and concentration is the current practice: the top 1% of NIH
principal investigators receive close to 10% of research project grant funding, a share that
has risen since 1985 (Katz and Matter, 2020; Lauer and Roychowdhury, 2021). Empirical studies
of funded portfolios generally find diminishing marginal returns to funding per researcher,
and reviews of this literature conclude that the evidence favors more dispersal than current
practice (Mongeon et al., 2016 [check]; Aagaard et al., 2020 [check]). Concentration is also
self-reinforcing: early grant winners accumulate advantages beyond the grant itself (Bol, de
Vaan, and van de Rijt, 2018), a dynamic that formal models show can reward questionable
research practices over a career (Heesen, 2024). This dispute typically concerns the
returns to funding: whether output per researcher increases or decreases as grants concentrate. In
our model, the answer also depends on the size of the funder's budget relative to
researchers' total resources: with a budget large enough to relieve most bottlenecks, the
optimal allocation itself reaches nearly every researcher, so little is lost by spreading
funds widely; with a tight budget, effectiveness requires concentrating funds on the
researchers with the largest capability-resource gaps (\S9). The value of any targeting, in turn, depends on what the
funder can know about researchers (\S6).

The efficacy of peer review is contested. Among funded NIH grants, percentile scores predict
subsequent productivity barely better than chance (Fang, Bowen, and Casadevall, 2016), and
reviewers assessing the same application agree only weakly (Pier et al., 2018 [check]).
Analyses with different designs find real but modest predictive power: better scores predict
more publications, citations, and patents, including under controls for investigator history
and institution (Li and Agha, 2015), and a funding natural experiment finds higher-scored
proposals outproduce lower-scored ones (Park, Lee, and Kim, 2015). Both sides of this dispute
measure review's accuracy: how well scores predict outcomes among funded applicants. Our
model instead measures review's value: its contribution to the funder's estimate of the
capability-resource gap. That value depends on the field; the same noisy signal can be worth
a great deal in one field and nothing in another (\S6), and no single number settles the
question.

Doubts about review have made lotteries a live proposal. Advocates argue that randomizing
among fundable applications saves review costs, reduces bias, and honestly reflects review's
limited ability to discriminate among fundable proposals (Fang and Casadevall, 2016 [check]; Avin, 2019a; Roumbanis, 2019
[check]); several funders have implemented partial lotteries: the Health Research Council of New
Zealand has allocated its Explorer Grants by lottery since 2013 (Liu et al., 2020), the Swiss
National Science Foundation randomizes among ties (Heyard et al., 2022), and a German funding
line runs a lottery before peer review (Luebber et al., 2025). The existing case prices what a lottery saves:
Gross and Bergstrom (2019) model grant competition as a contest in which scientists invest
costly effort in proposals, show that at low paylines this effort can rival the scientific
value of the research funded, and show that a partial lottery among proposals above a
threshold decouples the waste from the payline. What remains unpriced is what a lottery
forgoes. Our model supplies that half of the ledger: the value of targeting, which runs from
negligible, where budgets are ample, signals uninformative, or capability evenly spread, to
substantial, where capability is heavy-tailed and budgets tight (\S8). Together the two
halves price the choice field by field.

Formal models of science funding have concentrated on the incentive costs of competition and
on the value of exploration. Beyond the contest model above, Gross and Bergstrom (2025) find,
in a game-theoretic model of contests for scarce rewards, that intensifying competition
pushes scientists toward higher-risk, higher-return projects. Avin (2019b) simulates
funding strategies on an epistemic landscape and finds that allocation with an explicitly
random component significantly outperforms a peer-review-like strategy on large landscapes,
because review concentrates effort on regions already known to be productive. Sorzano and
Pueche-Granados (2026 [check]) treat allocation as a decision problem under heavy-tailed
uncertainty about project value, pairing an empirical bibliometric signal with a
biased-lottery mechanism. In each of these models, what funding buys is a single quantity, a
project's chance of success or a researcher's productivity, and none treats the timing of
spending as the funder's choice [check: Avin's simulation dynamics]. In our model, a
researcher's output increases in both capability and resources but is limited by the scarcer
of the two, and the funder observes neither directly.

Our formulation and analysis of the capability-resource gap allows us to address questions
raised by these literatures: how widely funds should be spread (\S9), what peer review is
worth (\S6), and when lotteries are defensible (\S8).

---

## Notes for Aydin

1. VERIFICATION LEDGER (updated 2026-08-14 after the check run):
   VERIFIED: Gross and Bergstrom 2019; Avin 2019b BJPS; Park, Lee, and Kim 2015; Pier 2018
   (PNAS 115(12):2952-2957); Bol et al. 2018 (PNAS 115(19):4887-4890); Mongeon et al. 2016
   (Research Evaluation 25(4):396-404; AUTHORS CORRECTED: Mongeon, Brodeur, Beaudry,
   Lariviere); Aagaard et al. 2020 (QSS 1(1):117-149; conclusion CONFIRMED, quotes: "rather
   strong inclination toward arguments in favor of increased dispersal", "most systems
   currently have moved too far toward concentration"); Sorzano and Pueche-Granados 2026
   (arXiv 2604.22793; characterization CONFIRMED: percentile-normalized bibliometric signal,
   biased-lottery framework, heavy tails, single-quantity researchers); NZ HRC Explorer
   Grants lottery since 2013 (Liu et al. 2020, Research Integrity and Peer Review 5:3);
   Heyard et al. 2022 (Statistics and Public Policy 9(1); pages still open); Avin 2019a
   Mavericks (SHPS 76:13-23); Roumbanis 2019 (ST&HV 44(6):994-1019).
   STILL OPEN: Fang and Casadevall 2016 mBio (fields consistent across secondary sources;
   DOI check at submission); the German-funding-line lottery study (Nature Communications
   2025, s41467-025-65660-9; confirm authors and that the scheme is Volkswagen's
   Experiment!); Avin's simulation dynamics for the none-treats-timing sentence (abstract
   and table of contents show per-cycle allocation with selection strategies, best/lottery/
   triage, which supports the claim, but the full text needs reading); Katz and Matter
   published version (working paper entry stands).
2. Precedents (water-filling, Cover and Thomas; poverty-gap transfer) moved to
   section-4-notes.md for placement beside the derivation.
3. The synthesis close replaces the contribution chain; the intro no longer gets repeated.
4. Warrant softenings: "orders of magnitude" removed (now negligible-to-substantial as the
   model's own comparative static); "weakest exactly where it is most often pressed" removed.
5. Third pass (clarity lens, paraphrase test on every sentence): close rewritten per Aydin's
   gloss (the gap lets us address the literatures' questions; "meet in one object" metaphor
   out); efficacy paragraph now contrasts accuracy (what the literature measures) with value
   (what our model measures), REVISING the previously approved "defines the quantity at
   issue" line, easy to revert; "the value of selection" collapsed into "the value of
   targeting" (one concept, one term) and the redundant sentence merged; "concedes honestly
   what review cannot discriminate" and "epistemic value of diversity" replaced with plain
   statements. All other sentences passed the paraphrase test.
6. 2026-08-14 second pass: "efficacy of peer review" is the term for the disputed property
   (funnel + topic sentence); study descriptions keep their factual verbs ("predict
   subsequent productivity"). Concentration positioning rewritten after the
   two-questions-without-antecedent failure; forward-only reading sweep done on \S2 and the
   intro (all remaining referring phrases have textual antecedents; P2's "two quantities"
   forward hook resolves in the next paragraph's first sentences and is deliberate).


## Final verification run (2026-08-14, from Aydin's folder)

- Avin 2019b BJPS, FULL TEXT: timing sentence VERIFIED ("Funding is represented as a process
  of selection. In every time step... the modelled funding mechanism selects from this pool";
  mechanisms differ only in whom they select). Characterization verified (random allocation
  beats the peer-review-like Estimated Potential mechanism on large landscapes; review biases
  toward conservative followers). Bonus: his agents are identical ("ignores natural ability";
  funds projects, not people), which strengthens our single-quantity contrast.
- Luebber et al. 2025 VERIFIED: Nature Communications 16:9824; lottery-first then peer review
  at "a large German funding organization" (unnamed in the paper; our sentence matches);
  68 percent lower estimated cost, more funded projects from female applicants.
- Heyard et al. 2022 VERIFIED: Statistics and Public Policy 9(1):110-121.
- Sorzano and Pueche-Granados VERIFIED from the PDF: two authors; deterministic and
  stochastic impact-based mechanisms converge to high concentration; biased lottery balances
  exploration and exploitation. Our characterization stands.
- Poverty-gap source chosen by the folder: Foster-Greer-Thorbecke. FGT 1984 (Econometrica
  52(3):761-766) added to the bib for \S4; the 2010 retrospective is on hand.
- STILL OPEN (not in the folder): Fang and Casadevall 2016 mBio original; Katz and Matter
  published version; Li and Agha original (the "under controls" clause); Park, Lee, and Kim
  original.
