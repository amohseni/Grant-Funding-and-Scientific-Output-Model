# Style and voice lessons from Aydin's intro edits (2026-08-06)

Derived from his line edits to intro v2. Each lesson: the edit, the principle behind it, and
how to apply it. These extend the paper-writing skill's rules; candidates for promotion via
skill-gaps once they hold across another section.

## 1. Calibrate cited evidence to the argument's own commitments

"No better than chance" became "barely better than chance." Two reasons. Accuracy: AUC 0.54 is
barely above chance, not chance. Coherence: a paper whose model assumes review carries signal
cannot lean on "review is random." The scope of the finding matters too: Fang et al. measured
funded grants only, which is compatible with review filtering out weak applications. Cite the
counter-evidence (Li and Agha 2015; Park, Lee, and Kim 2015) so the dispute is presented at
its true strength. Rule: before citing evidence against a position, check what the paper's own
machinery assumes; present disputes as disputes.

## 2. Motivation is experiential; the model is abstract

"Reduced to its core, the decision is this: a funder with a fixed budget faces a population of
researchers" reads as model description in the motivation slot. His fix puts the reader in the
funder's shoes, second person: "Consider the challenge of the grant funder: you have a fixed
budget and must choose how to allocate it..." Rule: the introduction's opening moves make the
reader feel the decision problem; abstraction begins only at "We model this decision situation
as follows."

## 3. Jargon only where it is a term of art doing work

"Two latent inputs" became "two ways." Latent, inputs, operational, full information: model
vocabulary, not needed in motivating prose. Terms of art that stay: capability-resource gap,
Bayesian funder, heavy-tailed distribution. The test: does the term carry technical content the
sentence needs, or does it just sound technical? Related cut: "funding frontier" moved from the
intro to the results section, along with "water-filling."

## 4. Spell out compressed technical references where clarity is at risk

"Under full information" became "In the ideal condition where the funder has complete
information regarding researchers' capabilities and resources." Not everywhere; wherever a
reader might not decode the shorthand. The expansion costs words and buys accessibility; at
first occurrence, pay the cost.

## 5. Method labels wait for the results sections

Water-filling is a fact about the derivation, not part of the motivation. The introduction
does not advertise technique. Name the method where the derivation appears.

## 6. Never claim novelty explicitly

"The rule, as a statement about science funding, is to our knowledge new" is out, always. Two
reasons: unprofessional (an author is expected to know their field, so "to our knowledge"
announces a deficiency), and a digression (novelty positioning is the literature section's
job, done by showing, not saying). The introduction motivates and sets the stage; it does not
stake priority.

## 7. No result numbers in the introduction

Extends the abstract rule. Empirical stage-setting numbers (budget sizes, concentration
shares) stay; the paper's own quantitative results appear qualitatively, with "we show" and
"we find" pointing forward.

## 8. State results as claims, method-neutral and comparative

"In simulation, a funder relying on publication records alone allocates far from the gap rule"
fails twice: it ties the claim to machinery not yet introduced, and it localizes to one method
what is also true elsewhere. His form: "We show that a funder relying on publication records
alone performs far worse than one that accounts for the gap." Rule: in the introduction, a
result is a claim the paper will establish, not a report of how it was established.

## 9. Do not foreclose the paper's own arc

"The rule is not operational, because the gap cannot be observed" issues a verdict the later
sections overturn (estimation makes it operational). His fix describes the epistemic situation
instead: "In reality, the funder lacks perfect access... a researcher's track record
underdetermines these quantities." Rule: characterize the obstacle in terms that leave room
for the solution the paper provides. "Underdetermines" is the precise word: the record is
evidence, just not enough.

## 10. Modality precision on what things are vs what they can do

"That is what peer review is, in this model: a noisy measurement" became "That is what peer
review can provide: a more or less noisy signal." Capacity claims over identity claims; and
"more or less noisy" is graded, matching the model's noise parameter, where "a noisy signal"
fixes a property. Say what the object can do at the strength the model licenses.

## 11. Decision-theoretic precision is welcome where it is the content

"The value of peer review for a Bayesian funder is the expected value of its contribution to
the funder's estimate of the capability-resource gap." Precise, technical, and his own
phrasing: the term-of-art exception in action. When the precise formulation IS the claim,
write it precisely; plain language elsewhere.

## 12. Spell out populations at first use

"Heavy-tailed fields" became "fields with heavy-tailed distributions of capability." The
compressed form can return after the spelled-out form has appeared.

## 13. Clarity beats compression where they conflict

Several of his edits lengthen sentences. "Omit needless words" means omit words doing no
work, not maximize density. A longer sentence whose every word aids decoding beats a shorter
one that must be re-read.

## 14. Fixed conventions (Aydin's standing forms)

- The roadmap paragraph always opens "We proceed as follows." and uses \S for "section"
  (\S2, \S3, ...) in LaTeX.
- Result-contrast sentences favor explicit parallel construction: "strategic timing of when
  one provides funds matters far less than strategic targeting of whom one provides funds."
  The parallelism (timing/targeting, when/whom) carries the contrast; a compressed form
  ("whom beats when") stays in the abstract, not the body. When such a sentence opens a
  paragraph, cut any later sentence restating the same contrast.

## Consistency items his edits raised (flagged, his call)

- He wrote "resource-capability gap" once; title and abstract say capability-resource gap.
  Treated as a slip; locked order kept.
- He wrote "options" for the two intuitive allocations (was "rules"). Adopted: options.
- He wrote "uniform seed grants" (was "uniform floors"); matches the model's uniform seed
  funding strategy. Adopted in intro v3, but the LOCKED ABSTRACT says "uniform funding
  floors"; one concept, one term requires changing one of them. Proposed: abstract line
  becomes "Uniform seed grants and lotteries are cheap exactly where review is worth little."
- His floors sentence had "always incur a cost... floors are costless where": smoothed to
  "that cost approaches zero where" to avoid the literal contradiction.
- He wrote "capacities" once; "capabilities" kept.

## 15. Introductions preview findings at sentence strength (calibrated 2026-08-14 against his three published intros)

HARKing runs ~600 words plus roadmap: hook, thesis with its determining variable, method in
one sentence, scope disclaimer, refutation schema, roadmap. Causation spends its length on
dialectic (credentials, precedent, obstacle, vignette), never on walking results. Neither
previews section content at paragraph length. Rule: each result gets one or two sentences at
full strength in the introduction; mechanisms, condition lists, and law structures live in
their body sections. An introduction paragraph that explains WHY a result holds is a
miniature section and gets compressed. Target scale: ~600-700 words plus roadmap.

## 16. Bans and preferences from the 2026-08-14 line-edit pass

- "Reads" banned wherever it stands for "accounts for," "integrates," or "tracks."
  Replacement of record: "accounts for."
- "Assumptions," never "commitments," for a model's substantive premises.
- "Our model," never "the model," in the paper's own prose.
- NO formulas in the introduction at all (supersedes the earlier one-formula allowance):
  the gap rule is described in words ("a grant increasing in their capability and
  decreasing in their existing resources").
- Never use a term of art the paper has not introduced: "a field's technology" cut because
  technology is nowhere defined; replaced with plain description.
- "Peer review supplies an additional signal of capability": publication records are also a
  signal of capability, so "what the record cannot supply" misdescribed the contrast.
- "In spite of this, basic dimensions of allocation are contested." (his P1 line);
  "Whether peer review is effective is disputed."; "Each option can be intuitively
  compelling."; "even a coarse signal captures most of this value."; "The same reasoning
  elucidates the cost of..."; "We show how strategic timing..." (all his phrasings, applied).

## 17. The forward-only reading test, and form never precedes content (2026-08-14, after a \S2 failure)

The failure: "This dispute is typically conducted as a dispute about the returns to funding,
holding the allocation mechanism fixed. Our model treats the two questions as one..." No two
questions had been posed; the aside was not a question; and the unification did not describe
anything the model does (we propose no mechanism; the true claim was budget sensitivity).
Two rules follow.

- Forward-only reading test: every referring phrase ("the two questions," "this tension,"
  "both approaches") must have an antecedent the reader can point to in the text already
  read. Before writing "the two X," confirm the text has posed two X. A deliberate forward
  hook is licensed only if the next paragraph's opening resolves it explicitly.
- Form never precedes content: rhetorical schemas (treats-two-as-one, X-is-really-Y,
  mirrored antitheses) are compressions of content that already exists in plain form. Write
  the claim plainly first; compress only if the schema fits the plain version exactly. A
  pivot sentence that sounds decisive but whose precise claim cannot be reconstructed from
  the preceding text is manufactured; delete it and state the finding. Template pressure
  (every paragraph "needs" a punchy positioning pivot) is the usual source.

Also this pass: "the efficacy of peer review" is the term for the disputed property,
throughout; study descriptions keep factual verbs.

## 18. The order of virtues (2026-08-14, Aydin's constitution; now hot rule 0 in the skill)

His words: "Prose is about clarity above all: clarity in thought and clarity in
communication. The beauty of writing is first and foremost in its realization of the form of
clear thought. Beyond that, simplicity, straightforwardness, and precision and truth in
every sentence." Lexical priority: truth and warrant; clarity; simplicity and precision;
only then style. Enforced by plain-first drafting (generation) and the paraphrase test
(audit); both now in core-style.md and the .skill package. Diagnosis this codifies: most of
his corrections this session traced to punch-seeking (jargon, density, confusing pivots);
the skill's own punch-rewarding rules (bold closes, verdict resets) are now explicitly
subordinated to clarity.

## 19. Numbers only with their generating context (2026-08-14 \S4 session; promoted to the skill 2026-09-04)

The \S4 failure: "improves output over no funding by 23.5 percent" with no parameters,
population, or provenance specified. A number the reader cannot evaluate is an assertion,
not evidence. Two rules, generalizing lesson 7 beyond the introduction. A quantitative
result appears only where its generating context (parameter values, population, provenance)
is specified or explicitly pointed to. And a section whose authority is analytic carries
only claims derivable from its formal results: the \S4 fix replaced the simulation numbers
with exact counterexamples and a limit corollary in the proof appendix, banked the numbers
for the quantitative sections, and in the process caught a false sentence ("the optimal
allocation converges toward uniform funding") that the derivability requirement exposed.
Also his qualification at promotion: no formulas or numbers in the abstract either (the
abstract and introduction now share one rule).

## 20. Guide, do not proclaim (2026-09-05, promoted at his instruction)

"Everything turns on the marginal value of a dollar" became "Consider the marginal value
of a dollar." The former is bombastic, the latter is clear and guiding. Rule:
importance-announcing openers ("Everything turns on X," "The whole question is X") are
bombast; a guiding construction directs the reader's attention to the same object plainly.
Now anti-pattern 23 in the skill.

## 21. The body carries the flow; footnotes carry the rest (2026-09-05, promoted at his instruction)

The water-filling and poverty-gap precedent sentences, and the pointers to appendix
examples, moved from body to footnotes. His rule: body prose is for the flow of
information and understanding; history, precedent, and qualifications that would break the
flow and are not strictly necessary go to footnotes. A citation whose only role is
precedent does not earn a body sentence. Now a core-style addendum block.

## 22. Formal objects in their own vocabulary (2026-09-05, promoted after a run of failures Aydin caught in one section)

The pattern behind "the budget disciplines the targets," "picture c rising," "a small
budget stops c early," "prices both at once," "output saturates at 2AK," "the ceiling
their capability sets," and "at capacity": figurative language imported into the
description of formal objects. The metaphor asserts content the formalism does not
license, and the damage runs from vagueness (disciplines) through category error (a
constant defined by an equation described as moving or being stopped) to near-falsehood
(ceiling and saturates assert an attained plateau; in truth expected output increases
with every dollar and is bounded above by 2AK, which it approaches but does not attain).
Rule: describe formal objects in the vocabulary of their definitions (increases,
decreases, is bounded above by, converges to, is determined by, is defined as).
Enforcement, the substitution test: for every sentence about a formal quantity, write the
statement of the formalism it translates; no statement, or surplus over it, means
rewrite. Sub-rule: "increases"/"decreases" for monotone change, never "rises"/"falls"/
"grows" (his instruction, applied through abstract and \S2 as well). This is the
generalization his ceiling correction demanded: it is interpretive fidelity applied to
every informal gloss, not only to result statements.
