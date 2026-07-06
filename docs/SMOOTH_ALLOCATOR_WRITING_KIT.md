# Writing kit — what the exact (smooth) allocator changed, and how to write it up

*Written 2026-07-06 after re-running the full manifest with the smooth allocator. Plain-language,
built to let you move straight into writing and presenting. Read the TL;DR first; the rest is detail
you can pull from as you draft.*

---

## TL;DR (read this first)

We replaced the old grant allocator — which filled each budget in coarse discrete steps — with an
**exact** one (continuous optimization, no step size). We did this because a diligence pass found the
old allocator was inventing a small spurious "forward advantage" at certain budget values (a
discretization artifact). The exact allocator removes that artifact **by construction**.

When we re-ran everything with the exact allocator, two things happened:

1. **Your central result got stronger and cleaner.** Forward (look-ahead) funding beats myopic funding
   **when knowledge compounds** (ε > 0), and the advantage **grows with the funding horizon** — rising
   to about **+2.0 extra publications** at strong compounding and 5 rounds. The headline figure is now a
   clean, monotonic surface instead of a slightly bumpy one.

2. **Several secondary "advantages" disappeared.** In regimes with little or no compounding — heavy
   tails, budget-scale, funder-scale, information-value, noise — the old allocator showed forward
   advantages of +0.15 to +0.50 that were **mostly discretization artifacts**. With the exact allocator
   they collapse to ~0. This is confirmed independently: even a *fine-grained* version of the old
   allocator agrees they should be ~0.

**What this means for the paper:** your core thesis holds and is actually sharper — *forward planning's
value is a compounding phenomenon that scales with horizon.* But the secondary regime claims need to be
rewritten or dropped. This is a genuine revision, not a number refresh. The good news: it's a cleaner,
more defensible story, and we caught it before you presented.

---

## The one-sentence result (for your abstract / title slide)

> **Forward-looking grant allocation outperforms myopic allocation precisely to the extent that
> knowledge compounds, with the advantage accumulating over the funding horizon; apparent advantages in
> non-compounding regimes were artifacts of discrete allocation and vanish under exact optimization.**

---

## What changed under the hood (methods, in plain terms)

- **Old allocator ("greedy"):** poured each round's budget in `n_steps` equal chunks, each handed to the
  researcher with the current highest marginal return. With `n_steps = 50` (the default), the chunk size
  was coarse. At special budget values the chunk grid landed differently for the forward vs. myopic
  schedules, opening a gap between two allocations that should be nearly identical. That gap is spurious.
- **New allocator ("smooth"):** solves the *exact* budget-constrained optimum — water-filling for a
  single round, projected-gradient ascent for the multi-round forward plan. No step size, no grid, so the
  artifact cannot occur.
- **We verified the new allocator is correct**, not just different: it reaches a higher value of the same
  objective than even a very fine old allocator (`n_steps = 800`), and satisfies the optimality
  (KKT) conditions — at moderate tails, at heavy tails, and at long horizons. (Details in "How to defend
  it" below.)

Ready-to-paste methods sentence:

> *Grants are allocated by exact constrained optimization of the expected-output objective —
> water-filling within a round and projected-gradient ascent across the planning horizon — replacing an
> earlier discrete greedy fill. This removes a discretization artifact that inflated the
> forward-vs-myopic contrast at isolated budget values, and yields the exact optimum (verified against a
> fine-grained greedy limit and the KKT conditions).*

---

## What survives vs. what collapses (the map)

"Advantage" = mean forward-minus-myopic publications, grant-signal contrast (S8−S5), 200 seeds.

| Sweep | Old (greedy, n_steps=50) | New (exact) | Verdict |
|---|---|---|---|
| **horizon_growth** (T × ε) — **the headline** | up to +1.66 | up to **+1.99** | **HOLDS, cleaner** |
| alpha_regime | up to +0.25 | up to +0.26 | HOLDS |
| horizon_growth (pubs contrast P) | up to +3.30 | up to +3.10 | HOLDS |
| tail_map (heavy tails) | +0.50 | ~0.006 | **artifact → gone** |
| regime_map | +0.22 | ~0.006 | **artifact → gone** |
| info_value (ε≈0) | +0.22 | ~0.002 | **artifact → gone** |
| horizon_noise | +0.25 | ~0.03 | **artifact → gone** |
| horizon_scale (T × b, low ε) | +0.39 | ~0.06 | **artifact → gone** |
| funder_scale | +0.15 | ~0.006 | **artifact → gone** |
| signal_value / signal_precision | small | small | ~unchanged (both ~0) |
| **horizon_long** (T=5–10, strong ε) | — | *pending re-run* | expected to HOLD (strong ε) |

The unifying rule: **advantage tracks ε (compounding).** Big where ε and T are large; ~0 where ε is
small — regardless of the other knobs. The collapsed sweeps are the ones run at the base ε=0.1 (weak
compounding) or ε≈0.

---

## The reframed story (your new narrative arc)

**Old framing (implied):** forward planning helps across many regimes.
**New framing (true and sharper):** forward planning helps *because and when knowledge compounds*, and
the benefit accumulates over the horizon.

Narrative beats:

1. **Setup.** A funder allocates a fixed budget over T rounds. Myopic = maximize this round. Forward =
   plan the whole horizon, anticipating how today's funding raises tomorrow's productivity (compounding).
2. **Core finding.** When knowledge compounds, forward planning front-loads to let the gains compound
   longer, and this beats myopic by a margin that **grows with both the compounding rate ε and the
   horizon T** — up to ~2 extra publications at strong ε and T=5. (Fig 1, the clean ε×T surface.)
3. **Mechanism.** The advantage is exactly the value of anticipated compounding: at ε=0 it is provably
   zero (forward = myopic), and it rises monotonically in ε.
4. **Boundary / honesty.** Where compounding is weak or absent, forward planning gives essentially no
   advantage. (This is where careful, exact optimization matters — coarse optimization can manufacture a
   spurious advantage; we show it disappears under exact optimization.)
5. **Methodological contribution.** Using an exact allocator both sharpens the main result and prevents a
   class of discretization artifacts — a cautionary, reusable point for this modeling literature.

This is a stronger paper than "forward helps everywhere," because the effect now has a clean mechanism
(compounding), a clean monotonic figure, and an honest boundary.

---

## Section-by-section guidance for the draft

- **Abstract / intro:** lead with the compounding×horizon result. Drop any claim that forward helps
  broadly across tails/noise/scale.
- **Main result (Fig 1):** the ε×T surface. Emphasize monotonicity in both axes and the ε=0 baseline.
- **Mechanism:** the ε=0 ⇒ no-advantage result is your cleanest causal anchor. State it.
- **Robustness / secondary sweeps:** reframe from "advantage persists across regimes" to "advantage is
  governed by compounding; other factors modulate little." Report the near-zero contrasts honestly — they
  *support* the mechanism (no compounding ⇒ no advantage).
- **Methods:** the exact-allocator paragraph above.
- **Limitations / methods note:** one paragraph on the discretization artifact and why exact optimization
  matters — this is a genuine methodological contribution, not an embarrassment.
- **Things to delete:** any sentence asserting a forward advantage in tail_map, funder_scale, regime_map,
  info_value, horizon_noise, or low-ε horizon_scale.

---

## Figure guide

Regenerated from the exact-allocator data into `figures/smooth_preview/` (your published `figures/` are
untouched). Compare each against the old one before swapping.

- **fig1_horizon_phase_diagram** — the headline. Should look cleaner (monotonic). **Keep, update numbers.**
- **fig2_signal_value_law, fig3_signal_precision** — contrasts were ~0 before and after; check the story
  still reads. Likely **soften** claims.
- **fig4_correlation_robustness** — modest, roughly holds.
- **fig5_seed_value** — seed-value contrast; check magnitude dropped.
- **fig6_tail_map** — **biggest change**: the forward-advantage heat largely goes to ~0. Reframe as
  "advantage is tail-independent / near-zero without compounding."
- **fig7_information_value** — at ε≈0, advantage ~0. **Reframe** as confirming the mechanism.
- **fig8_horizon_long** — *pending the final sweep* (T=5–10). This is the one to watch for the "does the
  advantage keep growing or saturate" point.

---

## Presentation talking points (5 beats)

1. "Funders act under uncertainty and over time. Does planning ahead beat funding the best bet each
   round?"
2. "Yes — **but only when knowledge compounds.** The payoff to planning is the value of anticipated
   compounding." (Show Fig 1 ε×T surface.)
3. "It grows with the horizon: ~2 extra publications by round 5 at strong compounding."
4. "Where knowledge doesn't compound, planning ahead gives you nothing — and we can prove it (ε=0 ⇒
   forward = myopic)."
5. "A methodological aside: we found and removed a discretization artifact that had inflated some effects.
   Exact optimization gives a cleaner, more honest result." (Credibility-builder.)

---

## How to defend it (if challenged)

The exact allocator is validated, not assumed:

- **Optimality:** it reaches a *higher* objective value than a very fine version of the old allocator
  (n_steps=800) everywhere tested, and satisfies the KKT optimality conditions (all funded slots share a
  common marginal to ~1e-5).
- **Regimes tested:** moderate tails (T=2), heavy tails (k_shape=1.3, T=2), and long horizon (T=5) — it
  holds in all.
- **The collapses are real, not planner failure:** for the collapsed cells, the exact allocator agrees
  with the *fine-grained* old allocator (n_steps=800), while only the coarse (n_steps=50) old allocator
  was inflated. Example (heavy-tail cell): coarse-old +0.50, fine-old +0.03, exact +0.005.
- **Positive control:** the information-value sweep runs at ε≈0 and correctly gives ~0 advantage — the
  exact prediction of the mechanism.

Scripts: `tests/gonogo_smooth.R`, `tests/validate_smooth_supp.R`, `tests/diff_smooth_vs_fixed.R`,
`tests/waterfill_round_t.R`, `tests/plan_forward_smooth.R`.

---

## What's still pending

- **`horizon_long` re-run (T=5–10 at ε ∈ {0.3, 0.85}).** Slowest sweep; finishing now. It extends the
  headline to longer horizons (does the advantage saturate?). Expected to **hold** because it's the
  strong-compounding regime where the effect is real. I'll slot its numbers and fig8 in when it lands.
- **Decisions for you (not done automatically):** whether to make the exact allocator the model default,
  re-point the canonical data folder to the new run, and edit `RESULTS.md` in place. All staged; one word
  and I'll do them.

---

## File map (where everything is)

- **New results:** `sweep_results/T_run_smooth/` (exact allocator). Old preserved at
  `sweep_results/T_run_fixed/` (greedy, superseded).
- **New figures:** `figures/smooth_preview/`. Old published figures untouched in `figures/`.
- **Cell-by-cell diff:** `tests/diff_smooth_vs_fixed.log` (run `Rscript tests/diff_smooth_vs_fixed.R`).
- **The allocator code:** `model.R` (search `waterfill_round_t`, `plan_forward_smooth`); selected via
  `allocator = "smooth"`. Default is still `"greedy"` (bit-identical to the old model) until you flip it.
- **Full technical handoff:** `docs/SMOOTH_ALLOCATOR_HANDOFF.md`.
