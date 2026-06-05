# Meta-Analysis Audit Tool

A single-file, **fully offline** browser tool that re-derives a published
meta-analysis from its study-level data and flags discrepancies between what was
*reported* and what the numbers *actually support*. You paste the per-study
effect sizes (and, optionally, the paper's reported pooled effect, CI, I², τ²,
Q, and weights); the tool recomputes everything under several estimators and
raises critical / warning findings where the published figures disagree.

**Live app:** open `index.html` (or the GitHub Pages link). No build step, no
network requests, no external CDN — system fonts are used.

## Layout

```
index.html   single-file UI (loads engine.js)
engine.js    pure statistical core — runs in Node and the browser
tests.js     Node test harness, 58 assertions
LICENSE      Apache-2.0
```

## What it checks

- **Pooled effect & CI** recomputed under Fixed-Effect, DerSimonian–Laird, REML,
  Paule–Mandel and Empirical-Bayes τ² estimators (each optionally with the
  Hartung–Knapp–Sidik–Jonkman small-sample CI), then matched against the
  published estimate. Ratio measures (OR/RR/HR) are pooled on the **log scale**
  and back-transformed.
- **Sign flips** — published effect pointing the opposite way to the recomputed one.
- **Confidence-interval width** — a published CI substantially narrower than the
  data support (a classic symptom of a wrong SE or an undisclosed HKSJ).
- **Heterogeneity** — I², τ² (against a Q-profile 95% CI), and Q.
- **Study weights** — per-study weight mismatches > 5%.
- **Undisclosed exclusion** — leave-one-out re-pooling that reproduces the
  published number, suggesting a silently dropped study.
- **Publication-bias diagnostics** — Egger's regression test for funnel-plot
  asymmetry, plus a funnel plot.
- **Prediction interval** — `t_{k-2}` PI compared against the reported one.

## Statistical core (`engine.js`)

| Function | What it does |
|---|---|
| `metaFE / metaDL / metaREML / metaPM / metaEB` | pool log-scale effects; τ² by FE/DerSimonian–Laird/REML/Paule–Mandel/Empirical-Bayes |
| `applyHKSJ` | Hartung–Knapp–Sidik–Jonkman CI with a variance floor and `t_{k-1}` |
| `calcPredictionInterval` | 95% prediction interval (`t_{k-2}`) |
| `tau2CI` | Q-profile 95% CI for τ² |
| `eggersTest` | Egger's regression test for funnel asymmetry |
| `leaveOneOut` | leave-one-out re-pooling |
| `pnorm / pchisq / qchisq / pt / qt / lgamma / gammainc / betacf` | distribution helpers |

τ² is `max(0, (Q − (k−1)) / C)` with `C = Σw − Σw²/Σw`; I² is
`max(0, (Q − (k−1)) / Q)`.

## Fixes applied during revival (2026-06-05)

1. **Offline:** removed the Google Fonts `<link>` (the only network resource);
   the page now loads nothing external (the remaining `http://www.w3.org/2000/svg`
   strings are SVG XML namespace identifiers, never fetched).
2. **Single source of truth:** extracted the pure statistical functions
   **verbatim** into `engine.js`, wired `<script src="engine.js">`, and deleted
   the inline duplicates. The DOM/audit/render layer stays inline and consumes
   the engine globals.
3. **Correctness fix — HKSJ variance floor.** The original Hartung–Knapp SE was
   `sqrt(Q*/((k−1)·Σw))` with **no floor**, so when `Q* < k−1` the adjusted CI
   *narrowed below* the ordinary random-effects CI (anti-conservative — the
   documented HKSJ failure mode). It now uses the standard multiplier
   `max(1, Q*/(k−1))`, locked by a test. The `t_{k-1}` critical value was already
   correct and is unchanged.
4. Renamed `metaanlysisaudittool.html` → `index.html`; added `.nojekyll`,
   `.gitignore`, this README, and the E156 protocol.

The remaining pooling math (τ², I², Q, weights, the REML/PM/EB iterations,
Egger's test, Q-profile τ² CI, the `t_{k-2}` PI) was verified against hand
computation and left unchanged.

## Tests

```
node tests.js
# 58 passed, 0 failed
```

Coverage includes the distribution helpers (`Φ(0)=0.5`, `Φ(1.96)=0.975`,
`qt(0.975, df)` vs the t-table, χ² CDF reference points), a fully hand-derived
two-study DerSimonian–Laird example (τ²≈0.006579, pooled OR≈0.7478, I²≈20.83%),
edge cases (k=1 passthrough, two-identical ⇒ τ²=I²=0, empty ⇒ null), the HKSJ
floor property (never narrower than RE; widens under heterogeneity; `t_{k-1}`),
Egger nulls for k<3 and symmetric funnels, leave-one-out cardinality, and a
prediction interval wider than the CI.

## Caveats

DerSimonian–Laird underestimates τ² for small *k* (REML / Paule–Mandel are
preferred for k<10); all five estimators are reported side-by-side. The
"best-matching estimator" heuristic minimises distance to the published values
and is a diagnostic aid, not a claim about what the original authors actually
ran. Egger's test has low power for k<10 and is degenerate when every study
shares the same SE. Treat findings as leads for closer inspection of a
meta-analysis, not as a verdict. Apache-2.0 licensed.
