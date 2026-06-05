# E156-PROTOCOL — Meta-Analysis Audit Tool

- **Project:** Meta-analysis-audit-tool (GitHub repo `Meta-analysis-audit-tool`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `metaanlysisaudittool.html` dump)
- **Type:** single-file offline browser tool + Node-testable engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline**: removed the Google Fonts CDN `<link>` (the only
  network resource); the page now loads nothing external. (Remaining
  `http://www.w3.org/2000/svg` strings are XML namespace identifiers, not fetches.)
- Extracted the statistical core **verbatim** into a pure `engine.js` (single
  source of truth); wired `<script src="engine.js">` and deleted the inline
  duplicates. The DOM/audit/render layer stays inline.
- Added `tests.js` (58 assertions, all passing) with a hand-derived two-study
  DerSimonian–Laird worked example.
- **Fixed a correctness bug** in `applyHKSJ`: the Hartung–Knapp variance had no
  floor, so a homogeneous set (`Q* < k−1`) produced a CI narrower than the
  ordinary random-effects CI. Now floored at `max(1, Q*/(k−1))`.
- Added Pages scaffold (`.nojekyll`, `.gitignore`, README); renamed
  `metaanlysisaudittool.html` → `index.html`.

## Body (E156 draft — CURRENT BODY)

When a meta-analysis reports a pooled effect, can its study-level numbers be
re-derived to expose discrepancies between what was published and what the data
support? This browser tool ingests per-study effect sizes plus the paper's
reported pooled estimate, confidence interval, I², τ², Q and weights. It pools
log-scale effects under fixed-effect, DerSimonian–Laird, REML, Paule–Mandel and
Empirical-Bayes τ² estimators, optionally with a Hartung–Knapp interval, and
contrasts each against the published figures. Findings flag sign flips,
implausibly narrow confidence intervals, heterogeneity mismatches, weight errors,
Egger funnel asymmetry, and leave-one-out patterns that hint at an undisclosed
exclusion. A revival audit found and fixed a Hartung–Knapp variance floor whose
absence let homogeneous data yield an anti-conservatively narrow interval, and
locked the core behind a 58-assertion suite checked against hand computation. The
result is a transparent recomputation aid that surfaces leads for closer scrutiny.
It detects discrepancies, not intent, and is not a verdict on any specific review.

SUBMITTED: [ ]
