# Eidetic power-map algebras

Active source: `eidetic_algebras.tex`.  Build from this directory:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error eidetic_algebras.tex
```

Layout after the 2026-09-16 restructure (38 pages):

- `eidetic_algebras.tex` — preamble, abstract, introduction (definition, D8/Q8,
  Monster statement, Theorems A–D, related constructions, outline), open problems,
  bibliography.
- `sections/algebra.tex` — §2 weighted power systems, transfer calculus,
  semisimplicity and bicommutant, balanced kernels, saturation criteria.
- `sections/goodlabels.tex` — §3 Galois action, fields of character values,
  Euler product, ramified local factors, block characters.
- `sections/sporadic.tex` — §4 sporadic census, good-label characters, rank-one
  saturation, Monster theorem, symmetries invisible to the character table.
- `sections/families.tex` — §5 prime exponent, cyclic prime powers, dihedral,
  PSL(2,q) rational classes, comparison panel.
- `sections/certification.tex` — §6 computational verification.
- `sections/dynamics.tex` — Appendix A functional-graph zeta packet.
- `sections/*_rows.tex` — table bodies.

Material cut in the restructure, and a full pre-restructure snapshot, are in
`Documents/beam_rh/eidetic_stripped/`.

The computational supplement is restored in [`ancillary/`](ancillary/README.md),
with fresh GAP/CTblLib exports, exact rational runners, block characters, and
replayable certificates. See [`ancillary/PAPER_DISCREPANCIES.md`](ancillary/PAPER_DISCREPANCIES.md):
independent reproduction corrected the Fi23 decomposition and found an omitted
PSL(2,29) example. Both corrections were verified independently from the raw
exports and are incorporated in the manuscript and PDF as of 2026-09-16
(Fi23 row: M91(Q) ⊕ M3(Q) ⊕ Q⁴, dim 8294, center 6). Reproduce with
`.venv/bin/python ancillary/reproduce.py`; `ancillary/compare_paper.py` reports
26/26 rows matching.
