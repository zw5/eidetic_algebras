# Eidetic power-map algebras

**Eidetic power-map algebras: balanced kernels, rational classes, and ramified completion** — Simon Velez.

- [`paper/eidetic_algebras.pdf`](paper/eidetic_algebras.pdf) — the manuscript (38 pp.).
  Sources in `paper/`; build with `latexmk -pdf eidetic_algebras.tex` from that directory.
- [`code/`](code/README.md) — computational supplement: GAP/CTblLib exports, rational
  word-closure runners, the 186-vector Monster certificate, all 26 sporadic computations,
  and a reproduction script with tests.

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r code/requirements.txt
.venv/bin/python code/reproduce.py
```
