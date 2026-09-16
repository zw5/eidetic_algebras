# Eidetic power-map algebra computational supplement

Computational supplement for the paper, regenerated on 2026-09-16 with fresh
GAP/CTblLib exports and the certificate methods described in the paper.

**Read [PAPER_DISCREPANCIES.md](PAPER_DISCREPANCIES.md).** Computation found an
incorrect Fi23 decomposition and an omitted PSL(2,29) example in the current
manuscript. The corrected data are in `data/census.json`; the uncorrected
historical algorithm's reproduced output is retained only for comparison.

## Reproduce

From the manuscript directory, using a Python environment with pip:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r code/requirements.txt
.venv/bin/python code/reproduce.py
```

The default command uses the included exact GAP exports. It regenerates the
Fi23 counterexample, sporadic decomposition and principal certificates, block
characters, ramified factors, Monster saturation certificate, dynamics packets,
family ladder, PSL examples, and 31-group comparison census; then it runs tests.
There is no dependency on chat logs or other research repositories. Do not run
certificate scripts with Python's `-O` option, which disables assertions.

To refresh the external inputs as well, install GAP 4.14.0 with CTblLib 1.3.9
and its dependencies, then run:

```sh
.venv/bin/python code/reproduce.py --gap /absolute/path/to/gap
```

A local optional bootstrap is provided:

```sh
bash code/bootstrap_gap.sh
.venv/bin/python code/reproduce.py --gap "$PWD/.tools/gap-4.14.0/gap"
```

The bootstrap downloads upstream GAP/packages and builds GAP locally; it does
not install system-wide. The frozen exports identify GAP and CTblLib versions.
`requirements-lock.txt` records the Python packages used in the restoration;
`requirements.txt` is the portable dependency list. Validation used Python 3.14.

Run just the tests:

```sh
cd code
../.venv/bin/python -m pytest -q
```

## Entrypoints and outputs

| Entrypoint | Output / purpose |
|---|---|
| `restore_census.py` | `data/census.json`, `data/principal_certificates.json`: corrected exact 26-group census |
| `monster_certificate.py` | `data/monster_certificate.json`: 186 words, nonzero minors at two primes, symmetry generators, root profiles, scalar separation |
| `characters.py` | `data/good_label_characters.json`: exact central actions, conductors, quadratic discriminants |
| `local_factors.py` | `data/ramified_local_factors.json`: all 173 bad-prime slots by two independent routes, 80 abelian controls |
| `dynamics.py` | `data/dynamics.json`: all 258 prime-power maps and 26 distinct unlabelled packets |
| `families.py` | `data/families.json`: 21 dihedral/quaternion rungs and 8 cyclic controls |
| `psl_examples.py` | `data/psl_examples.json`: all 17 prime powers from 3 through 32 |
| `investigate_fi23.py` | exact 91-dimensional orbit and explicit missing invariant vector |
| `compare_paper.py` | `data/paper_comparison.json`: current printed dimensions versus computed dimensions |
| `export_enriched.py`, `export_psl.g` | regenerate class names, weights, power maps, unit actions, stored table automorphisms |
| `runners/experiments/finite_group_adams_star_census_20260724/` | recovered general rational closure runner, exporter, tests and 31-group output |
| `runners/experiments/sporadic_adams_star_census_20260724/` | recovered residual closure/center factorization routines, extended to retain exact idempotents |

`data/census.json` is authoritative for the restored sporadic results. The
older `.../sporadic_adams_star_census_20260724/output/sporadic_adams_star_census.json`
is explicitly historical and includes the false Fi23 row; it is not consumed
by the corrected census, characters, or Monster verifier.

## Certificate boundaries

The main matrix blocks are certified by rank-one generation from the identity
class. A nonzero modular minor supplies a lower bound over Q; exact invariance
of the complementary rational carrier supplies the upper bound. For Fi23 the
92-cell equitable partition does not give a 92-dimensional irreducible block:
the extra exact invariant vector is moved to the residual carrier.

Residual word bases are lifted to rational matrices, checked for rank and
closure over Q, and split by rational central idempotents. Center fields and
Schur index one follow from their ideal dimensions, module dimensions, and
the coprime split-degree/multiplicity test. The two modular primes are
1,000,003 and 1,000,033; bad denominators are rejected.

The restored local-factor calculation exports the full local unit actions and
therefore verifies all 173 slots. It does not reproduce the historical omission
of five slots caused by sparse input maps. The optional floating-point
Laplacian diagnostics in recovered legacy outputs are not used as proof.

This package reproduces finite computations. It does not certify every
noncomputational proof or the paper's literature/novelty claims.

Data source: [GAP](https://www.gap-system.org/) and
[CTblLib](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/).
See [RECOVERY_PROVENANCE.md](RECOVERY_PROVENANCE.md) for recovered versus
reconstructed source provenance.
