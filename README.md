# Eidetic power-map algebras

This repo contains a paper I developed with the assistance of GPT 5.5. It describes Eidetic algebras of finite groups, which are: given a finite group G, the eidetic algebra is created by taking the rational class functions with the inner product weighted by class sizes, and the power maps $\Psi_p f(C) = f(C^p)$ for the primes $p$ dividing the exponent.

In this paper we prove the algebra is semisimple and has its own bicommutant. Its commutant is the space of balanced kernels on the weighted power graph, so the algebra is only identical to the full matrix algebra IFF the balanced kernels are scalars. 

It is not a function of the character table. $D_8$ and $Q_8$ have the same character table, and their algebras have dimensions 17 and 10.

For primes $p$ that are indivisible by the exponent, $\Psi_p$ is commutative with the algebra and can be expressed as the Frobenius at $p$ on the étale algebra of character-value fields. The product $\prod_p \det(I - \Psi_p p^{-s})^{-1}$ is the Dedekind zeta function of that algebra, and the paper shows the factors at ramified primes as well.

We compute the Eidetic algebras of all 26 sporadic simple groups, and find every simple factor has Schur index one, with all rational centers being either $\mathbb{Q}$, $\mathbb{Q}(\zeta_3)$, and one group having $\mathbb{Q}(\zeta_5)$. The groups with a non-rational center are exactly the six pariah groups. 

For the Monster group, 
$$\mathcal{P}{\mathrm{ram}}(\mathbb{M}) \cong M{170}(\mathbb{Q}) \oplus M_5(\mathbb{Q}) \oplus M_3(\mathbb{Q}) \oplus M_2(\mathbb{Q})^{4} \oplus \mathbb{Q}^{8},$$ 
of dimension 28958. The algebra is the centralizer of the automorphism group of the Monster's weighted power system, which is elementary abelian of rank 14, with thirteen generators as Galois conjugations of classes. The fourteenth is the permutation $(16B\ 16C)(32A\ 32B)$. 

$Fi_{23}$ also shows a different kind of hidden symmetry. Its algebra has a one-dimensional block spanned by a vector with four distinct coefficients on four classes of order 12.

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
