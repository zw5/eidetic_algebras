# Findings from actual code restoration and independent verification

## Fi23: a false full-matrix assumption in the old runner

The current paper prints

    P_ram(Fi23) = M_92(Q) + M_3(Q) + Q^3
    dim A = 8476, dim Z(A) = dim A' = 5.

The restored independent certificate instead proves

    P_ram(Fi23) = M_91(Q) + M_3(Q) + Q^4
    dim A = 8294, dim Z(A) = dim A' = 6.

All module multiplicities are one. The represented dimension is
91 + 3 + 4 = 98. The algebra dimension is 91^2 + 3^2 + 4 = 8294.

An exact extra invariant line is spanned by

    v = delta_12H - (1/3) delta_12L - (1/4) delta_12M + (1/12) delta_12O.

Every ramified prime generator and its weighted adjoint preserves this line.
For primes 2 and 3 the eigenvalue is zero; for the other prime divisors of the
exponent it is one. The vector is orthogonal to the identity-class cyclic
module. Exact rational row reduction gives dimension 91 for that module, and
nonzero 91-vector minors replay at both modular primes.

The recovered historical algorithm called the number of stable weighted
refinement cells the principal full-matrix degree. Its 92-cell partition
is equitable, but that alone does not prove fullness. It reproduced the
printed erroneous number until the independent cyclicity test was added.
The corrected runner moves v into the seven-dimensional residual carrier and
certifies that entire residual algebra over Q, including its central splitting.
The legacy direct Fi23 shortcut now raises an error instead of silently
returning the unverified row.

Required paper changes include:

- Sporadic census: replace the Fi23 row by (98, 8294, 6, 6),
  M_91(Q) + M_3(Q) + Q^4.
- Good-label table: replace the single trivial 92-dimensional block by a
  trivial 91-dimensional block and a trivial 1-dimensional block. The rational
  class count remains 92. Other characters in that row are unchanged.
- Certification text: explain identity-class cyclicity and the Fi23 residual
  line, rather than deriving fullness from equitable refinement alone.

The other 25 sporadic decomposition rows agree with the regenerated results.
The Monster decomposition, 186-vector certificate, centers census, and Schur
index-one conclusion all survive the independent checks.

Evidence: `data/fi23_independent_check.json`, `data/fi23_orbit_mod.json`,
`data/principal_certificates.json`, `data/census.json`,
`data/good_label_characters.json`, and `data/paper_comparison.json`.

## PSL(2,29): an omitted example

The Families section says it lists all odd-prime rational-class sizes for
3 <= q <= 32, but omits q = 29. Its own general formula gives torus orders
14 and 15. Element orders d = 7 and d = 14 give phi(d)/2 = 3.
Fresh GAP power maps confirm both rational classes of size 3.

Add PSL(2,29): d = 7,14, r = 3 to that list. The general theorem is unaffected.
See `data/psl_examples.json` and `psl_examples.py`.

## Release status

The code is restored and reproducible. The manuscript/PDF still needs these
corrections and a matching code-availability statement; this restoration did
not silently edit the mathematical text. Read computational verification as
verification of the corrected outputs, not endorsement of the unchanged PDF.
