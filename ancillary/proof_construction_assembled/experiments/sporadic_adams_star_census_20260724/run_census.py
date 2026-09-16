"""Compute Adams-star invariants for all 26 sporadic simple groups.

GAP's Character Table Library supplies exact conjugacy-class power maps,
so no permutation or matrix realization of the large groups is required.
The stable diagonal quotient gives the principal full matrix block.  The
remaining within-cell carrier is small enough for exact rational closure.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import gcd
from pathlib import Path
import shutil
import subprocess
import sys

import flint
from flint import fmpq, fmpq_mat
import numpy as np
import sympy as sp

from proof_construction_assembled.experiments.finite_group_adams_star_census_20260724.run_census import (
    MODULAR_PRIMES,
    adjoint_matrix,
    center_dimension,
    commutant_dimension,
    component_count,
    complex_type_solutions,
    laplacian_certificate,
    pullback_matrix,
    star_closure,
)


HERE = Path(__file__).resolve().parent
REPOSITORY_ROOT = HERE.parents[2]
DEFAULT_GAP = Path(shutil.which("gap") or "gap")
DEFAULT_DATA = HERE / "output" / "sporadic_power_maps.json"
DEFAULT_OUTPUT = HERE / "output" / "sporadic_adams_star_census.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance(gap_data: Path) -> dict[str, object]:
    source = json.loads(gap_data.read_text())
    return {
        "python": sys.version.split()[0],
        "python_flint": flint.__version__,
        "flint": flint.__FLINT_VERSION__,
        "numpy": np.__version__,
        "sympy": sp.__version__,
        "gap": source.get("software", {}).get("gap"),
        "ctbllib": source.get("software", {}).get("ctbllib"),
        "raw_input_sha256": sha256(gap_data),
        "runner_sha256": sha256(Path(__file__)),
        "gap_export_sha256": sha256(HERE / "export_sporadic_tables.g"),
    }


def export_gap_data(gap: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        [str(gap), "-q", str(HERE / "export_sporadic_tables.g")],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    output.write_text(json.dumps(payload, indent=2) + "\n")


def stable_diagonal_partition(class_sizes: list[int], maps: list[list[int]]):
    """Compute a weighted equitable partition; fullness needs a cyclic certificate."""

    size = len(class_sizes)
    zero_based = [[target - 1 for target in mapping] for mapping in maps]
    cells = [0] * size
    history = []
    while True:
        block_count = max(cells) + 1
        signatures = []
        for target in range(size):
            signature: list[object] = [cells[target]]
            for mapping in zero_based:
                signature.append(cells[mapping[target]])
                incoming = [0] * block_count
                for source, image in enumerate(mapping):
                    if image == target:
                        incoming[cells[source]] += class_sizes[source]
                signature.extend(
                    Fraction(value, class_sizes[target]) for value in incoming
                )
            signatures.append(tuple(signature))
        lookup: dict[tuple[object, ...], int] = {}
        refined = []
        for signature in signatures:
            lookup.setdefault(signature, len(lookup))
            refined.append(lookup[signature])
        history.append(len(lookup))
        if refined == cells:
            return cells, history
        cells = refined


def partition_cells(labels: list[int]) -> list[list[int]]:
    cells: dict[int, list[int]] = defaultdict(list)
    for index, label in enumerate(labels):
        cells[label].append(index)
    return [cells[index] for index in sorted(cells)]


def quotient_connected(labels: list[int], maps: list[list[int]]) -> bool:
    count = max(labels) + 1
    parent = list(range(count))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for mapping in maps:
        for source, target in enumerate(mapping):
            union(labels[source], labels[target - 1])
    return len({find(value) for value in range(count)}) == 1


def exact_pullback(mapping: list[int]) -> fmpq_mat:
    size = len(mapping)
    return fmpq_mat([
        [1 if column == mapping[row] - 1 else 0 for column in range(size)]
        for row in range(size)
    ])


def exact_adjoint(matrix: fmpq_mat, class_sizes: list[int]) -> fmpq_mat:
    size = matrix.nrows()
    return fmpq_mat([
        [
            fmpq(int(matrix[column, row]) * class_sizes[column], class_sizes[row])
            for column in range(size)
        ]
        for row in range(size)
    ])


def residual_basis(
    class_sizes: list[int], cells: list[list[int]]
) -> tuple[fmpq_mat, list[int]]:
    """Return a Plancherel-orthogonal within-cell basis and pivot rows."""

    size = len(class_sizes)
    columns: list[list[fmpq]] = []
    pivots: list[int] = []
    for cell in cells:
        reference = cell[0]
        for index in cell[1:]:
            column = [fmpq(0) for _ in range(size)]
            column[index] = fmpq(1, class_sizes[index])
            column[reference] = fmpq(-1, class_sizes[reference])
            columns.append(column)
            pivots.append(index)
    if not columns:
        return fmpq_mat(size, 0), pivots
    return fmpq_mat([
        [columns[column][row] for column in range(len(columns))]
        for row in range(size)
    ]), pivots


def residual_restriction(
    matrix: fmpq_mat,
    basis: fmpq_mat,
    pivots: list[int],
    class_sizes: list[int],
) -> fmpq_mat:
    image = matrix * basis
    rank = len(pivots)
    restricted = fmpq_mat([
        [class_sizes[pivots[row]] * image[pivots[row], column]
         for column in range(rank)]
        for row in range(rank)
    ])
    if basis * restricted != image:
        raise ArithmeticError("stable diagonal residual is not invariant")
    return restricted


def matrix_mod_prime(matrix: fmpq_mat, prime: int) -> np.ndarray:
    if any(int(x.denominator) % prime == 0 for x in matrix.entries()):
        raise ValueError("bad reduction prime divides a denominator")
    return np.asarray([
        [
            int(matrix[row, column].numerator)
            * pow(int(matrix[row, column].denominator) % prime, prime - 2, prime)
            % prime
            for column in range(matrix.ncols())
        ]
        for row in range(matrix.nrows())
    ], dtype=np.int64)


def exact_closure_certificate(
    exact_generators: tuple[fmpq_mat, ...],
    modular_basis_matrices: tuple[np.ndarray, ...],
    expected_dimension: int,
) -> dict:
    """Lift a modular word basis and certify closure over Q."""

    from proof_construction_assembled.experiments.finite_group_adams_star_census_20260724.run_census import (
        ModularBasis,
    )

    size = exact_generators[0].nrows()
    prime = MODULAR_PRIMES[0]
    modular_basis = ModularBasis(size * size, prime)
    identity_exact = fmpq_mat([
        [1 if row == column else 0 for column in range(size)]
        for row in range(size)
    ])
    identity_modular = np.eye(size, dtype=np.int64)
    exact_basis: list[fmpq_mat] = []
    frontier: list[tuple[np.ndarray, fmpq_mat]] = []
    for modular, exact in (
        (identity_modular, identity_exact),
        *zip(modular_basis_matrices, exact_generators, strict=True),
    ):
        if modular_basis.add(modular):
            exact_basis.append(exact)
            frontier.append((modular, exact))
    while frontier:
        new = []
        for modular_left, exact_left in frontier:
            for modular_right, exact_right in zip(
                modular_basis_matrices, exact_generators, strict=True
            ):
                modular_product = modular_left @ modular_right % prime
                exact_product = exact_left * exact_right
                if modular_basis.add(modular_product):
                    exact_basis.append(exact_product)
                    new.append((modular_product, exact_product))
        frontier = new
    exact_rank = fmpq_mat(
        [list(matrix.entries()) for matrix in exact_basis]
    ).rank()
    closure_rows = [
        list((left * right).entries())
        for left in exact_basis
        for right in exact_generators
    ]
    closure_rank = fmpq_mat(
        [list(matrix.entries()) for matrix in exact_basis] + closure_rows
    ).rank()
    return {
        "word_basis_size": len(exact_basis),
        "exact_word_basis_rank": exact_rank,
        "exact_closure_rank": closure_rank,
        "expected_dimension": expected_dimension,
        "exact_over_Q": (
            len(exact_basis)
            == exact_rank
            == closure_rank
            == expected_dimension
        ),
    }


def exact_word_basis(
    exact_generators: tuple[fmpq_mat, ...],
    modular_generators: tuple[np.ndarray, ...],
) -> list[fmpq_mat]:
    from proof_construction_assembled.experiments.finite_group_adams_star_census_20260724.run_census import (
        ModularBasis,
    )

    size = exact_generators[0].nrows()
    prime = MODULAR_PRIMES[0]
    modular_basis = ModularBasis(size * size, prime)
    identity_exact = fmpq_mat([
        [1 if row == column else 0 for column in range(size)]
        for row in range(size)
    ])
    exact_basis: list[fmpq_mat] = []
    frontier: list[tuple[np.ndarray, fmpq_mat]] = []
    for modular, exact in (
        (np.eye(size, dtype=np.int64), identity_exact),
        *zip(modular_generators, exact_generators, strict=True),
    ):
        if modular_basis.add(modular):
            exact_basis.append(exact)
            frontier.append((modular, exact))
    while frontier:
        new = []
        for modular_left, exact_left in frontier:
            for modular_right, exact_right in zip(
                modular_generators, exact_generators, strict=True
            ):
                modular_product = modular_left @ modular_right % prime
                exact_product = exact_left * exact_right
                if modular_basis.add(modular_product):
                    exact_basis.append(exact_product)
                    new.append((modular_product, exact_product))
        frontier = new
    return exact_basis


def to_sympy(value: fmpq) -> sp.Rational:
    return sp.Rational(int(value.numerator), int(value.denominator))


def central_isotypic_dimensions(
    exact_basis: list[fmpq_mat],
    exact_generators: tuple[fmpq_mat, ...],
    expected_center_dimension: int,
) -> dict:
    """Use a generic exact central element to recover d_j m_j."""

    size = exact_generators[0].nrows()
    constraints = []
    for generator in exact_generators:
        commutators = [
            element * generator - generator * element
            for element in exact_basis
        ]
        for row in range(size):
            for column in range(size):
                constraints.append([
                    to_sympy(commutator[row, column])
                    for commutator in commutators
                ])
    nullspace = sp.Matrix(constraints).nullspace()
    if len(nullspace) != expected_center_dimension:
        raise ArithmeticError("exact center nullspace has the wrong dimension")

    factor_receipt = None
    for attempt in range(1, 6):
        coefficients = [
            sum(
                sp.Integer((center_index + 1) ** attempt)
                * nullspace[center_index][basis_index]
                for center_index in range(len(nullspace))
            )
            for basis_index in range(len(exact_basis))
        ]
        central = sp.zeros(size)
        for coefficient, element in zip(coefficients, exact_basis, strict=True):
            if coefficient:
                central += coefficient * sp.Matrix([
                    [to_sympy(element[row, column]) for column in range(size)]
                    for row in range(size)
                ])
        polynomial = central.charpoly().as_poly()
        _, factors = sp.factor_list(polynomial)
        root_count = sum(factor.degree() for factor, _ in factors)
        if root_count == expected_center_dimension:
            variable = polynomial.gens[0]
            minimal_polynomial = sp.prod(
                factor.as_expr() for factor, _ in factors
            )
            factor_receipt = []
            identity = sp.eye(size)
            for factor, exponent in factors:
                quotient = sp.cancel(
                    minimal_polynomial / factor.as_expr()
                )
                inverse = sp.invert(
                    quotient, factor.as_expr(), domain=sp.QQ
                )
                idempotent_polynomial = sp.rem(
                    quotient * inverse, minimal_polynomial,
                    domain=sp.QQ,
                )
                coefficients = sp.Poly(
                    idempotent_polynomial, variable
                ).all_coeffs()
                idempotent = sp.zeros(size)
                for coefficient in coefficients:
                    idempotent = idempotent * central + coefficient * identity
                ideal_rows = []
                for element in exact_basis:
                    product_matrix = idempotent * sp.Matrix([
                        [
                            to_sympy(element[row, column])
                            for column in range(size)
                        ]
                        for row in range(size)
                    ])
                    ideal_rows.append([
                        entry for entry in product_matrix
                    ])
                ideal_dimension = sp.Matrix(ideal_rows).rank()
                factor_receipt.append({
                    "degree": factor.degree(),
                    "characteristic_exponent": exponent,
                    "defining_polynomial": str(factor.as_expr()),
                    "rational_simple_factor_dimension": ideal_dimension,
                    "idempotent": [[str(idempotent[i,j]) for j in range(size)] for i in range(size)],
                })
            dimensions = sorted(
                exponent
                for factor, exponent in factors
                for _ in range(factor.degree())
            )
            return {
                "generic_center_attempt": attempt,
                "central_matrix": [[str(central[i,j]) for j in range(size)] for i in range(size)],
                "factor_degrees_and_exponents": factor_receipt,
                "isotypic_dimensions_over_C": dimensions,
            }
    raise ArithmeticError(
        f"generic center did not separate all blocks: {factor_receipt}"
    )


def rational_factors_from_center(
    complex_type: tuple[tuple[int, int], ...],
    center_receipt: dict,
) -> list[dict] | None:
    """Recover rational matrix factors from exact central idempotent ideals."""

    factors = []
    recovered_complex_type = []
    for receipt in center_receipt["factor_degrees_and_exponents"]:
        center_degree = receipt["degree"]
        ideal_dimension = receipt["rational_simple_factor_dimension"]
        quotient = ideal_dimension // center_degree
        matrix_size = int(quotient**0.5)
        if (
            matrix_size * matrix_size != quotient
            or receipt["characteristic_exponent"] % matrix_size
        ):
            return None
        multiplicity = receipt["characteristic_exponent"] // matrix_size
        if gcd(matrix_size, multiplicity) != 1:
            return None
        recovered_complex_type.extend(
            [(matrix_size, multiplicity)] * center_degree
        )
        factors.append({
            "center_degree": center_degree,
            "center_defining_polynomial": receipt["defining_polynomial"],
            "center_field": identify_center_field(
                receipt["defining_polynomial"]
            ),
            "matrix_size_over_center": matrix_size,
            "module_multiplicity": multiplicity,
            "schur_index": 1,
        })
    if sorted(recovered_complex_type) != sorted(complex_type):
        return None
    return factors


def identify_center_field(defining_polynomial: str) -> str:
    variable = sp.Symbol("x")
    polynomial = sp.Poly(
        sp.sympify(defining_polynomial.replace("lambda", "x")),
        variable,
    )
    degree = polynomial.degree()
    if degree == 1:
        return "Q"
    if degree == 2:
        discriminant = int(sp.discriminant(polynomial.as_expr(), variable))
        factors = sp.factorint(abs(discriminant))
        squarefree = -1 if discriminant < 0 else 1
        for prime, exponent in factors.items():
            if exponent % 2:
                squarefree *= prime
        return f"Q(sqrt({squarefree}))"
    if degree == 4:
        cyclotomic = sp.Poly(
            variable**4 + variable**3 + variable**2 + variable + 1,
            variable,
        )
        from sympy.polys.numberfields import field_isomorphism

        if field_isomorphism(
            sp.CRootOf(polynomial.as_expr(), 0),
            sp.CRootOf(cyclotomic.as_expr(), 0),
            fast=False,
        ) is not None:
            return "Q(zeta_5)"
    return f"number_field_degree_{degree}"


def residual_algebra(
    class_sizes: list[int],
    maps: list[list[int]],
    cells: list[list[int]],
    extra_vector=None,
) -> dict:
    basis, pivots = residual_basis(class_sizes, cells)
    extra_left = None
    if extra_vector is not None:
        basis = fmpq_mat([[basis[i,j] for j in range(basis.ncols())] + [fmpq(extra_vector[i])] for i in range(basis.nrows())])
        rr, rank = basis.transpose().rref()
        pivots = [next(j for j in range(rr.ncols()) if rr[i,j]) for i in range(rank)]
        assert rank == basis.ncols()
        extra_left = fmpq_mat([[basis[i,j] for j in range(rank)] for i in pivots]).inv()
    residual_dimension = len(pivots)
    if residual_dimension == 0:
        return {
            "carrier_dimension": 0,
            "algebra_dimension": 0,
            "center_dimension": 0,
            "commutant_dimension": 0,
            "complex_type": [],
            "rational_wedderburn_factors": [],
            "rational_classification_certified": True,
            "rational_split_certified": True,
            "exact_certificate": {"exact_over_Q": True},
        }

    pullbacks = tuple(exact_pullback(mapping) for mapping in maps)
    adjoints = tuple(exact_adjoint(matrix, class_sizes) for matrix in pullbacks)
    def restrict(matrix):
        if extra_left is None:
            return residual_restriction(matrix, basis, pivots, class_sizes)
        image = matrix * basis
        small = extra_left * fmpq_mat([[image[i,j] for j in range(image.ncols())] for i in pivots])
        if basis * small != image:
            raise ArithmeticError("extended residual is not invariant")
        return small
    exact_generators = tuple(restrict(matrix) for matrix in pullbacks + adjoints)
    receipts = []
    reference = None
    first_modular_generators = None
    for prime in MODULAR_PRIMES:
        modular_generators = tuple(
            matrix_mod_prime(matrix, prime) for matrix in exact_generators
        )
        closure = star_closure(modular_generators, prime)
        receipt = {
            "prime": prime,
            "algebra_dimension": len(closure),
            "center_dimension": center_dimension(
                closure, modular_generators, prime
            ),
            "commutant_dimension": commutant_dimension(
                modular_generators, prime
            ),
        }
        receipts.append(receipt)
        key = (
            receipt["algebra_dimension"],
            receipt["center_dimension"],
            receipt["commutant_dimension"],
        )
        if reference is None:
            reference = key
            first_modular_generators = modular_generators
        elif key != reference:
            raise ArithmeticError("residual modular receipts disagree")
    algebra_dimension, center_dim, commutant_dim = reference
    types = complex_type_solutions(
        residual_dimension, algebra_dimension, commutant_dim, center_dim
    )
    certificate = exact_closure_certificate(
        exact_generators, first_modular_generators, algebra_dimension
    )
    if not certificate["exact_over_Q"]:
        raise ArithmeticError("exact residual closure certificate failed")
    exact_basis = exact_word_basis(
        exact_generators, first_modular_generators
    )
    center_receipt = central_isotypic_dimensions(
        exact_basis, exact_generators, center_dim
    )
    target = center_receipt["isotypic_dimensions_over_C"]
    types = [
        solution
        for solution in types
        if sorted(
            degree * multiplicity for degree, multiplicity in solution
        ) == target
    ]
    rational_factors = (
        rational_factors_from_center(types[0], center_receipt)
        if len(types) == 1 else None
    )
    return {
        "carrier_dimension": residual_dimension,
        "algebra_dimension": algebra_dimension,
        "center_dimension": center_dim,
        "commutant_dimension": commutant_dim,
        "modular_receipts": receipts,
        "complex_type_solution_count": len(types),
        "complex_type": [
            {"matrix_size": degree, "multiplicity": multiplicity}
            for degree, multiplicity in types[0]
        ] if len(types) == 1 else None,
        "generic_center_certificate": center_receipt,
        "rational_wedderburn_factors": rational_factors,
        "rational_classification_certified": rational_factors is not None,
        "rational_split_certified": (
            rational_factors is not None
            and all(factor["center_degree"] == 1 for factor in rational_factors)
        ),
        "exact_certificate": certificate,
    }


def analyze_group(row: dict, extra_vector=None) -> dict:
    if row["label"] == "Fi23" and extra_vector is None:
        raise ValueError("Fi23 needs the independently verified extra invariant vector; run restore_census.py")
    labels, history = stable_diagonal_partition(
        row["class_sizes"], row["power_maps"]
    )
    cells = partition_cells(labels)
    principal_size = len(cells) - (1 if extra_vector is not None else 0)
    if not quotient_connected(labels, row["power_maps"]):
        raise ArithmeticError(f"{row['label']} quotient power graph disconnected")
    residual = residual_algebra(
        row["class_sizes"], row["power_maps"], cells, extra_vector
    )
    if residual["complex_type"] is None:
        raise ArithmeticError(f"{row['label']} residual complex type is ambiguous")
    full_type = [
        {"matrix_size": principal_size, "multiplicity": 1},
        *residual["complex_type"],
    ]
    if not residual["rational_classification_certified"]:
        raise ArithmeticError(
            f"{row['label']} has an unresolved rational Wedderburn factor"
        )
    rational_factors = [{
        "center_degree": 1,
        "center_defining_polynomial": "lambda",
        "center_field": "Q",
        "matrix_size_over_center": principal_size,
        "module_multiplicity": 1,
        "schur_index": 1,
    }, *residual["rational_wedderburn_factors"]]
    algebra_dimension = principal_size**2 + residual["algebra_dimension"]
    center_dim = 1 + residual["center_dimension"]
    commutant_dim = 1 + residual["commutant_dimension"]
    noncommutative_residual = [
        block["matrix_size"]
        for block in residual["complex_type"]
        if block["matrix_size"] > 1
    ]
    polymatroid = scalable_polymatroid_certificate(row)
    return {
        "label": row["label"],
        "table_name": row["table_name"],
        "group_order": row["order"],
        "class_count": row["class_count"],
        "group_exponent": row["exponent"],
        "directions": row["directions"],
        "direction_count": len(row["directions"]),
        "stable_diagonal_history": history,
        "principal_block_size": principal_size,
        "residual_carrier_dimension": residual["carrier_dimension"],
        "algebra_dimension": algebra_dimension,
        "center_dimension": center_dim,
        "commutant_dimension": commutant_dim,
        "multiplicity_free": center_dim == commutant_dim,
        "complex_block_multiplicity_type": full_type,
        "rational_wedderburn_block_sizes": [
            block["matrix_size"] for block in full_type
        ],
        "rational_wedderburn_factors": rational_factors,
        "rational_classification_certified": True,
        "rational_split_certified": all(
            factor["center_degree"] == 1 for factor in rational_factors
        ),
        "residual_algebra": residual,
        "residual_noncommutative_block_sizes": noncommutative_residual,
        "residual_noncommutative_block_count": len(noncommutative_residual),
        "residual_noncommutative_algebra_dimension": sum(
            size * size for size in noncommutative_residual
        ),
        "residual_noncommutative_carrier_dimension": sum(
            noncommutative_residual
        ),
        "coherence_laplacian": laplacian_certificate(
            row["class_sizes"], row["power_maps"]
        ),
        "coherence_polymatroid": polymatroid,
    }


def prime_base(value: int) -> int:
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return divisor
        divisor += 1
    return value


def scalable_polymatroid_certificate(row: dict) -> dict:
    """Exact saturation data without an exponential submodularity audit."""

    size = row["class_count"]
    maps = row["power_maps"]
    directions = row["directions"]
    singleton_ranks = [
        size - component_count(size, [mapping]) for mapping in maps
    ]
    full_rank = size - component_count(size, maps)
    by_prime: dict[int, list[int]] = defaultdict(list)
    for index, direction in enumerate(directions):
        by_prime[prime_base(direction)].append(index)
    minimum_sets = []
    for selection in product(*by_prime.values()):
        chosen_maps = [maps[index] for index in selection]
        if size - component_count(size, chosen_maps) == full_rank:
            minimum_sets.append([directions[index] for index in selection])
    if not minimum_sets:
        raise ArithmeticError("no prime-covering set saturates coherence")
    return {
        "full_rank": full_rank,
        "nullity": size - full_rank,
        "singleton_ranks": singleton_ranks,
        "minimum_saturating_cardinality": len(by_prime),
        "minimum_saturating_subset_count": len(minimum_sets),
        "minimum_saturating_direction_sets": minimum_sets,
        "submodularity_exact_from_linear_representation": True,
        "is_matroid_rank": all(rank <= 1 for rank in singleton_ranks),
    }


def main() -> None:
    """Reconstructed CLI; computational functions above recovered from logs."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--gap", type=Path, default=DEFAULT_GAP)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reuse-gap-data", action="store_true")
    parser.add_argument("--groups", nargs="*")
    args = parser.parse_args()
    if not args.reuse_gap_data or not args.data.exists():
        export_gap_data(args.gap, args.data)
    source = json.loads(args.data.read_text())
    rows = []
    for row in source["groups"]:
        if args.groups and row["label"] not in args.groups:
            continue
        print(f"analyzing {row['label']} ({row['class_count']} classes)", flush=True)
        rows.append(analyze_group(row))
        result = {"schema": "sporadic-adams-star-restored-v1",
                  "provenance": provenance(args.data), "groups": rows,
                  "complete": len(rows) == (len(args.groups) if args.groups else len(source["groups"]))}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Wrote {len(rows)} groups to {args.output}")

if __name__ == "__main__":
    from restore_census import main as restored_main
    restored_main()
