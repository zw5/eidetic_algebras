"""Compute the finite-group Adams-star invariant package.

For every group G, E(G) consists of all prime powers dividing exp(G).
GAP supplies exact conjugacy-class power maps. Algebra closure, center,
and commutant dimensions are certified by exact rational FLINT ranks;
two good finite fields provide discovery and consistency receipts.
Graph and polymatroid statements are exact combinatorial calculations.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
from itertools import combinations
import json
from math import isqrt
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Iterable

import flint
import numpy as np
from flint import fmpq, fmpq_mat


HERE = Path(__file__).resolve().parent
REPOSITORY_ROOT = HERE.parents[2]
DEFAULT_GAP = Path(shutil.which("gap") or "gap")
DEFAULT_DATA = HERE / "output" / "gap_groups.json"
DEFAULT_OUTPUT = HERE / "output" / "finite_group_adams_star_census.json"
MODULAR_PRIMES = (1_000_003, 1_000_033)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def provenance(gap_data: Path) -> dict[str, object]:
    source = json.loads(gap_data.read_text())
    return {
        "python": sys.version.split()[0],
        "python_flint": flint.__version__,
        "flint": flint.__FLINT_VERSION__,
        "numpy": np.__version__,
        "gap": source.get("software", {}).get("gap"),
        "ctbllib": source.get("software", {}).get("ctbllib"),
        "raw_input_sha256": sha256(gap_data),
        "runner_sha256": sha256(Path(__file__)),
        "gap_export_sha256": sha256(HERE / "export_gap_groups.g"),
    }


class ModularBasis:
    def __init__(self, width: int, prime: int):
        self.width = width
        self.prime = prime
        self.rows: dict[int, np.ndarray] = {}

    def add(self, values) -> bool:
        row = np.asarray(values, dtype=np.int64).reshape(-1) % self.prime
        if len(row) != self.width:
            raise ValueError("wrong row width")
        for column, pivot in sorted(self.rows.items()):
            if row[column]:
                row = (row - row[column] * pivot) % self.prime
        nonzero = np.flatnonzero(row)
        if not len(nonzero):
            return False
        column = int(nonzero[0])
        row = row * pow(int(row[column]), self.prime - 2, self.prime) % self.prime
        for old_column, old_row in tuple(self.rows.items()):
            if old_row[column]:
                self.rows[old_column] = (
                    old_row - old_row[column] * row
                ) % self.prime
        self.rows[column] = row
        return True


class DSU:
    def __init__(self, size: int):
        self.parent = list(range(size))
        self.count = size

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[right] = left
            self.count -= 1


def export_gap_data(gap: Path, script: Path, output: Path) -> None:
    if not gap.exists():
        raise FileNotFoundError(f"GAP executable not found: {gap}")
    completed = subprocess.run(
        [str(gap), "-q", str(script)],
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(completed.stdout)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n")


def pullback_matrix(mapping: list[int], prime: int) -> np.ndarray:
    size = len(mapping)
    matrix = np.zeros((size, size), dtype=np.int64)
    for source, target_one_based in enumerate(mapping):
        matrix[source, target_one_based - 1] = 1
    return matrix % prime


def adjoint_matrix(
    pullback: np.ndarray,
    class_sizes: list[int],
    prime: int,
) -> np.ndarray:
    sizes = np.asarray([value % prime for value in class_sizes], dtype=np.int64)
    if np.any(sizes == 0):
        raise ValueError("bad reduction prime divides a class size")
    inverses = np.asarray(
        [pow(int(value), prime - 2, prime) for value in sizes],
        dtype=np.int64,
    )
    # D^{-1} P^T D in class coordinates.
    return (
        inverses[:, None] * pullback.T % prime
    ) * sizes[None, :] % prime


def star_closure(generators: tuple[np.ndarray, ...], prime: int):
    size = generators[0].shape[0]
    basis = ModularBasis(size * size, prime)
    matrices = []
    frontier = []
    for matrix in (np.eye(size, dtype=np.int64),) + generators:
        matrix = matrix % prime
        if basis.add(matrix):
            matrices.append(matrix)
            frontier.append(matrix)
    while frontier:
        new = []
        for left in frontier:
            for right in generators:
                product = left @ right % prime
                if basis.add(product):
                    matrices.append(product)
                    new.append(product)
        frontier = new
    return tuple(matrices)


def exact_pullback_matrix(mapping: list[int]) -> fmpq_mat:
    size = len(mapping)
    rows = [[0] * size for _ in range(size)]
    for source, target_one_based in enumerate(mapping):
        rows[source][target_one_based - 1] = 1
    return fmpq_mat(rows)


def exact_adjoint_matrix(
    pullback: fmpq_mat,
    class_sizes: list[int],
) -> fmpq_mat:
    size = pullback.nrows()
    return fmpq_mat(
        [
            [
                fmpq(int(pullback[column, row]) * class_sizes[column],
                     class_sizes[row])
                for column in range(size)
            ]
            for row in range(size)
        ]
    )


def exact_star_closure_certificate(
    maps: list[list[int]],
    class_sizes: list[int],
    expected_dimension: int,
    expected_center_dimension: int,
    expected_commutant_dimension: int,
    discovery_prime: int = MODULAR_PRIMES[0],
) -> dict:
    """Verify the discovered star closure over Q using FLINT.

    The modular basis chooses a compact word basis.  Its exact rational
    lifts must have the same rank, and the rank must remain unchanged after
    adjoining every basis-by-generator product.
    """

    exact_pullbacks = tuple(exact_pullback_matrix(mapping) for mapping in maps)
    exact_adjoints = tuple(
        exact_adjoint_matrix(matrix, class_sizes)
        for matrix in exact_pullbacks
    )
    exact_generators = exact_pullbacks + exact_adjoints
    size = exact_generators[0].nrows()
    modular_generators = tuple(
        pullback_matrix(mapping, discovery_prime) for mapping in maps
    ) + tuple(
        adjoint_matrix(
            pullback_matrix(mapping, discovery_prime),
            class_sizes,
            discovery_prime,
        )
        for mapping in maps
    )

    modular_basis = ModularBasis(size * size, discovery_prime)
    exact_basis: list[fmpq_mat] = []
    frontier: list[tuple[np.ndarray, fmpq_mat]] = []
    all_closure_rows: list[list[fmpq]] = []
    exact_identity = fmpq_mat(
        [[1 if row == column else 0 for column in range(size)]
         for row in range(size)]
    )
    modular_identity = np.eye(size, dtype=np.int64)
    for modular, exact in (
        (modular_identity, exact_identity),
        *zip(modular_generators, exact_generators, strict=True),
    ):
        all_closure_rows.append(list(exact.entries()))
        if modular_basis.add(modular):
            exact_basis.append(exact)
            frontier.append((modular % discovery_prime, exact))

    while frontier:
        new = []
        for modular_left, exact_left in frontier:
            for modular_right, exact_right in zip(
                modular_generators, exact_generators, strict=True
            ):
                modular_product = (
                    modular_left @ modular_right % discovery_prime
                )
                exact_product = exact_left * exact_right
                all_closure_rows.append(list(exact_product.entries()))
                if modular_basis.add(modular_product):
                    exact_basis.append(exact_product)
                    new.append((modular_product, exact_product))
        frontier = new

    basis_rank = fmpq_mat(
        [list(matrix.entries()) for matrix in exact_basis]
    ).rank()
    closure_rank = fmpq_mat(all_closure_rows).rank()

    center_constraints = []
    for generator in exact_generators:
        commutators = [
            element * generator - generator * element
            for element in exact_basis
        ]
        for row in range(size):
            for column in range(size):
                center_constraints.append(
                    [commutator[row, column] for commutator in commutators]
                )
    center_constraint_rank = fmpq_mat(center_constraints).rank()
    exact_center_dimension = len(exact_basis) - center_constraint_rank

    commutant_constraints = []
    for generator in exact_generators:
        for row in range(size):
            for column in range(size):
                equation = [fmpq(0) for _ in range(size * size)]
                for middle in range(size):
                    equation[row * size + middle] += generator[middle, column]
                    equation[middle * size + column] -= generator[row, middle]
                commutant_constraints.append(equation)
    commutant_constraint_rank = fmpq_mat(commutant_constraints).rank()
    exact_commutant_dimension = size * size - commutant_constraint_rank
    return {
        "discovery_prime": discovery_prime,
        "word_basis_size": len(exact_basis),
        "exact_word_basis_rank": basis_rank,
        "exact_basis_plus_generator_products_rank": closure_rank,
        "expected_dimension": expected_dimension,
        "exact_center_dimension": exact_center_dimension,
        "expected_center_dimension": expected_center_dimension,
        "exact_commutant_dimension": exact_commutant_dimension,
        "expected_commutant_dimension": expected_commutant_dimension,
        "exact_over_Q": (
            len(exact_basis)
            == basis_rank
            == closure_rank
            == expected_dimension
            and exact_center_dimension == expected_center_dimension
            and exact_commutant_dimension == expected_commutant_dimension
        ),
    }


def commutant_dimension(generators, prime: int) -> int:
    size = generators[0].shape[0]
    constraints = ModularBasis(size * size, prime)
    for matrix in generators:
        nonzero_columns = [
            np.flatnonzero(matrix[:, column]) for column in range(size)
        ]
        nonzero_rows = [
            np.flatnonzero(matrix[row, :]) for row in range(size)
        ]
        for row in range(size):
            for column in range(size):
                equation = np.zeros(size * size, dtype=np.int64)
                for middle in nonzero_columns[column]:
                    equation[row * size + middle] += matrix[middle, column]
                for middle in nonzero_rows[row]:
                    equation[middle * size + column] -= matrix[row, middle]
                constraints.add(equation)
    return size * size - len(constraints.rows)


def center_dimension(algebra, generators, prime: int) -> int:
    size = generators[0].shape[0]
    constraints = ModularBasis(len(algebra), prime)
    for generator in generators:
        commutators = [
            (element @ generator - generator @ element) % prime
            for element in algebra
        ]
        for row in range(size):
            for column in range(size):
                constraints.add(
                    [commutator[row, column] for commutator in commutators]
                )
    return len(algebra) - len(constraints.rows)


def component_count(size: int, maps: Iterable[list[int]]) -> int:
    dsu = DSU(size)
    for mapping in maps:
        for source, target_one_based in enumerate(mapping):
            dsu.union(source, target_one_based - 1)
    return dsu.count


def subset_rank(size: int, maps: list[list[int]], mask: int) -> int:
    selected = [
        mapping for index, mapping in enumerate(maps) if mask & (1 << index)
    ]
    return size - component_count(size, selected)


def polymatroid_certificate(size: int, maps: list[list[int]]) -> dict:
    direction_count = len(maps)
    ranks = {
        mask: subset_rank(size, maps, mask)
        for mask in range(1 << direction_count)
    }
    full_rank = ranks[(1 << direction_count) - 1]
    minimum = min(
        mask.bit_count() for mask, rank in ranks.items() if rank == full_rank
    )
    minimum_count = sum(
        rank == full_rank and mask.bit_count() == minimum
        for mask, rank in ranks.items()
    )
    minimum_masks = [
        mask
        for mask, rank in ranks.items()
        if rank == full_rank and mask.bit_count() == minimum
    ]
    submodular = all(
        ranks[left] + ranks[right]
        >= ranks[left | right] + ranks[left & right]
        for left in ranks
        for right in ranks
    )
    return {
        "full_rank": full_rank,
        "nullity": size - full_rank,
        "singleton_ranks": [
            ranks[1 << index] for index in range(direction_count)
        ],
        "minimum_saturating_cardinality": minimum,
        "minimum_saturating_subset_count": minimum_count,
        "minimum_saturating_masks": minimum_masks,
        "submodularity_exhaustively_checked": submodular,
        "is_matroid_rank": all(
            ranks[1 << index] <= 1 for index in range(direction_count)
        ),
        "rank_distribution": {
            str(rank): sum(value == rank for value in ranks.values())
            for rank in sorted(set(ranks.values()))
        },
    }


def laplacian_certificate(
    class_sizes: list[int],
    maps: list[list[int]],
) -> dict:
    size = len(class_sizes)
    weights = np.asarray(class_sizes, dtype=np.float64)
    laplacian = np.zeros((size, size), dtype=np.float64)
    identity = np.eye(size)
    for mapping in maps:
        matrix = np.zeros((size, size), dtype=np.float64)
        for source, target_one_based in enumerate(mapping):
            target = target_one_based - 1
            matrix[source, target] = np.sqrt(weights[source] / weights[target])
        displacement = identity - matrix
        laplacian += displacement.T @ displacement
    eigenvalues = np.linalg.eigvalsh((laplacian + laplacian.T) / 2)
    exact_nullity = component_count(size, maps)
    positive = eigenvalues[exact_nullity:]
    return {
        "exact_nullity": exact_nullity,
        "exact_rank": size - exact_nullity,
        "spectral_gap_float": float(positive[0]) if len(positive) else None,
        "largest_eigenvalue_float": float(eigenvalues[-1]),
        "trace_float": float(np.trace(laplacian)),
    }


def complex_type_solutions(
    dimension: int,
    algebra_dimension: int,
    commutant_dimension_value: int,
    center_dimension_value: int,
    limit: int = 200,
):
    """Infer possible complex block sizes d_j and multiplicities m_j.

    For A_C = sum M_{d_j}(C) acting on sum (C^{d_j})^{m_j}:
      dim A = sum d_j^2,
      dim A' = sum m_j^2,
      dim V = sum d_j m_j,
      dim Z(A_C) = number of blocks.
    """

    pairs = [
        (degree, multiplicity)
        for degree in range(1, dimension + 1)
        for multiplicity in range(1, dimension + 1)
        if degree * multiplicity <= dimension
        and degree * degree <= algebra_dimension
        and multiplicity * multiplicity <= commutant_dimension_value
    ]
    solutions = []

    def search(
        slot: int,
        start: int,
        remaining_dimension: int,
        remaining_algebra: int,
        remaining_commutant: int,
        chosen: list[tuple[int, int]],
    ):
        if len(solutions) >= limit:
            return
        remaining_slots = center_dimension_value - slot
        if remaining_slots == 0:
            if (
                remaining_dimension == 0
                and remaining_algebra == 0
                and remaining_commutant == 0
            ):
                solutions.append(tuple(chosen))
            return
        if min(remaining_dimension, remaining_algebra, remaining_commutant) < remaining_slots:
            return
        for index in range(start, len(pairs)):
            degree, multiplicity = pairs[index]
            product_dimension = degree * multiplicity
            if (
                product_dimension > remaining_dimension
                or degree * degree > remaining_algebra
                or multiplicity * multiplicity > remaining_commutant
            ):
                continue
            search(
                slot + 1,
                index,
                remaining_dimension - product_dimension,
                remaining_algebra - degree * degree,
                remaining_commutant - multiplicity * multiplicity,
                chosen + [(degree, multiplicity)],
            )

    search(
        0,
        0,
        dimension,
        algebra_dimension,
        commutant_dimension_value,
        [],
    )
    return solutions


def analyze_group(row: dict) -> dict:
    size = row["class_count"]
    maps = row["power_maps"]
    prime_receipts = []
    reference = None
    for prime in MODULAR_PRIMES:
        pullbacks = tuple(pullback_matrix(mapping, prime) for mapping in maps)
        adjoints = tuple(
            adjoint_matrix(matrix, row["class_sizes"], prime)
            for matrix in pullbacks
        )
        generators = pullbacks + adjoints
        algebra = star_closure(generators, prime)
        receipt = {
            "prime": prime,
            "algebra_dimension": len(algebra),
            "center_dimension": center_dimension(algebra, generators, prime),
            "commutant_dimension": commutant_dimension(generators, prime),
        }
        prime_receipts.append(receipt)
        key = (
            receipt["algebra_dimension"],
            receipt["center_dimension"],
            receipt["commutant_dimension"],
        )
        if reference is None:
            reference = key
        elif key != reference:
            raise ArithmeticError(f"modular receipts disagree for {row['label']}")

    algebra_dimension, center_dim, commutant_dim = reference
    rational_closure = exact_star_closure_certificate(
        maps,
        row["class_sizes"],
        algebra_dimension,
        center_dim,
        commutant_dim,
    )
    if not rational_closure["exact_over_Q"]:
        raise ArithmeticError(
            f"rational closure failed for {row['label']}: {rational_closure}"
        )
    type_solutions = complex_type_solutions(
        size, algebra_dimension, commutant_dim, center_dim
    )
    type_payload = [
        [
            {"matrix_size": degree, "multiplicity": multiplicity}
            for degree, multiplicity in solution
        ]
        for solution in type_solutions[:10]
    ]
    polymatroid = polymatroid_certificate(size, maps)
    polymatroid["minimum_saturating_direction_sets"] = [
        [
            row["directions"][index]
            for index in range(len(row["directions"]))
            if mask & (1 << index)
        ]
        for mask in polymatroid.pop("minimum_saturating_masks")
    ]
    return {
        "label": row["label"],
        "family": row["family"],
        "group_order": row["order"],
        "class_count": size,
        "group_exponent": row["exponent"],
        "directions": row["directions"],
        "direction_count": len(row["directions"]),
        "algebra_dimension": algebra_dimension,
        "ambient_endomorphism_dimension": size * size,
        "algebra_density": algebra_dimension / (size * size),
        "center_dimension": center_dim,
        "commutant_dimension": commutant_dim,
        "double_commutant_dimension": algebra_dimension,
        "multiplicity_free": commutant_dim == center_dim,
        "modular_receipts": prime_receipts,
        "rational_closure_certificate": rational_closure,
        "complex_type_solution_count_capped": len(type_solutions),
        "complex_block_multiplicity_type": (
            type_payload[0] if len(type_solutions) == 1 else None
        ),
        "first_complex_type_possibilities": type_payload,
        "coherence_laplacian": laplacian_certificate(
            row["class_sizes"], maps
        ),
        "coherence_polymatroid": polymatroid,
    }


def monster_flagship() -> dict:
    path = REPOSITORY_ROOT / "data" / "census.json"
    rows = json.loads(path.read_text())["groups"]
    data = next(row for row in rows if row["label"] == "Monster")
    block_sizes = sorted(
        data["rational_wedderburn_block_sizes"], reverse=True
    )
    polymatroid = data["coherence_polymatroid"]
    return {
        "label": "Monster",
        "family": "sporadic_simple_flagship",
        "group_order": str(data["group_order"]),
        "class_count": data["class_count"],
        "direction_count": data["direction_count"],
        "directions": data["directions"],
        "algebra_dimension": data["algebra_dimension"],
        "ambient_endomorphism_dimension": data["class_count"] ** 2,
        "algebra_density": (
            data["algebra_dimension"] / data["class_count"] ** 2
        ),
        "center_dimension": data["center_dimension"],
        "commutant_dimension": data["commutant_dimension"],
        "double_commutant_dimension": data["algebra_dimension"],
        "multiplicity_free": data["multiplicity_free"],
        "rational_wedderburn_block_sizes": block_sizes,
        "coherence_laplacian": {
            "exact_nullity": data["coherence_laplacian"]["exact_nullity"],
            "exact_rank": data["coherence_laplacian"]["exact_rank"],
        },
        "coherence_polymatroid": {
            "full_rank": polymatroid["full_rank"],
            "minimum_saturating_cardinality": (
                polymatroid["minimum_saturating_cardinality"]
            ),
            "minimum_saturating_subset_count": (
                polymatroid["minimum_saturating_subset_count"]
            ),
            "submodularity_exact_from_linear_representation": (
                polymatroid[
                    "submodularity_exact_from_linear_representation"
                ]
            ),
            "is_matroid_rank": polymatroid["is_matroid_rank"],
        },
    }


def summarize(rows: list[dict]) -> dict:
    finite = [row for row in rows if row["label"] != "Monster"]
    return {
        "finite_group_count": len(finite),
        "family_count": len({row["family"] for row in finite}),
        "multiplicity_free_groups": [
            row["label"] for row in finite if row["multiplicity_free"]
        ],
        "full_matrix_algebra_groups": [
            row["label"]
            for row in finite
            if row["algebra_dimension"]
            == row["ambient_endomorphism_dimension"]
        ],
        "connected_coherence_groups": [
            row["label"]
            for row in finite
            if row["coherence_laplacian"]["exact_nullity"] == 1
        ],
        "matroid_rank_groups": [
            row["label"]
            for row in finite
            if row["coherence_polymatroid"]["is_matroid_rank"]
        ],
        "largest_nonmonster_algebra": max(
            finite, key=lambda row: row["algebra_dimension"]
        )["label"],
        "highest_nonmonster_density": max(
            finite, key=lambda row: row["algebra_density"]
        )["label"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gap", type=Path, default=DEFAULT_GAP)
    parser.add_argument("--gap-data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reuse-gap-data", action="store_true")
    args = parser.parse_args()
    if not args.reuse_gap_data or not args.gap_data.exists():
        export_gap_data(args.gap, HERE / "export_gap_groups.g", args.gap_data)
    source = json.loads(args.gap_data.read_text())
    rows = [analyze_group(row) for row in source["groups"]]
    rows.append(monster_flagship())
    result = {
        "schema": "finite-group-adams-star-census-v2",
        "provenance": provenance(args.gap_data),
        "exponent_set_rule": "all prime powers p^a dividing exp(G)",
        "modular_primes": MODULAR_PRIMES,
        "groups": rows,
        "summary": summarize(rows),
        "claim_boundary": (
            "For every non-Monster census row, the algebra dimension and closure, "
            "center dimension, and commutant dimension are certified over Q by "
            "exact FLINT rational ranks; the two good-prime computations are "
            "independent discovery and consistency receipts. Displayed block "
            "types describe the complexification and are reported only when the "
            "four exact dimension equations have a unique positive-integral "
            "solution. A rational division-algebra refinement requires separate "
            "central factorization; the Monster split-Q decomposition has its "
            "own exact certificate."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": "completed",
        "output": str(args.output),
        "group_count": len(rows),
        "summary": result["summary"],
    }, indent=2))


if __name__ == "__main__":
    main()
