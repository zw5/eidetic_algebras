from __future__ import annotations

import json
from math import gcd
from pathlib import Path

from runners.experiments.finite_group_adams_star_census_20260724.run_census import (
    DEFAULT_OUTPUT,
    analyze_group,
    complex_type_solutions,
    polymatroid_certificate,
)


def test_complex_type_solver_on_split_multiplicity_free_example() -> None:
    solutions = complex_type_solutions(
        dimension=5,
        algebra_dimension=13,
        commutant_dimension_value=2,
        center_dimension_value=2,
    )
    assert solutions == [((2, 1), (3, 1))]


def test_polymatroid_rank_for_cyclic_four_power_maps() -> None:
    maps = [
        [1, 3, 1, 3],  # squaring on C4 in a suitable element order
        [1, 1, 1, 1],  # fourth power
    ]
    certificate = polymatroid_certificate(4, maps)
    assert certificate["submodularity_exhaustively_checked"]
    assert certificate["full_rank"] == 3
    assert certificate["minimum_saturating_cardinality"] == 1


def test_c2_exact_adams_star_package() -> None:
    row = {
        "label": "C2",
        "family": "cyclic",
        "order": 2,
        "class_count": 2,
        "exponent": 2,
        "directions": [2],
        "class_sizes": [1, 1],
        "class_orders": [1, 2],
        "power_maps": [[1, 1]],
    }
    result = analyze_group(row)
    assert result["algebra_dimension"] == 4
    assert result["center_dimension"] == 1
    assert result["commutant_dimension"] == 1
    assert result["multiplicity_free"]
    assert result["coherence_laplacian"]["exact_rank"] == 1
    assert result["rational_closure_certificate"]["exact_over_Q"]


def _cyclic_prime_power_row(prime: int, height: int) -> dict[str, object]:
    order = prime**height
    directions = [prime**power for power in range(1, height + 1)]
    return {
        "label": f"C{order}",
        "family": "cyclic_prime_power",
        "order": order,
        "class_count": order,
        "exponent": order,
        "directions": directions,
        "class_sizes": [1] * order,
        "class_orders": [
            1 if value == 0 else order // gcd(order, value)
            for value in range(order)
        ],
        "power_maps": [
            [((direction * value) % order) + 1 for value in range(order)]
            for direction in directions
        ],
    }


def test_cyclic_prime_power_family_formula_in_new_exact_cases() -> None:
    c9 = analyze_group(_cyclic_prime_power_row(3, 2))
    assert c9["algebra_dimension"] == 14
    assert c9["complex_block_multiplicity_type"] == [
        {"matrix_size": 1, "multiplicity": 4},
        {"matrix_size": 2, "multiplicity": 1},
        {"matrix_size": 3, "multiplicity": 1},
    ]
    assert c9["commutant_dimension"] == 18
    assert c9["rational_closure_certificate"]["exact_over_Q"]

    c16 = analyze_group(_cyclic_prime_power_row(2, 4))
    assert c16["algebra_dimension"] == 39
    assert c16["complex_block_multiplicity_type"] == [
        {"matrix_size": 1, "multiplicity": 4},
        {"matrix_size": 2, "multiplicity": 2},
        {"matrix_size": 3, "multiplicity": 1},
        {"matrix_size": 5, "multiplicity": 1},
    ]
    assert c16["commutant_dimension"] == 22
    assert c16["rational_closure_certificate"]["exact_over_Q"]


def test_saved_census_has_exact_rational_receipts_and_unique_complex_types() -> None:
    assert Path(DEFAULT_OUTPUT).exists()
    census = json.loads(Path(DEFAULT_OUTPUT).read_text())
    finite_rows = [
        row for row in census["groups"] if "rational_closure_certificate" in row
    ]
    assert len(finite_rows) == 31
    assert all(
        row["rational_closure_certificate"]["exact_over_Q"] for row in finite_rows
    )
    assert all(row["complex_type_solution_count_capped"] == 1 for row in finite_rows)


def test_saved_census_records_reproducible_provenance() -> None:
    census = json.loads(Path(DEFAULT_OUTPUT).read_text())
    assert census["schema"] == "finite-group-adams-star-census-v2"
    provenance = census["provenance"]
    assert provenance["gap"] == "4.14.0"
    assert provenance["ctbllib"] == "1.3.9"
    import flint
    assert provenance["flint"] == flint.__FLINT_VERSION__
    assert all(
        len(provenance[key]) == 64
        for key in (
            "raw_input_sha256",
            "runner_sha256",
            "gap_export_sha256",
        )
    )


def test_saved_census_distinguishes_symmetric_group_block_geometries() -> None:
    census = json.loads(Path(DEFAULT_OUTPUT).read_text())
    rows = {row["label"]: row for row in census["groups"]}
    assert rows["S3"]["complex_block_multiplicity_type"] == [
        {"matrix_size": 3, "multiplicity": 1}
    ]
    assert rows["S6"]["complex_block_multiplicity_type"] == [
        {"matrix_size": 1, "multiplicity": 1},
        {"matrix_size": 3, "multiplicity": 1},
        {"matrix_size": 7, "multiplicity": 1},
    ]
    assert rows["S8"]["complex_block_multiplicity_type"] == [
        {"matrix_size": 1, "multiplicity": 1},
        {"matrix_size": 21, "multiplicity": 1},
    ]
