"""Pure exact helpers recovered from the 2026-07-27 Monster certificate."""
from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterable
from math import prod
from sympy import factorint, primitive_root
from sympy.ntheory.modular import crt


def _canonical_colors(signatures: Iterable[object]) -> list[int]:
    labels: dict[object, int] = {}
    out: list[int] = []
    for signature in signatures:
        if signature not in labels:
            labels[signature] = len(labels)
        out.append(labels[signature])
    return out

def _rank_mod(rows: list[list[int]], prime: int) -> int:
    matrix = [[entry % prime for entry in row] for row in rows]
    if not matrix:
        return 0
    row = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (i for i in range(row, len(matrix)) if matrix[i][column]), None
        )
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        inverse = pow(matrix[row][column], -1, prime)
        matrix[row] = [(inverse * entry) % prime for entry in matrix[row]]
        for i in range(len(matrix)):
            if i == row or matrix[i][column] == 0:
                continue
            scale = matrix[i][column]
            matrix[i] = [
                (left - scale * right) % prime
                for left, right in zip(matrix[i], matrix[row])
            ]
        row += 1
        if row == len(matrix):
            break
    return row

def _refined_cells(
    registry: MonsterRegistry, primes: list[int], maps: list[list[int]]
) -> tuple[list[int], list[list[int]], int]:
    """Canonical color refinement for unary colors and directed functions."""
    count = len(registry.names)
    colors = _canonical_colors(
        (registry.orders[i], registry.class_sizes[i]) for i in range(count)
    )
    rounds = 0
    while True:
        signatures = []
        for i in range(count):
            signature: list[object] = [colors[i]]
            for power_map in maps:
                signature.append(colors[power_map[i]])
                incoming = Counter(
                    colors[j] for j in range(count) if power_map[j] == i
                )
                signature.append(tuple(sorted(incoming.items())))
            signatures.append(tuple(signature))
        new_colors = _canonical_colors(signatures)
        rounds += 1
        if new_colors == colors:
            break
        colors = new_colors
    cells_by_color: dict[int, list[int]] = defaultdict(list)
    for i, color in enumerate(colors):
        cells_by_color[color].append(i)
    cells = list(cells_by_color.values())
    return colors, cells, rounds

@dataclass(frozen=True)
class SwapSystem:
    pair_cells: tuple[tuple[int, int], ...]
    components: tuple[tuple[int, ...], ...]
    pinned_components: tuple[tuple[int, ...], ...]

def _swap_system(cells: list[list[int]], maps: list[list[int]]) -> SwapSystem:
    """Solve all power-map compatibility equations among two-point cells.

    Color refinement proves that every automorphism fixes each cell setwise.
    Every non-singleton cell is a pair, so a candidate automorphism is a
    bit-vector of pair swaps.  A power map either equates two swap bits or
    pins one bit to zero.  Union-find therefore solves the full group.
    """
    assert all(len(cell) in (1, 2) for cell in cells)
    pairs = tuple(tuple(cell) for cell in cells if len(cell) == 2)
    pair_of = {point: i for i, pair in enumerate(pairs) for point in pair}
    partner = {a: b for a, b in pairs} | {b: a for a, b in pairs}

    parent = list(range(len(pairs)))
    pinned: set[int] = set()

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for power_map in maps:
        for point in range(len(power_map)):
            source_pair = pair_of.get(point)
            target = power_map[point]
            target_pair = pair_of.get(target)
            if source_pair is None:
                if target_pair is not None:
                    pinned.add(target_pair)
                continue
            other_target = power_map[partner[point]]
            if target_pair is None:
                if other_target != target:
                    pinned.add(source_pair)
                continue
            if pair_of.get(other_target) != target_pair or other_target == target:
                pinned.add(source_pair)
                pinned.add(target_pair)
            else:
                union(source_pair, target_pair)

    components_by_root: dict[int, list[int]] = defaultdict(list)
    for i in range(len(pairs)):
        components_by_root[find(i)].append(i)
    pinned_roots = {find(i) for i in pinned}
    free = tuple(
        tuple(component)
        for root, component in components_by_root.items()
        if root not in pinned_roots
    )
    fixed = tuple(
        tuple(component)
        for root, component in components_by_root.items()
        if root in pinned_roots
    )

    # Exhaustively verify the generator of every free component.
    for component in free:
        permutation = list(range(sum(len(cell) for cell in cells)))
        for pair_index in component:
            left, right = pairs[pair_index]
            permutation[left], permutation[right] = right, left
        for power_map in maps:
            assert all(
                permutation[power_map[i]] == power_map[permutation[i]]
                for i in range(len(permutation))
            )

    return SwapSystem(pairs, free, fixed)

def _local_unit_generators(prime: int, exponent: int) -> list[tuple[int, int]]:
    modulus = prime**exponent
    if prime == 2:
        if exponent == 1:
            return []
        if exponent == 2:
            return [(modulus - 1, 2)]
        return [(modulus - 1, 2), (5, 2 ** (exponent - 2))]
    generator = int(primitive_root(modulus))
    return [(generator, (prime - 1) * prime ** (exponent - 1))]

def _global_unit_generators(modulus: int) -> list[tuple[int, int]]:
    factors = factorint(modulus)
    moduli = [prime**exponent for prime, exponent in factors.items()]
    generators: list[tuple[int, int]] = []
    for position, (prime, exponent) in enumerate(factors.items()):
        for local, order in _local_unit_generators(prime, exponent):
            residues = [1] * len(moduli)
            residues[position] = local
            lifted, _ = crt(moduli, residues)
            generators.append((int(lifted), order))
    assert prod(order for _, order in generators) == int(
        prod(
            prime ** (exponent - 1) * (prime - 1)
            for prime, exponent in factors.items()
        )
    )
    return generators
