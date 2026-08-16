#!/usr/bin/env python3
"""Extract S^3 seeds from a Lutz census and certify them.

The source census is not copied into the release.  This script reads the two
Lutz files in ``seeds/10/``, filters the S^3 entries, and applies the strict
Pachner-move predicates from the public verifier.  It writes a compact seed
list and one no-insertion simplification certificate per seed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Iterable


HERE = Path(__file__).resolve().parent
CERTIFICATES = HERE.parents[1]
sys.path.insert(0, str(CERTIFICATES))

from verify_pachner_certificate import (  # noqa: E402
    Facet,
    PachnerMove,
    apply_move,
    canonical_complex_sha256,
    certificate_document,
    facets_jsonable,
    move_from_alpha_beta,
    move_legality,
    normalize_facets,
    summarize_complex,
    vertices,
    verify_certificate,
)


ID_LINE = re.compile(r"^manifold_3_\d+_(\d+)=")
S3_ID = re.compile(r"manifold_3_\d+_(\d+): S\^3")
TETRAHEDRON = re.compile(r"\[\s*([0-9,\s]+?)\s*\]")


def enumerate_moves(facets: set[Facet], removed_count: int) -> list[PachnerMove]:
    """Enumerate legal 2-3, 3-2, or 4-1 moves by the strict criterion."""
    if removed_count not in (2, 3, 4):
        raise ValueError("only 2-3, 3-2, and 4-1 moves are supported")
    index: dict[Facet, set[Facet]] = defaultdict(set)
    alpha_size = 5 - removed_count
    for facet in facets:
        for alpha in combinations(sorted(facet), alpha_size):
            index[frozenset(alpha)].add(facet)
    moves: list[PachnerMove] = []
    for alpha, incident in index.items():
        if len(incident) != removed_count:
            continue
        support = set().union(*incident)
        if len(support) != 5:
            continue
        beta = frozenset(support) - alpha
        if len(beta) != removed_count:
            continue
        move = move_from_alpha_beta(3, alpha, beta)
        if set(move.removed) != incident:
            continue
        if move_legality(facets, move)[0]:
            moves.append(move)
    return sorted(
        moves,
        key=lambda move: (
            move.move_type,
            tuple(sorted(move.alpha)),
            tuple(sorted(move.beta)),
        ),
    )


def simplex_boundary(facets: set[Facet]) -> bool:
    vertex_set = vertices(facets)
    return len(vertex_set) == 5 and facets == {
        frozenset(facet) for facet in combinations(sorted(vertex_set), 4)
    }


def has_legal_3_2_move(facets: set[Facet]) -> bool:
    """Fast seed test equivalent to the strict 3-2 move predicate."""
    edge_to_facets: dict[tuple[int, int], list[Facet]] = defaultdict(list)
    triangular_faces: set[tuple[int, int, int]] = set()
    for facet in facets:
        for edge in combinations(sorted(facet), 2):
            edge_to_facets[edge].append(facet)
        triangular_faces.update(combinations(sorted(facet), 3))
    for edge, incident in edge_to_facets.items():
        if len(incident) != 3:
            continue
        support = set().union(*incident)
        if len(support) != 5:
            continue
        opposite_triangle = tuple(sorted(support - set(edge)))
        if opposite_triangle not in triangular_faces:
            return True
    return False


def monotone_path(facets: set[Facet]) -> list[PachnerMove] | None:
    """Find a 4-1/3-2-only path; every recursive step lowers f_3."""
    failed: set[str] = set()

    def visit(current: set[Facet]) -> list[PachnerMove] | None:
        if simplex_boundary(current):
            return []
        key = canonical_complex_sha256(current)
        if key in failed:
            return None
        for move in enumerate_moves(current, 4) + enumerate_moves(current, 3):
            suffix = visit(apply_move(current, move))
            if suffix is not None:
                return [move, *suffix]
        failed.add(key)
        return None

    return visit(set(facets))


def no_insertion_path(
    facets: set[Facet], max_2_3_moves: int = 3
) -> list[PachnerMove] | None:
    """Search using reductions and at most ``max_2_3_moves`` 2-3 bridges.

    The budget makes the search finite while allowing a temporary increase in
    tetrahedra when a seed has no immediately available reduction.  No 1-4
    move is ever considered.
    """
    for budget in range(max_2_3_moves + 1):
        failed: set[tuple[str, int]] = set()

        def visit(current: set[Facet], remaining_up_moves: int) -> list[PachnerMove] | None:
            if simplex_boundary(current):
                return []
            key = (canonical_complex_sha256(current), remaining_up_moves)
            if key in failed:
                return None
            for move in enumerate_moves(current, 4) + enumerate_moves(current, 3):
                suffix = visit(apply_move(current, move), remaining_up_moves)
                if suffix is not None:
                    return [move, *suffix]
            if remaining_up_moves:
                for move in enumerate_moves(current, 2):
                    suffix = visit(apply_move(current, move), remaining_up_moves - 1)
                    if suffix is not None:
                        return [move, *suffix]
            failed.add(key)
            return None

        result = visit(set(facets), budget)
        if result is not None:
            return result
    return None


def s3_ids(type_file: Path) -> set[int]:
    return {int(match) for match in S3_ID.findall(type_file.read_text())}


def lutz_records(census_file: Path, allowed_ids: set[int]) -> Iterable[tuple[int, set[Facet]]]:
    block: list[str] = []

    def parse(block_lines: list[str]) -> tuple[int, set[Facet]] | None:
        raw = "".join(block_lines)
        match = ID_LINE.match(raw)
        if match is None:
            raise ValueError(f"unrecognised census record beginning {raw[:80]!r}")
        identifier = int(match.group(1))
        if identifier not in allowed_ids:
            return None
        tetrahedra = [
            tuple(map(int, text.split(",")))
            for text in TETRAHEDRON.findall(raw)
        ]
        return identifier, normalize_facets(tetrahedra, 3)

    with census_file.open() as source:
        for line in source:
            if ID_LINE.match(line) and block:
                record = parse(block)
                if record is not None:
                    yield record
                block = [line]
                continue
            if line.strip():
                block.append(line)
                continue
            if not block:
                continue
            record = parse(block)
            block = []
            if record is not None:
                yield record
    if block:
        record = parse(block)
        if record is not None:
            yield record


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-directory", type=Path, default=HERE.parents[2] / "seeds/10")
    ap.add_argument("--output-directory", type=Path, default=HERE)
    ap.add_argument("--census-file", default="3_manifolds_10_all.txt")
    ap.add_argument("--type-file", default="3_manifolds_10_s3.txt")
    ap.add_argument("--vertex-count", type=int, default=10)
    ap.add_argument("--max-2-3-moves", type=int, default=3)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    source_dir = args.source_directory
    census = source_dir / args.census_file
    type_file = source_dir / args.type_file
    sphere_ids = s3_ids(type_file)
    seeds: list[tuple[int, set[Facet]]] = []
    for identifier, facets in lutz_records(census, sphere_ids):
        if not has_legal_3_2_move(facets):
            seeds.append((identifier, facets))

    seed_rows = []
    certificates: dict[str, dict] = {}
    for identifier, facets in seeds:
        path = no_insertion_path(facets, args.max_2_3_moves)
        if path is None:
            raise RuntimeError(f"no no-insertion simplification found for seed {identifier}")
        move_rows = [move.to_json(step=i) for i, move in enumerate(path)]
        final = set(facets)
        for move in path:
            final = apply_move(final, move)
        if not simplex_boundary(final):
            raise AssertionError(f"seed {identifier} did not reach a simplex boundary")
        cert = certificate_document(
            3,
            facets,
            move_rows,
            final,
            {
                "source_id": f"manifold_3_10_{identifier}",
                "purpose": f"no-insertion simplification of a {args.vertex_count}-vertex S^3 seed",
                "one_four_moves": 0,
                "seed_definition": "no legal 3-2 move",
            },
        )
        verify_certificate(facets, cert, paranoid=True)
        seed_rows.append(
            {
                "id": identifier,
                "source_id": f"manifold_3_10_{identifier}",
                "facets": facets_jsonable(facets),
                "summary": summarize_complex(facets),
                "legal_2_3_moves": len(enumerate_moves(facets, 2)),
                "legal_3_2_moves": 0,
            }
        )
        certificates[str(identifier)] = cert
        if not args.quiet:
            print(f"seed {identifier}: {len(path)} moves", flush=True)

    print(f"identified {len(seeds)} {args.vertex_count}-vertex S^3 seeds", flush=True)
    if args.dry_run:
        return
    args.output_directory.mkdir(parents=True, exist_ok=True)
    (args.output_directory / "seed-list.json").write_text(
        json.dumps(
            {
                "format": f"lutz-{args.vertex_count}-vertex-s3-seed-list-v1",
                "source": {
                    "census_file": census.name,
                    "census_sha256": sha256(census),
                    "s3_type_file": type_file.name,
                    "s3_type_sha256": sha256(type_file),
                    "s3_entries": len(sphere_ids),
                },
                "seed_definition": "no legal 3-2 bistellar move",
                "seeds": seed_rows,
            },
            indent=2,
        )
        + "\n"
    )
    (args.output_directory / "simplification-certificates.json").write_text(
        json.dumps(
            {
                "format": "n10-s3-seed-certificates-v1",
                "certificates": certificates,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
