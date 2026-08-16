#!/usr/bin/env python3
"""Build the public high-valence twelve-vertex seed artifacts.

``seeds/12/12seeds.txt`` has 134 high-valence closed 3-manifold seed facet
lists.  The companion homology calculation filters this raw list to the 37
entries with S^3-type homology over Q and six finite fields.  For each of
those 37 entries, this script records a no-1-4 Pachner path to the boundary of
a 4-simplex; the path itself is the final, independent S^3 certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERTIFICATES = HERE.parents[1]
sys.path.insert(0, str(CERTIFICATES))
sys.path.insert(0, str(HERE.parent / "n10"))

from build_n10_s3_seed_certificates import (  # noqa: E402
    enumerate_moves,
    no_insertion_path,
    simplex_boundary,
)
from verify_pachner_certificate import (  # noqa: E402
    apply_move,
    certificate_document,
    check_connected,
    check_orientable,
    facets_jsonable,
    normalize_facets,
    summarize_complex,
    verify_certificate,
    vertices,
)


RECORD = re.compile(r"^## (\d+), (.*?)(?=^## \d+,|\Z)", re.MULTILINE | re.DOTALL)
TETRAHEDRON = re.compile(r"\[([0-9,]+)\]")
DEGREES = re.compile(r"^##  deg = (.+)$", re.MULTILINE)
S3_BETTI = [1, 0, 0, 1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_records(source: Path):
    text = source.read_text()
    for match in RECORD.finditer(text):
        identifier = int(match.group(1))
        body = match.group(2)
        facets = normalize_facets(
            [tuple(map(int, item.split(","))) for item in TETRAHEDRON.findall(body)],
            3,
        )
        degree_match = DEGREES.search(body)
        degree_sequence = (
            [int(value) for value in degree_match.group(1).split(",")]
            if degree_match
            else None
        )
        yield identifier, facets, degree_sequence


def edge_valences(facets: set) -> dict[tuple[int, int], int]:
    counts: Counter[tuple[int, int]] = Counter()
    for facet in facets:
        counts.update(combinations(sorted(facet), 2))
    return dict(counts)


def vertex_links_are_spheres(facets: set) -> bool:
    """Check the local condition for a closed combinatorial 3-manifold."""
    for vertex in vertices(facets):
        link = {facet - {vertex} for facet in facets if vertex in facet}
        summary = summarize_complex(link)
        if not summary["closed_pseudomanifold"] or not summary["connected"]:
            return False
        if summary["euler_characteristic"] != 2 or not check_orientable(link):
            return False
    return True


def checked_seed_row(identifier: int, facets: set, degree_sequence: list[int] | None) -> dict:
    summary = summarize_complex(facets)
    if len(vertices(facets)) != 12:
        raise ValueError(f"record {identifier} does not have twelve vertices")
    if not summary["closed_pseudomanifold"]:
        raise ValueError(f"record {identifier} is not a closed pseudomanifold")
    if not summary["connected"] or not check_connected(facets):
        raise ValueError(f"record {identifier} is disconnected")
    if not vertex_links_are_spheres(facets):
        raise ValueError(f"record {identifier} is not a combinatorial 3-manifold")
    orientable = check_orientable(facets)
    valences = edge_valences(facets)
    if min(valences.values()) < 4:
        raise ValueError(f"record {identifier} has an edge of valence below four")
    if enumerate_moves(facets, 3):
        raise ValueError(f"record {identifier} has a legal 3-2 move")
    return {
        "id": identifier,
        "facets": facets_jsonable(facets),
        "summary": summary,
        "orientable": orientable,
        "tetrahedron_degrees": degree_sequence,
        "edge_valence_range": [min(valences.values()), max(valences.values())],
        "legal_3_2_moves": 0,
    }


def load_homology_candidates(path: Path, source_sha256: str) -> tuple[dict[int, dict], dict]:
    document = json.loads(path.read_text())
    if document.get("format") != "n12-homology-over-fields-v1":
        raise ValueError(f"{path} is not an n12 homology-over-fields document")
    if document["source"]["sha256"] != source_sha256:
        raise ValueError("homology result was not computed from the requested seed file")
    entries = {row["id"]: row for row in document["results"]}
    if len(entries) != 134:
        raise ValueError("homology result must contain 134 distinct records")
    candidates = {
        identifier: row
        for identifier, row in entries.items()
        if row["s3_homology_over_tested_fields"]
    }
    if len(candidates) != 37:
        raise ValueError(f"expected 37 S3-homology candidates, found {len(candidates)}")
    for row in candidates.values():
        if row["rational_betti_numbers"] != S3_BETTI:
            raise ValueError("homology candidate has incorrect rational Betti numbers")
        if any(betti != S3_BETTI for betti in row["finite_field_betti_numbers"].values()):
            raise ValueError("homology candidate has incorrect finite-field Betti numbers")
    return candidates, document


def make_certificate(identifier: int, facets: set, homology_row: dict, max_2_3_moves: int) -> dict:
    path = no_insertion_path(facets, max_2_3_moves)
    if path is None:
        raise RuntimeError(f"no no-insertion simplification found for record {identifier}")
    final = set(facets)
    move_rows = []
    maximum_vertices = len(vertices(final))
    for step, move in enumerate(path):
        final = apply_move(final, move)
        maximum_vertices = max(maximum_vertices, len(vertices(final)))
        move_rows.append(move.to_json(step=step))
    if maximum_vertices > 12:
        raise AssertionError(f"record {identifier} left G(12)")
    if not simplex_boundary(final):
        raise AssertionError(f"record {identifier} did not reach a simplex boundary")
    certificate = certificate_document(
        3,
        facets,
        move_rows,
        final,
        {
            "source_id": f"12seeds.txt#{identifier}",
            "purpose": "no-insertion simplification of a high-valence twelve-vertex S^3 seed",
            "one_four_moves": 0,
            "maximum_vertices": maximum_vertices,
            "seed_definition": "no legal 3-2 move",
            "rational_betti_numbers": homology_row["rational_betti_numbers"],
            "finite_field_betti_numbers": homology_row["finite_field_betti_numbers"],
        },
    )
    verify_certificate(facets, certificate, paranoid=True)
    return certificate


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, default=HERE.parents[2] / "seeds/12/12seeds.txt")
    ap.add_argument(
        "--homology-results",
        type=Path,
        default=HERE.parents[2] / "seeds/12/homology-results.json",
    )
    ap.add_argument("--output-directory", type=Path, default=HERE)
    ap.add_argument("--max-2-3-moves", type=int, default=3)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    source_hash = sha256(args.source)
    homology_candidates, homology_document = load_homology_candidates(
        args.homology_results, source_hash
    )
    rows_by_id: dict[int, dict] = {}
    facets_by_id: dict[int, set] = {}
    for identifier, facets, degree_sequence in load_records(args.source):
        rows_by_id[identifier] = checked_seed_row(identifier, facets, degree_sequence)
        facets_by_id[identifier] = facets
    if len(rows_by_id) != 134:
        raise ValueError(f"expected 134 entries in {args.source}, found {len(rows_by_id)}")
    if set(homology_candidates) - set(rows_by_id):
        raise ValueError("homology candidates include IDs missing from the source")

    candidate_rows = []
    certificates: dict[str, dict] = {}
    for identifier in sorted(homology_candidates):
        row = dict(rows_by_id[identifier])
        row["rational_betti_numbers"] = homology_candidates[identifier]["rational_betti_numbers"]
        row["finite_field_betti_numbers"] = homology_candidates[identifier][
            "finite_field_betti_numbers"
        ]
        candidate_rows.append(row)
        certificates[str(identifier)] = make_certificate(
            identifier,
            facets_by_id[identifier],
            homology_candidates[identifier],
            args.max_2_3_moves,
        )
        if not args.quiet:
            print(f"seed {identifier}: {certificates[str(identifier)]['num_moves']} moves", flush=True)

    print(
        f"identified 134 high-valence seeds and certified {len(candidate_rows)} S3 seeds",
        flush=True,
    )
    if args.dry_run:
        return
    args.output_directory.mkdir(parents=True, exist_ok=True)
    source_metadata = {
        "file": args.source.name,
        "sha256": source_hash,
        "entries": len(rows_by_id),
    }
    raw_document = {
        "format": "high-valence-12-vertex-seed-list-v1",
        "source": source_metadata,
        "scope": "all high-valence 12-vertex seeds listed in the source file",
        "seed_definition": "no legal 3-2 bistellar move",
        "checks": {
            "closed_connected_3_manifold": True,
            "vertex_links_are_2_spheres": True,
            "minimum_edge_valence": 4,
            "legal_3_2_moves": 0,
        },
        "orientable_entries": [row["id"] for row in rows_by_id.values() if row["orientable"]],
        "nonorientable_entries": [row["id"] for row in rows_by_id.values() if not row["orientable"]],
        "s3_homology_over_tested_fields_entries": sorted(homology_candidates),
        "seeds": [rows_by_id[identifier] for identifier in sorted(rows_by_id)],
    }
    candidate_document = {
        "format": "high-valence-12-vertex-s3-seed-list-v1",
        "source": source_metadata,
        "selection": {
            "method": "S3-type Betti numbers over Q and the finite fields in homology-results.json",
            "homology_results_file": args.homology_results.name,
            "homology_results_sha256": sha256(args.homology_results),
            "coefficient_fields": homology_document["coefficient_fields"],
            "entries": len(candidate_rows),
        },
        "certificate_status": "Every listed seed has a verified no-1-4 path to a simplex boundary.",
        "seeds": candidate_rows,
    }
    certificate_document_out = {
        "format": "n12-high-valence-s3-seed-certificates-v1",
        "source": source_metadata,
        "candidate_seed_list": "s3-seed-list.json",
        "certificates": certificates,
    }
    (args.output_directory / "seed-list.json").write_text(json.dumps(raw_document, indent=2) + "\n")
    (args.output_directory / "s3-seed-list.json").write_text(
        json.dumps(candidate_document, indent=2) + "\n"
    )
    (args.output_directory / "simplification-certificates.json").write_text(
        json.dumps(certificate_document_out, indent=2) + "\n"
    )
    (args.output_directory / "homology-results.json").write_text(
        json.dumps(homology_document, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
