#!/usr/bin/env python3
"""Verify the public n=12 high-valence seed-list artifact.

Every property is checked against ``seed-list.json`` alone.  The checker
reads no file outside this repository.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERTIFICATES = HERE.parents[1]
sys.path.insert(0, str(CERTIFICATES))
sys.path.insert(0, str(HERE.parent / "n10"))

from build_n10_s3_seed_certificates import enumerate_moves  # noqa: E402
from build_n12_high_valence_seed_certificates import edge_valences, vertex_links_are_spheres
from verify_pachner_certificate import (
    check_connected,
    check_orientable,
    normalize_facets,
    summarize_complex,
    vertices,
)


def main() -> None:
    document = json.loads((HERE / "seed-list.json").read_text())
    assert document["format"] == "high-valence-12-vertex-seed-list-v1"
    assert len(document["seeds"]) == 134
    for row in document["seeds"]:
        facets = normalize_facets(row["facets"], 3)
        assert len(vertices(facets)) == 12
        assert summarize_complex(facets)["closed_pseudomanifold"]
        assert check_connected(facets)
        assert vertex_links_are_spheres(facets)
        assert row["orientable"] == check_orientable(facets)
        valences = edge_valences(facets)
        assert min(valences.values()) >= 4
        assert enumerate_moves(facets, 3) == []
        assert row["edge_valence_range"] == [min(valences.values()), max(valences.values())]
        assert row["legal_3_2_moves"] == 0
    assert document["orientable_entries"] == [
        row["id"] for row in document["seeds"] if row["orientable"]
    ]
    assert document["nonorientable_entries"] == [
        row["id"] for row in document["seeds"] if not row["orientable"]
    ]
    print(f"verified {len(document['seeds'])} high-valence n=12 seeds")


if __name__ == "__main__":
    main()
