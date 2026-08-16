#!/usr/bin/env python3
"""Replay every one-step unflippable-sphere certificate."""
from __future__ import annotations

import json
import sys
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from verify_pachner_certificate import load_facets, vertices, verify_certificate  # noqa: E402


EXPECTED = {
    "u16-plus": 1229,
    "u20-plus": 4293,
    "u1-21-plus": 2658,
    "u2-21-plus": 1840,
}


def simplex_boundary(facets: set[frozenset[int]]) -> bool:
    vertex_set = vertices(facets)
    return len(vertex_set) == 5 and facets == {
        frozenset(face) for face in combinations(sorted(vertex_set), 4)
    }


def main() -> None:
    for name, expected_moves in EXPECTED.items():
        directory = HERE / name
        facets = load_facets(directory / "input.json", 3)
        certificate = json.loads((directory / "simplification-certificate.json").read_text())
        final, counts = verify_certificate(facets, certificate, paranoid=True)
        assert sum(counts.values()) == expected_moves
        assert counts["1-4"] == 1
        assert simplex_boundary(final)
        print(f"PASS: {name}: {expected_moves} moves; {dict(sorted(counts.items()))}")


if __name__ == "__main__":
    main()
