#!/usr/bin/env python3
"""Replay every public no-insertion certificate for an n=12 S^3 seed."""
from __future__ import annotations

import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERTIFICATES = HERE.parents[1]
sys.path.insert(0, str(CERTIFICATES))
sys.path.insert(0, str(HERE.parent / "n10"))

from build_n10_s3_seed_certificates import simplex_boundary  # noqa: E402
from verify_pachner_certificate import normalize_facets, verify_certificate  # noqa: E402


def main() -> None:
    seeds = json.loads((HERE / "s3-seed-list.json").read_text())
    certificates = json.loads((HERE / "simplification-certificates.json").read_text())
    homology = json.loads((HERE / "homology-results.json").read_text())
    assert seeds["format"] == "high-valence-12-vertex-s3-seed-list-v1"
    assert certificates["format"] == "n12-high-valence-s3-seed-certificates-v1"
    assert len(seeds["seeds"]) == 37
    assert len(certificates["certificates"]) == 37
    expected_ids = homology["s3_homology_over_tested_fields_ids"]
    assert [row["id"] for row in seeds["seeds"]] == expected_ids
    for row in seeds["seeds"]:
        identifier = str(row["id"])
        certificate = certificates["certificates"][identifier]
        assert certificate["metadata"]["one_four_moves"] == 0
        assert certificate["metadata"]["maximum_vertices"] == 12
        assert all(move["type"] != "1-4" for move in certificate["moves"])
        final, counts = verify_certificate(
            normalize_facets(row["facets"], 3), certificate, paranoid=True
        )
        assert counts["1-4"] == 0
        assert simplex_boundary(final)
    print(f"replayed {len(seeds['seeds'])} high-valence n=12 S^3 seed certificates")


if __name__ == "__main__":
    main()
