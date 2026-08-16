#!/usr/bin/env python3
"""Strictly replay every eleven-vertex S^3 seed certificate."""
from __future__ import annotations

import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
CERTIFICATES = HERE.parents[1]
sys.path.insert(0, str(CERTIFICATES))
sys.path.insert(0, str(HERE.parent / "n10"))

from verify_pachner_certificate import normalize_facets, verify_certificate  # noqa: E402
from build_n10_s3_seed_certificates import enumerate_moves, simplex_boundary  # noqa: E402


def main() -> None:
    seed_document = json.loads((HERE / "seed-list.json").read_text())
    certificate_document = json.loads((HERE / "simplification-certificates.json").read_text())
    seeds = seed_document["seeds"]
    certificates = certificate_document["certificates"]
    identifiers = {str(row["id"]) for row in seeds}
    assert identifiers == set(certificates)
    assert len(seeds) == 222

    for row in seeds:
        identifier = str(row["id"])
        facets = normalize_facets(row["facets"], 3)
        assert not enumerate_moves(facets, 3), identifier
        certificate = certificates[identifier]
        assert certificate["metadata"]["one_four_moves"] == 0
        assert all(move["type"] != "1-4" for move in certificate["moves"])
        final, _ = verify_certificate(facets, certificate, paranoid=True)
        assert simplex_boundary(final), identifier
    print(f"PASS: replayed {len(seeds)} eleven-vertex S^3 seed certificates without 1-4 moves")


if __name__ == "__main__":
    main()
