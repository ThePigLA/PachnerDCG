#!/usr/bin/env python3
"""Regenerate the 11-vertex S^3 seed release from the local census files."""
from __future__ import annotations

import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "n10"))

from build_n10_s3_seed_certificates import main  # noqa: E402


if __name__ == "__main__":
    sys.argv = [
        sys.argv[0],
        "--source-directory",
        str(HERE.parents[2] / "seeds/11"),
        "--census-file",
        "3_manifolds_11_seeds.txt",
        "--type-file",
        "3_manifolds_11_seeds_type.txt",
        "--vertex-count",
        "11",
        "--output-directory",
        str(HERE),
        "--quiet",
        *sys.argv[1:],
    ]
    main()
