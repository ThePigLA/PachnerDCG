# Ten-vertex S^3 seed certificates

`seed-list.json` contains the 42 seed triangulations among the 247,882
labelled ten-vertex 3-spheres in the Lutz census. Here a seed means a
triangulation with no legal `3-2` bistellar move.

`simplification-certificates.json` maps each Lutz identifier to a strict
Pachner certificate ending at the boundary of a 4-simplex. None of the 42
paths uses a `1-4` insertion: together they use 210 `4-1`, 39 `3-2`, and 10
`2-3` moves. The longest path has 14 moves.

The original full census is deliberately not duplicated here. Its two source
filenames and SHA-256 hashes are recorded in `seed-list.json`; the compact
release retains every selected facet list and every move needed to check the
claim independently.

## Verification

From the repository root, replay all 42 certificates with:

```bash
python3 certificates/seed-lists/n10/verify_n10_s3_seed_certificates.py
```

The test confirms that every listed input has no legal `3-2` move, replays
every certificate using the generic strict verifier, checks that no `1-4`
move occurs, and confirms the terminal complex is a 4-simplex boundary.

`build_n10_s3_seed_certificates.py` reproducibly regenerates the two JSON
files from `seeds/10/3_manifolds_10_all.txt` and
`seeds/10/3_manifolds_10_s3.txt`.
