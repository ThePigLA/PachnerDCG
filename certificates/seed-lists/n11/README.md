# Eleven-vertex S^3 seed certificates

`seed-list.json` contains the 222 (S^3) entries in the 11-vertex seed
census. A seed is a triangulation with no legal `3-2` bistellar move.

`simplification-certificates.json` maps every source identifier to a strict
Pachner certificate ending at the boundary of a 4-simplex. No path uses a
`1-4` insertion. Together the paths use 1,332 `4-1`, 268 `3-2`, and 49
`2-3` moves; the longest certificate has 17 moves.

The source filenames and SHA-256 hashes are recorded in `seed-list.json`.
The release retains every selected facet list and every move needed to verify
the result independently, without copying the full input census.

## Verification

From the repository root, replay every certificate with:

```bash
python3 certificates/seed-lists/n11/verify_n11_s3_seed_certificates.py
```

The test confirms seedhood, replays every path using the generic strict
verifier, checks that no `1-4` move occurs, and confirms every terminal
complex is a 4-simplex boundary.

`build_n11_s3_seed_certificates.py` regenerates the release files from
`seeds/11/3_manifolds_11_seeds.txt` and
`seeds/11/3_manifolds_11_seeds_type.txt`.
