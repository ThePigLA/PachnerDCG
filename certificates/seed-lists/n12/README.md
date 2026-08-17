# High-valence twelve-vertex seed list

`seed-list.json` is the raw 134-entry high-valence (12)-vertex seed list
used in the computation.  It was extracted from the raw enumeration output
`12seeds.txt`, whose name, entry count, and SHA-256 digest are recorded in the
`source` field of the JSON artifact as a provenance record.  That raw file is
not redistributed here and nothing needs it: `seed-list.json` carries the
facet lists themselves, and the checker below re-derives every asserted
property from them.

`homology-results.json` records the rational and finite-field homology
screening.  `s3-seed-list.json` contains its 37 (S^3)-homology candidates,
and `simplification-certificates.json` gives a direct no-(1\)-(4) path to
the boundary of a 4-simplex for every one.  These paths independently certify
that the 37 included complexes are (S^3)'s.

Run the independent checker from the repository root:

```bash
python3 certificates/seed-lists/n12/verify_n12_high_valence_seed_list.py
python3 certificates/seed-lists/n12/verify_n12_s3_seed_certificates.py
```

For every entry, it verifies the facet list is a connected, closed
3-dimensional pseudomanifold on twelve vertices; that all edge valences are
at least four; and that no legal \(3\)-\(2\) move exists.  Thus the artifact
certifies the claimed high-valence seed property.  The checker also determines
and records orientability.

The raw list is not a list of 134 3-spheres: 11 entries are nonorientable and
many others have homology inconsistent with (S^3).  They are retained in the
raw artifact to document the high-valence enumeration, but no path to
\(\mathcal T(12)\) is asserted for them.
