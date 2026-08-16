# One-step certificates for the known unflippable 3-spheres

This directory contains the four certificates used for the one-step merger
result in `paper/main-gpt.tex`. Each `simplification-certificate.json` starts
with the indicated unflippable 3-sphere, performs exactly one `1-4` move, and
ends at the boundary of the 4-simplex using only `2-3`, `3-2`, and `4-1`
moves. No 4-dimensional complex or certificate is used here.

| complex | input vertices | certificate moves |
| --- | ---: | ---: |
| `u16-plus` (`U(16)^+`) | 16 | 1,229 |
| `u20-plus` (`U(20)^+`) | 20 | 4,293 |
| `u1-21-plus` (`U_1(21)^+`) | 21 | 2,658 |
| `u2-21-plus` (`U_2(21)^+`) | 21 | 1,840 |

## Contents of each complex directory

- `input.json` is the original unflippable 3-sphere as a tetrahedron list.
- `simplification-certificate.json` is the complete proof certificate.

## Verification

From `certificates/`, replay any certificate using the strict verifier:

```bash
python3 verify_pachner_certificate.py \
  unflippable-complexes/u20-plus/input.json \
  unflippable-complexes/u20-plus/simplification-certificate.json \
  --paranoid
```

The verifier reconstructs each move and checks its exact live star, the
absence of the opposite face, and the final-complex hash. It also checks the
closed connected orientable pseudomanifold conditions when `--paranoid` is
used. Replay every certificate at once with:

```bash
python3 unflippable-complexes/test_unflippable_certificates.py
```

The mathematical proof objects are the four complete JSON certificates. The
randomized exploratory searches used to find them are not part of this public
release and are not needed to verify them.
