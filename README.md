# PachnerDCG

Computational data and verification code for *Component Trees of Pachner
Graphs with Applications to Triangulated 3-Spheres* by Vance Faber and Michael
Murphy.

The paper studies fixed-vertex Pachner graphs and their components under facet
subdivision.  Alongside the general component-stabilization theorem, this
repository provides independently checkable computational certificates for
small triangulated 3-spheres and code supporting the surface stabilization
experiments.

## Contents

| Directory | Contents |
| --- | --- |
| [`paper/`](paper/) | LaTeX source for the manuscript. |
| [`certificates/`](certificates/) | Self-contained Pachner-move certificates and strict replay tools. |
| [`surfaces/`](surfaces/) | Code, summaries, and checkpoints for the surface flip-graph experiments. |

The certificate release includes:

- the 42 ten-vertex and 222 eleven-vertex \(S^3\) seed lists, with
  no-\(1\)-\(4\) simplifications to the boundary of a 4-simplex;
- the 134-entry high-valence twelve-vertex seed list, with the 37 certified
  \(S^3\) entries and their simplification certificates; and
- one-facet-subdivision certificates for the four known unflippable
  3-spheres \(U(16)\), \(U(20)\), \(U_1(21)\), and \(U_2(21)\).

See [`certificates/README.md`](certificates/README.md) for the data layout and
certificate semantics.

## Requirements

The public certificate checkers and the surface code use the Python standard
library and have been tested with Python 3.11. No package installation is
needed for the commands below.

## Verify the 3-sphere certificates

From the repository root, run:

```bash
python3 certificates/unflippable-complexes/test_unflippable_certificates.py
python3 certificates/seed-lists/n10/verify_n10_s3_seed_certificates.py
python3 certificates/seed-lists/n11/verify_n11_s3_seed_certificates.py
python3 certificates/seed-lists/n12/verify_n12_high_valence_seed_list.py
python3 certificates/seed-lists/n12/verify_n12_s3_seed_certificates.py
```

Each replay checks every move's exact local legality, the intermediate
closed-manifold conditions, and the claimed terminal complex.

## Surface experiments

The surface implementation and its experiment summaries are described in
[`surfaces/README.md`](surfaces/README.md). The lightweight self-checks can be
run with:

```bash
python3 surfaces/test_coalesce.py
python3 surfaces/test_fibers.py
python3 surfaces/test_fiberstate.py
python3 surfaces/survey/test_survey.py
```

The large third-party Lutz census files are intentionally not stored in Git.
The surface README identifies the required inputs and the commands used to
reproduce the recorded experiment summaries.

## Citation and archival release

The development repository is hosted at
[github.com/ThePigLA/PachnerDCG](https://github.com/ThePigLA/PachnerDCG).
The archival release is deposited on Zenodo, which issues two DOIs for it:

| DOI | Resolves to | Use it when |
| --- | --- | --- |
| _concept DOI: to be added_ | whichever version is newest | citing this work in general |
| [10.5281/zenodo.21980010](https://doi.org/10.5281/zenodo.21980010) | this snapshot, permanently | reproducing a specific result |

Prefer the concept DOI: it is stable across future revisions of the data and
code, so a citation made today keeps resolving to the current release.  Cite
the version DOI when a claim depends on the exact bytes of one snapshot.

Zenodo mints the concept DOI automatically the first time the deposit is
published, and it will be recorded here and in
[`CITATION.cff`](CITATION.cff) at that point.  Until the deposit is
published, neither DOI resolves.

The Zenodo release contains the immutable repository snapshot and the
additional archival data required for full reproduction.  Machine-readable
citation metadata for this repository, and for the accompanying paper, is in
[`CITATION.cff`](CITATION.cff).

## License

This repository is dedicated to the public domain under the
[CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/)
dedication.  You may use, copy, modify, distribute, and build upon its
contents for any purpose, without asking permission.
