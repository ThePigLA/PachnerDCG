# Computational certificates

This directory is the public archival release for the computational evidence
cited in `paper/main.tex`.

The present release contains the complete one-facet-subdivision certificates
for the four known unflippable triangulated 3-spheres; the complete lists and
no-insertion simplification certificates for the 42 ten-vertex and 222
eleven-vertex S^3 seeds; and the raw 134-entry high-valence twelve-vertex seed
list with complete no-insertion certificates for its 37 S^3 seeds. See
[`unflippable-complexes/`](unflippable-complexes/),
[`seed-lists/n10/`](seed-lists/n10/), and
[`seed-lists/n11/`](seed-lists/n11/), and
[`seed-lists/n12/`](seed-lists/n12/) for data and replay instructions.
The generic strict verifier is
[`verify_pachner_certificate.py`](verify_pachner_certificate.py).

Each certificate is self-contained: it records its input-complex hash, every
bistellar move, its final complex, and its final-complex hash. It can be
checked independently without reproducing the randomized searches that found
the paths.

The n=12 certificates use at most twelve vertices and no (1\)-(4) moves,
so each directly certifies membership in \(\mathcal T(12)\).  The remaining
97 raw high-valence entries are not asserted to be 3-spheres.
