# Computational certificates

This directory is the public archival release for the computational evidence
cited in `paper/main-gpt.tex`.

The present release contains the complete one-facet-subdivision certificates
for the four known unflippable triangulated 3-spheres. See
[`unflippable-complexes/`](unflippable-complexes/) for the data and replay
instructions. The generic strict verifier is
[`verify_pachner_certificate.py`](verify_pachner_certificate.py).

Each certificate is self-contained: it records its input-complex hash, every
bistellar move, its final complex, and its final-complex hash. It can be
checked independently without reproducing the randomized searches that found
the paths.
