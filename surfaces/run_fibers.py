"""
Compute the fibers of s_n : Pi_M(n) -> Pi_M(n+1) from a Lutz-format file.

    python run_fibers.py <file> <chi> <orientable:0|1> [mode] [trials]

    mode = walk       (default) randomized degree-reduction walks.
                      Cheap, low memory, certifies positive merges.
    mode = exhaustive Full component exhaustion with hash-only bookkeeping.
                      The only mode that can prove classes do NOT merge.

Example (nonorientable genus 4, n = 9):
    python run_fibers.py manifold_lex_d2_n9_o0.txt -2 0

Example (orientable genus 6, n = 12):
    python run_fibers.py manifold_lex_d2_n12_o1.txt -10 1 walk 60
"""

import sys
import time

from coalesce import (parse_triangulations_file, select_surface, canonical_form,
                      components_at_level, subdivide_facet, vertices_of,
                      f_vector)
from fibers import build_base_index, walk_fibers, exhaustive_fibers


def main(path, chi, orientable, mode="walk", trials=40):
    print(f"Parsing {path} ...")
    raw = parse_triangulations_file(path)
    print(f"  {len(raw)} blocks extracted")

    target = select_surface(raw, chi=chi, orientable=orientable)
    print(f"  {len(target)} closed surfaces with chi={chi}, orientable={orientable}")
    if not target:
        return
    print(f"  f-vector {f_vector(target[0])}")

    level = {canonical_form(t) for t in target}
    n = len(vertices_of(next(iter(level))))
    print(f"  {len(level)} isomorphism classes at n={n}")

    print(f"\nComponents of F(n={n}) ...")
    t0 = time.time()
    comps = components_at_level(level)
    comps.sort(key=len, reverse=True)
    print(f"  {len(comps)} components, sizes {[len(c) for c in comps]} "
          f"({time.time() - t0:.1f}s)")

    base_index = build_base_index(comps)

    # one facet subdivision of one representative per component suffices
    seeds = []
    for c in comps:
        rep = min(c)
        seeds.append(subdivide_facet(rep, rep[0], n))
    print(f"\nSeeds at n={n+1}: {len(seeds)}")

    t0 = time.time()
    if mode == "exhaustive":
        dsu, explored, seen = exhaustive_fibers(seeds, base_index, len(comps))
        print(f"  explored {explored}, seen {seen} ({time.time() - t0:.1f}s)")
    else:
        dsu, wit = walk_fibers(seeds, base_index, len(comps), trials=trials)
        print(f"  {len(wit)} merge witnesses found ({time.time() - t0:.1f}s)")

    fibers = dsu.classes()
    fibers.sort(key=len, reverse=True)
    print(f"\nFibers of s_{n}: {len(fibers)} distinct components at n={n+1}")
    print(f"  fiber sizes: {[len(f) for f in fibers]}")
    if len(fibers) == 1:
        print(f"\nAll {len(comps)} components of F({n}) coalesce at n={n+1}.")
        print("Absorption time = 1.")
    else:
        print(f"\n{len(fibers)} classes remain after one stabilization.")
        if mode != "exhaustive":
            print("NOTE: walk mode cannot prove non-merging. Re-run the")
            print("      remaining representatives in exhaustive mode, or with")
            print("      more trials, before drawing any conclusion.")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1], int(sys.argv[2]), bool(int(sys.argv[3])),
         sys.argv[4] if len(sys.argv) > 4 else "walk",
         int(sys.argv[5]) if len(sys.argv) > 5 else 40)
