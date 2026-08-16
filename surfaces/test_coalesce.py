import random
import time
from coalesce import (normalize, canonical_form, diagonal_flips, subdivide_facet,
                      f_vector, euler_characteristic, is_closed_surface,
                      is_orientable, components_at_level, do_seeds_coalesce)

# ---- known small triangulations -------------------------------------------

# boundary of the tetrahedron: S^2, 4 vertices
S2_4 = normalize([(0,1,2),(0,1,3),(0,2,3),(1,2,3)])

# octahedron: S^2, 6 vertices, non-edges 0-5, 1-4, 2-3
OCTA = normalize([(0,1,2),(0,2,4),(0,4,3),(0,3,1),
                  (5,2,1),(5,4,2),(5,3,4),(5,1,3)])

# stacked 6-vertex S^2: subdivide a facet of the 5-vertex bipyramid
BIPYR = normalize([(0,1,3),(1,2,3),(2,0,3),(0,1,4),(1,2,4),(2,0,4)])

# 6-vertex RP^2 (antipodal quotient of the icosahedron)
RP2_6 = normalize([(0,1,2),(0,2,3),(0,3,4),(0,4,5),(0,1,5),
                   (1,2,4),(2,3,5),(3,4,1),(4,5,2),(5,1,3)])

# 7-vertex torus (Moebius): {i,i+1,i+3} and {i,i+2,i+3} mod 7
T2_7 = normalize([tuple(sorted((i % 7, (i+1) % 7, (i+3) % 7))) for i in range(7)] +
                 [tuple(sorted((i % 7, (i+2) % 7, (i+3) % 7))) for i in range(7)])


def relabel(tri, perm):
    return normalize([tuple(sorted(perm[v] for v in t)) for t in tri])


def check(name, cond):
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}")
    assert cond, name


print("basic surface recognition")
check("S2_4 is a closed surface", is_closed_surface(S2_4))
check("S2_4 chi == 2", euler_characteristic(S2_4) == 2)
check("OCTA f-vector (6,12,8)", f_vector(OCTA) == (6, 12, 8))
check("OCTA chi == 2", euler_characteristic(OCTA) == 2)
check("RP2_6 f-vector (6,15,10)", f_vector(RP2_6) == (6, 15, 10))
check("RP2_6 chi == 1", euler_characteristic(RP2_6) == 1)
check("RP2_6 is a closed surface", is_closed_surface(RP2_6))
check("RP2_6 nonorientable", not is_orientable(RP2_6))
check("T2_7 f-vector (7,21,14)", f_vector(T2_7) == (7, 21, 14))
check("T2_7 chi == 0", euler_characteristic(T2_7) == 0)
check("T2_7 orientable", is_orientable(T2_7))
check("OCTA orientable", is_orientable(OCTA))

print("canonical form is relabeling-invariant")
random.seed(11)
for name, tri in [("OCTA", OCTA), ("RP2_6", RP2_6), ("T2_7", T2_7), ("BIPYR", BIPYR)]:
    base = canonical_form(tri)
    n = len({v for t in tri for v in t})
    same = True
    for _ in range(30):
        p = list(range(n))
        random.shuffle(p)
        if canonical_form(relabel(tri, {i: p[i] for i in range(n)})) != base:
            same = False
            break
    check(f"{name} invariant under 30 random relabelings", same)

print("canonical form separates non-isomorphic complexes")
check("OCTA != BIPYR", canonical_form(OCTA) != canonical_form(BIPYR))

print("flips")
check("RP2_6 admits no flip (neighborly)", diagonal_flips(RP2_6) == [])
check("T2_7 admits no flip (neighborly)", diagonal_flips(T2_7) == [])
check("S2_4 admits no flip", diagonal_flips(S2_4) == [])
check("OCTA admits some flips", len(diagonal_flips(OCTA)) > 0)

print("flip components of S^2 at n=6 (should be one component of size 2)")
STACKED6 = subdivide_facet(BIPYR, BIPYR[0], 5)
check("STACKED6 f-vector (6,12,8)", f_vector(STACKED6) == (6, 12, 8))
check("STACKED6 not isomorphic to OCTA",
      canonical_form(STACKED6) != canonical_form(OCTA))
comps = components_at_level([OCTA, STACKED6])
check("one component", len(comps) == 1)
check("exactly 2 classes", len(comps[0]) == 2)

print("flip components of S^2 at n=7")
seeds7 = set()
for t in (OCTA, STACKED6):
    for f in t:
        seeds7.add(subdivide_facet(t, f, 6))
comps7 = components_at_level(seeds7)
check("connected at n=7", len(comps7) == 1)
print(f"       |F_{{S^2}}(7)| = {len(comps7[0])}  (expected 5)")
check("5 classes", len(comps7[0]) == 5)

print("coalescence driver on a case with a known answer")
# RP^2 at n=6 is a single isolated vertex; at n=7 everything must be one component
seeds = [subdivide_facet(RP2_6, RP2_6[0], 6)]
ok, explored, dsu, seen = do_seeds_coalesce(seeds, verbose=False)
check("trivially coalesced (single seed)", ok)

print("timing: canonical_form")
t0 = time.time()
N = 2000
for _ in range(N):
    canonical_form(T2_7)
dt = time.time() - t0
print(f"       {N} calls on a 14-triangle complex in {dt:.2f}s "
      f"({1e6*dt/N:.0f} us/call)")

print("\nall tests passed")
