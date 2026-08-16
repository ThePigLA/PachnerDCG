import sys
import time
from pathlib import Path


# This test lives one directory below the shared surface helpers.  Keep the
# local survey modules first on sys.path, then make ``surfaces/fibers.py``
# available when the test is invoked from the repository root.
SURFACES = Path(__file__).resolve().parents[1]
if str(SURFACES) not in sys.path:
    sys.path.append(str(SURFACES))

from coalesce import (normalize, canonical_form, subdivide_facet,
                      components_at_level, f_vector, diagonal_flips)
from fibers import degrees, contract_degree3
from survey import (is_subdivision, component_is_newborn, newborn_analysis,
                    components_of_level, stabilization_map, survey_surface,
                    min_vertices, surface_name)

S2_4 = normalize([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
RP2_6 = normalize([(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 1, 5),
                   (1, 2, 4), (2, 3, 5), (3, 4, 1), (4, 5, 2), (5, 1, 3)])
T2_7 = normalize([tuple(sorted((i % 7, (i + 1) % 7, (i + 3) % 7))) for i in range(7)] +
                 [tuple(sorted((i % 7, (i + 2) % 7, (i + 3) % 7))) for i in range(7)])


def check(name, cond):
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}")
    assert cond, name


def closure(seed, lo, hi):
    cur = canonical_form(seed)
    n = f_vector(cur)[0]
    out, frontier = {}, None
    while n < lo:
        cur = canonical_form(subdivide_facet(cur, cur[0], n))
        n += 1
    frontier = {cur}
    while n <= hi:
        allt = sorted(set().union(*components_at_level(frontier)))
        out[n] = allt
        if n == hi:
            break
        frontier = {canonical_form(subdivide_facet(t, t[0], n)) for t in allt}
        n += 1
    return out


print("is_subdivision detects the subdividing vertex")
for base in (S2_4, RP2_6, T2_7):
    n = f_vector(base)[0]
    check(f"{f_vector(base)} minimal: not a subdivision",
          is_subdivision(base) is None)
    up = canonical_form(subdivide_facet(base, base[0], n))
    v = is_subdivision(up)
    check(f"{f_vector(up)} is a subdivision (witness vertex {v})",
          v is not None)
    check("witness has degree 3", degrees(up)[v] == 3)
    check("3--1 move is simplicial", contract_degree3(up, v) is not None)

print("\nminimal levels: every component newborn (no degree-3 vertex possible)")
# At the minimum vertex number a degree-3 vertex would yield a smaller
# triangulation of the same surface, so all components must be newborn.
for seed, name in ((S2_4, "S^2"), (RP2_6, "RP^2"), (T2_7, "T^2")):
    lvl = [canonical_form(seed)]
    comps, _ = components_of_level(lvl)
    na = newborn_analysis(comps, progress_every=0)
    check(f"{name} at minimum: {sum(a['newborn'] for a in na)}/{len(comps)} newborn",
          all(a["newborn"] for a in na))

print("\nabove minimum: S^2 levels 5..12 have no newborn component")
s2 = closure(S2_4, 5, 12)
for n, lvl in s2.items():
    comps, _ = components_of_level(lvl)
    na = newborn_analysis(comps, progress_every=0)
    nb = sum(a["newborn"] for a in na)
    check(f"S^2 n={n}: {len(lvl)} triangulations, {nb} newborn", nb == 0)

print("\nabove minimum: RP^2 levels 7..9 have no newborn component")
rp2 = closure(RP2_6, 7, 9)
for n, lvl in rp2.items():
    comps, _ = components_of_level(lvl)
    na = newborn_analysis(comps, progress_every=0)
    nb = sum(a["newborn"] for a in na)
    check(f"RP^2 n={n}: {len(lvl)} triangulations, {nb} newborn", nb == 0)

print("\nthe two newborn computations agree (independent code paths)")
# fiber-based: components at n+1 not hit by s_n
# degree-based: components at n+1 with no degree-3 vertex anywhere
agreements = 0
for levels in (s2, rp2):
    ns = sorted(levels)
    for n in ns[:-1]:
        comps, _ = components_of_level(levels[n])
        comps1, idx1 = components_of_level(levels[n + 1])
        img = set(stabilization_map(comps, n, idx1))
        by_fiber = sorted(j for j in range(len(comps1)) if j not in img)
        na1 = newborn_analysis(comps1, progress_every=0)
        by_degree = sorted(a["id"] for a in na1 if a["newborn"])
        check(f"n={n}->{n+1}: fibers {by_fiber} == degree-3 {by_degree}",
              by_fiber == by_degree)
        agreements += 1
check(f"{agreements} independent agreements", agreements >= 9)

print("\nneighborly triangulations are newborn (no vertex of degree 3)")
# T^2 at n=7 is neighborly: every vertex has degree 6
d = degrees(T2_7)
check(f"T^2_7 degrees all 6", set(d.values()) == {6})
check("T^2_7 not a subdivision", is_subdivision(T2_7) is None)
# RP^2 at n=6 is neighborly: every vertex has degree 5
check(f"RP^2_6 degrees all 5", set(degrees(RP2_6).values()) == {5})
check("RP^2_6 not a subdivision", is_subdivision(RP2_6) is None)

print("\npositive case: a non-subdivision whose COMPONENT is still not newborn")
# The octahedron sits at n=6, above S^2's minimum of 4, and has every vertex
# of degree 4 -- so it is not itself a stellar subdivision.  But it is
# flip-connected to the 6-vertex stacked sphere, which is one, so its
# component is in the image of s_5 and is NOT newborn.  This separates the
# per-triangulation test from the per-component test.
OCTA = normalize([(0,1,2),(0,2,4),(0,4,3),(0,3,1),
                  (5,2,1),(5,4,2),(5,3,4),(5,1,3)])
BIPYR = normalize([(0,1,3),(1,2,3),(2,0,3),(0,1,4),(1,2,4),(2,0,4)])
STACKED6 = canonical_form(subdivide_facet(BIPYR, BIPYR[0], 5))
OCTA = canonical_form(OCTA)
check("octahedron degrees all 4", set(degrees(OCTA).values()) == {4})
check("octahedron is NOT itself a subdivision", is_subdivision(OCTA) is None)
check("stacked 6-sphere IS a subdivision", is_subdivision(STACKED6) is not None)
nb_alone, wit_alone, md = component_is_newborn([OCTA])
check("octahedron alone would be newborn (min degree 4)",
      nb_alone is True and md == 4)
nb_comp, wit, _ = component_is_newborn([OCTA, STACKED6])
check("its real component is NOT newborn", nb_comp is False)
check("witness is the stacked sphere", wit is not None and wit[0] == STACKED6)
# and the fiber computation must agree
comps5, _ = components_of_level(s2[5])
comps6, idx6 = components_of_level(s2[6])
img = set(stabilization_map(comps5, 5, idx6))
check("fiber view: no newborn at n=6",
      [j for j in range(len(comps6)) if j not in img] == [])

print("\nshort-circuit: large components settle fast")
comps12, _ = components_of_level(s2[12])
t0 = time.time()
na = newborn_analysis(comps12, progress_every=0)
dt = time.time() - t0
check(f"7595-member component classified in {dt*1000:.1f} ms", dt < 1.0)

print("\nsurvey records carry the new fields")
recs = survey_surface(2, True, {n: s2[n] for n in (10, 11)}, verbose=False)
r = recs[0]
for key in ("at_minimum_n", "min_vertices", "n_newborn_here",
            "newborn_here_sizes", "min_degree_by_component",
            "newborn_by_degree_criterion", "newborn_criteria_agree"):
    check(f"record has {key}", key in r)
check("cross-check agrees", r["newborn_criteria_agree"] is True)
check("S^2 min_vertices is 4", r["min_vertices"] == 4)
check("n=10 is not minimal", r["at_minimum_n"] is False)

print("\nall tests passed")
