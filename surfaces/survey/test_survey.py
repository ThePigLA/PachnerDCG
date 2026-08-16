import io
import os
import time

from coalesce import (normalize, canonical_form, subdivide_facet,
                      components_at_level, diagonal_flips, f_vector,
                      euler_characteristic, is_orientable)
from survey import (surface_name, min_vertices, components_of_level,
                    stabilization_map, verify_stabilization_welldefined,
                    bfs_distances, coalescence_distance, survey_surface,
                    print_table, print_detail, load_levels)

S2_4 = normalize([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
RP2_6 = normalize([(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 1, 5),
                   (1, 2, 4), (2, 3, 5), (3, 4, 1), (4, 5, 2), (5, 1, 3)])
T2_7 = normalize([tuple(sorted((i % 7, (i + 1) % 7, (i + 3) % 7))) for i in range(7)] +
                 [tuple(sorted((i % 7, (i + 2) % 7, (i + 3) % 7))) for i in range(7)])


def check(name, cond):
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}")
    assert cond, name


def closure_levels(seed, lo, hi):
    """
    Complete levels lo..hi by flip closure, valid only for surfaces whose
    flip graph is connected at each level (verified below against published
    counts for S^2 and RP^2).
    """
    cur = canonical_form(seed)
    n = len(f_vector(cur)) and f_vector(cur)[0]
    levels = {}
    while n < lo:
        cur = canonical_form(subdivide_facet(cur, cur[0], n))
        n += 1
    frontier = {cur}
    while n <= hi:
        comps = components_at_level(frontier)
        allt = sorted(set().union(*comps))
        levels[n] = allt
        if n == hi:
            break
        frontier = {canonical_form(subdivide_facet(t, t[0], n)) for t in allt}
        n += 1
    return levels


print("surface naming and Heawood bounds")
check("S^2", surface_name(2, True) == "S^2")
check("T^2", surface_name(0, True) == "T^2")
check("Sigma_6", surface_name(-10, True) == "Sigma_6")
check("RP^2", surface_name(1, False) == "RP^2")
check("Klein", surface_name(0, False) == "Klein")
check("N_4", surface_name(-2, False) == "N_4")
check("S^2 needs 4", min_vertices(2, True) == 4)
check("RP^2 needs 6", min_vertices(1, False) == 6)
check("T^2 needs 7", min_vertices(0, True) == 7)
check("Klein needs 8 (exception)", min_vertices(0, False) == 8)
check("N_3 needs 9 (exception)", min_vertices(-1, False) == 9)
check("N_4 needs 9", min_vertices(-2, False) == 9)
check("Sigma_6 needs 12", min_vertices(-10, True) == 12)

print("\ncomplete levels for S^2 (published counts)")
s2 = closure_levels(S2_4, 4, 10)
KNOWN_S2 = {4: 1, 5: 1, 6: 2, 7: 5, 8: 14, 9: 50, 10: 233}
for n, k in KNOWN_S2.items():
    check(f"n={n}: {len(s2[n])} == {k}", len(s2[n]) == k)

print("\ncomplete levels for RP^2 (published counts)")
rp2 = closure_levels(RP2_6, 6, 9)
KNOWN_RP2 = {6: 1, 7: 3, 8: 16, 9: 134}
for n, k in KNOWN_RP2.items():
    check(f"n={n}: {len(rp2[n])} == {k}", len(rp2[n]) == k)

print("\ncomponents_of_level rejects incomplete input")
bad = [t for t in s2[7]][:2]
try:
    components_of_level(bad)
    check("raised on incomplete level", False)
except ValueError:
    check("raised on incomplete level", True)

print("\ncomponents on complete levels")
for n in range(4, 11):
    comps, idx = components_of_level(s2[n])
    check(f"S^2 n={n}: 1 component of size {len(comps[0])}",
          len(comps) == 1 and len(comps[0]) == KNOWN_S2[n])

print("\nlift-fiber lemma verified on real data")
total_checks = total_bad = 0
for n in range(4, 10):
    comps, _ = components_of_level(s2[n])
    _, idx1 = components_of_level(s2[n + 1])
    c, b = verify_stabilization_welldefined(comps, n, idx1)
    total_checks += c
    total_bad += b
for n in range(6, 9):
    comps, _ = components_of_level(rp2[n])
    _, idx1 = components_of_level(rp2[n + 1])
    c, b = verify_stabilization_welldefined(comps, n, idx1)
    total_checks += c
    total_bad += b
check(f"{total_checks} facet subdivisions, {total_bad} violations",
      total_bad == 0 and total_checks == 1214)

print("\nflip distances")
d = bfs_distances(s2[7][0], set(s2[7]))
check(f"all 5 classes of F_S2(7) reached", len(d) == 5)
check("source at distance 0", d[s2[7][0]] == 0)
ecc = max(d.values())
check(f"eccentricity {ecc} >= 1", ecc >= 1)
d10 = bfs_distances(s2[10][0], set(s2[10]))
check(f"all 233 classes of F_S2(10) reached", len(d10) == 233)
print(f"       eccentricity of that vertex in F_S2(10): {max(d10.values())}")

print("\ncoalescence distance is 0 on singleton fibers")
comps7, _ = components_of_level(s2[7])
dist, pair = coalescence_distance([0], comps7, 7, set(s2[8]))
check("singleton fiber has distance 0", dist == 0 and pair is None)

print("\nfull survey on S^2")
recs = survey_surface(2, True, s2, verbose=False)
check("one record per level", len(recs) == 7)
for r in recs:
    check(f"n={r['n']}: 1 component", r["n_components"] == 1)
    if "image_size" in r:
        check(f"n={r['n']}: image is 1 component", r["image_size"] == 1)
        check(f"n={r['n']}: largest fiber 1", r["largest_fiber"] == 1)
        check(f"n={r['n']}: 0 newborn", r["n_newborn"] == 0)
        check(f"n={r['n']}: 0 violations", r["welldefined_violations"] == 0)

print("\nfull survey on RP^2")
recs2 = survey_surface(1, False, rp2, verbose=False)
check("one record per level", len(recs2) == 4)
check("all levels connected", all(r["n_components"] == 1 for r in recs2))

print("\nnewborn detection on a synthetic truncation")
# Drop level 7 down to a single component's worth and pretend level 8 is
# complete: every component of level 8 not hit by s_7 must be reported.
fake = {7: s2[7], 8: s2[8]}
comps7, _ = components_of_level(fake[7])
comps8, idx8 = components_of_level(fake[8])
img = stabilization_map(comps7, 7, idx8)
newborn = [j for j in range(len(comps8)) if j not in set(img)]
check("S^2 has no newborn components at n=8", newborn == [])

print("\ntable and detail rendering")
buf = io.StringIO()
print_table(recs + recs2, stream=buf)
txt = buf.getvalue()
check("header present", "surface" in txt and "newborn" in txt)
check("S^2 rows present", txt.count("S^2") == 7)
check("RP^2 rows present", txt.count("RP^2") == 4)
buf2 = io.StringIO()
print_detail(recs2, stream=buf2)
check("detail mentions lift-fiber check", "lift-fiber check" in buf2.getvalue())

print("\nload_levels round trip through a Lutz-format file")
with open("/tmp/lvl9.txt", "w") as fh:
    for i, t in enumerate(rp2[9]):
        fh.write(f"#{i+1}= " + str([list(x) for x in t]).replace(" ", "") + "\n")
lv = load_levels(["/tmp/lvl9.txt"], verbose=False)
check("one surface found", len(lv) == 1)
key = list(lv)[0]
check("chi=1, nonorientable", key == (1, False))
check("134 at n=9", len(lv[key][9]) == 134)

print("\nall tests passed")
