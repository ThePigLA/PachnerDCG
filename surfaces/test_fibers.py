import random
import time
import coalesce
from coalesce import (normalize, canonical_form, subdivide_facet,
                      components_at_level, diagonal_flips, f_vector,
                      euler_characteristic, is_closed_surface, is_orientable)
from fibers import (pack, unpack, digest, degrees, degree3_vertices,
                    contract_degree3, build_base_index, project,
                    degree_reduction_walk, walk_fibers, exhaustive_fibers, DSU)

S2_4 = normalize([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
RP2_6 = normalize([(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 1, 5),
                   (1, 2, 4), (2, 3, 5), (3, 4, 1), (4, 5, 2), (5, 1, 3)])
T2_7 = normalize([tuple(sorted((i % 7, (i + 1) % 7, (i + 3) % 7))) for i in range(7)] +
                 [tuple(sorted((i % 7, (i + 2) % 7, (i + 3) % 7))) for i in range(7)])


def check(name, cond):
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}")
    assert cond, name


def grow(tri, target_n):
    """Subdivide facets until the vertex count reaches target_n."""
    cur = tri
    n = len({v for t in cur for v in t})
    while n < target_n:
        cur = subdivide_facet(cur, cur[0], n)
        n += 1
    return cur


print("pack / unpack / digest")
check("round trip", unpack(pack(RP2_6)) == normalize(RP2_6))
check("digest is deterministic", digest(RP2_6) == digest(RP2_6))
check("digest separates", digest(RP2_6) != digest(T2_7))
check("digest is 128 bit", digest(RP2_6) < 2 ** 128)

print("3--1 move inverts 1--3 subdivision")
random.seed(3)
ok = True
for base in (S2_4, RP2_6, T2_7):
    n = len({v for t in base for v in t})
    for facet in base:
        up = subdivide_facet(base, facet, n)
        check_deg = degrees(up)[n]
        if check_deg != 3:
            ok = False
        d3 = degree3_vertices(up)
        if n not in d3:
            ok = False
        back = contract_degree3(up, n)
        if canonical_form(back) != canonical_form(base):
            ok = False
check("new vertex has degree 3 and 3--1 recovers the original", ok)

print("projection identifies the correct component downstairs")
# S^2 at n=10 has 233 classes in ONE component.  Treat each class as its own
# label; every subdivision must project back to the label it came from.
level10 = set().union(*components_at_level({canonical_form(grow(S2_4, 10))}))
check("233 classes at n=10", len(level10) == 233)
index10 = build_base_index([{t} for t in sorted(level10)])
labels = sorted(level10)
hits_ok = True
for i, t in enumerate(labels[:40]):
    up = subdivide_facet(t, t[0], 10)
    if i not in project(up, index10):
        hits_ok = False
        break
check("each subdivision projects to its own label", hits_ok)

print("degree-reduction walk manufactures degree-3 vertices")
rng = random.Random(7)
big = grow(T2_7, 12)          # 12 vertices, chi = 0, 24 triangles
check("f-vector (12,36,24)", f_vector(big) == (12, 36, 24))
found = 0
for _ in range(20):
    start = big
    for _ in range(30):        # randomize away from the stacked start
        fl = diagonal_flips(start)
        if not fl:
            break
        start = rng.choice(fl)
    _, path = degree_reduction_walk(start, {}, rng, max_steps=120)
    if any(degree3_vertices(p) for p in path):
        found += 1
check(f"degree-3 vertices reached in {found}/20 randomized walks", found >= 18)

print("walk_fibers merges everything when the graph upstairs is connected")
# Ground truth: F_{S^2}(11) is connected, so all 233 labels must merge.
seeds = [subdivide_facet(t, t[0], 10) for t in labels]
t0 = time.time()
dsu, wit = walk_fibers(seeds, index10, len(labels), trials=6, max_steps=150,
                       seed=1, verbose=False)
dt = time.time() - t0
check(f"all 233 labels merged into {dsu.n_classes} class ({dt:.1f}s)",
      dsu.n_classes == 1)
check("witnesses recorded", len(wit) >= 232)

print("exhaustive_fibers agrees, on a smaller case")
# RP^2 at n=8: 16 classes in one component; F_{RP^2}(9) is connected.
rp2_8 = set().union(*components_at_level({canonical_form(grow(RP2_6, 8))}))
check("16 classes at n=8", len(rp2_8) == 16)
idx8 = build_base_index([{t} for t in sorted(rp2_8)])
seeds8 = [subdivide_facet(t, t[0], 8) for t in sorted(rp2_8)]
dsu8, explored, seen = exhaustive_fibers(seeds8, idx8, 16, verbose=False)
check(f"merged to {dsu8.n_classes} class (explored {explored}, seen {seen})",
      dsu8.n_classes == 1)

dsu8b, witb = walk_fibers(seeds8, idx8, 16, trials=6, max_steps=100, seed=2,
                          verbose=False)
check("walk mode agrees with exhaustive mode",
      dsu8b.n_classes == dsu8.n_classes)

print("memory: cache stays bounded")
coalesce._CANON_CACHE.clear()
cur = grow(T2_7, 11)
for _ in range(400):
    for nb in diagonal_flips(cur):
        canonical_form(nb)
    fl = diagonal_flips(cur)
    if not fl:
        break
    cur = fl[0]
check(f"cache size {len(coalesce._CANON_CACHE)} <= {coalesce._CANON_CACHE_MAX}",
      len(coalesce._CANON_CACHE) <= coalesce._CANON_CACHE_MAX)

print("\nall tests passed")
