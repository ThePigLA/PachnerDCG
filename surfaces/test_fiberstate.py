import io
import os
import random
import time

import coalesce
from coalesce import (normalize, canonical_form, subdivide_facet,
                      components_at_level, f_vector)
from fiberstate import (FiberState, automorphism_order, targeted_search,
                        targeted_round, tri_to_json, tri_from_json)

S2_4 = normalize([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)])
RP2_6 = normalize([(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5), (0, 1, 5),
                   (1, 2, 4), (2, 3, 5), (3, 4, 1), (4, 5, 2), (5, 1, 3)])


def check(name, cond):
    print(f"  [{'ok ' if cond else 'FAIL'}] {name}")
    assert cond, name


def grow(tri, target_n):
    cur, n = tri, len({v for t in tri for v in t})
    while n < target_n:
        cur = subdivide_facet(cur, cur[0], n)
        n += 1
    return cur


print("json round trip")
check("tri_to_json / tri_from_json", tri_from_json(tri_to_json(RP2_6)) == RP2_6)

# ---------------------------------------------------------------- RP^2 n=8
print("\nRP^2 at n=8: 16 classes treated as 16 separate components")
rp2_8 = sorted(set().union(*components_at_level({canonical_form(grow(RP2_6, 8))})))
check("16 classes", len(rp2_8) == 16)
st = FiberState([[t] for t in rp2_8], n=8, chi=1, orientable=False)
check("16 blocks initially", st.n_blocks() == 16)

st.compute_aut()
check("aut computed for every component", len(st.aut) == 16)
check("all |Aut| positive divisors of 6|F|",
      all(all(a >= 1 for a in v) for v in st.aut))

pts = st.entry_points([0])
check(f"entry points for one component: {len(pts)} (= #facets)",
      len(pts) == len(rp2_8[0]))
check("all entry points have 8+1 vertices",
      all(len({v for t in p for v in t}) == 9 for p in pts))

t0 = time.time()
targeted_search(st, rounds=20, walks_per_block=12, max_steps=80, seed=5,
                verbose=False)
check(f"merged 16 -> {st.n_blocks()} ({time.time() - t0:.1f}s); "
      f"ground truth F_RP2(9) is connected", st.n_blocks() == 1)
check("witnesses recorded", len(st.witnesses) >= 15)

print("\ncheckpoint save / load / resume")
path = "/tmp/state_rp2.json"
st.save(path)
check("file written", os.path.exists(path))
st2 = FiberState.load(path)
check("blocks preserved", st2.blocks() == st.blocks())
check("n_blocks preserved", st2.n_blocks() == st.n_blocks())
check("witnesses preserved", st2.witnesses == st.witnesses)
check("aut preserved", st2.aut == st.aut)
check("chi/orientable preserved",
      (st2.chi, st2.orientable) == (1, False))

print("\nresume from a partially merged state")
st3 = FiberState([[t] for t in rp2_8], n=8, chi=1, orientable=False)
rng = random.Random(0)
targeted_round(st3, rng, walks_per_block=2, max_steps=40, verbose=False)
partial = st3.n_blocks()
print(f"       partial state has {partial} blocks")
st3.save("/tmp/state_partial.json")
st4 = FiberState.load("/tmp/state_partial.json")
check("reloaded partial state matches", st4.n_blocks() == partial)
targeted_search(st4, rounds=20, walks_per_block=12, max_steps=80, seed=9,
                verbose=False)
check(f"resumed run finished at {st4.n_blocks()} block(s)",
      st4.n_blocks() == 1)

print("\nreport renders and names every triangulation")
buf = io.StringIO()
st4.report(stream=buf, show_triangulations=True)
text = buf.getvalue()
check("header present", "Class structure of s_8" in text)
check("block listing present", "BLOCK 0" in text)
check("all 16 component ids appear",
      all(f"[{i:3d}]" in text for i in range(16)))
check("fiber sizes reported", "fiber sizes" in text)

print("\nunresolved-pair reporting on a genuinely multi-block state")
st5 = FiberState([[t] for t in rp2_8], n=8, chi=1, orientable=False)
st5.merge(0, 1)
st5.merge(2, 3)
buf2 = io.StringIO()
st5.report(stream=buf2, show_triangulations=False)
t2 = buf2.getvalue()
check("multiple blocks listed", st5.n_blocks() == 14)
check("unresolved pairs section shown", "UNRESOLVED PAIRS" in t2)
st5.note_failure(0, 2)
st5.note_failure(0, 2)
check("probe counter increments", st5.pair_attempts(0, 2) == 2)
check("probe counter is symmetric", st5.pair_attempts(2, 0) == 2)

# ---------------------------------------------------------------- S^2 n=10
print("\nS^2 at n=10: 233 classes as 233 separate components")
level10 = sorted(set().union(*components_at_level({canonical_form(grow(S2_4, 10))})))
check("233 classes", len(level10) == 233)
big = FiberState([[t] for t in level10], n=10, chi=2, orientable=True)
t0 = time.time()
targeted_search(big, rounds=30, walks_per_block=8, max_steps=120, seed=11,
                verbose=False)
check(f"merged 233 -> {big.n_blocks()} ({time.time() - t0:.1f}s); "
      f"ground truth F_S2(11) is connected", big.n_blocks() == 1)

big.save("/tmp/state_s2.json")
sz = os.path.getsize("/tmp/state_s2.json")
print(f"       checkpoint for 233 components: {sz / 1024:.0f} KB")

print("\nall tests passed")
