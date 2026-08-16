"""
Checkpointable class structure for the stabilization map s_n.

What this adds over `fibers.py`
-------------------------------
* The partition of Pi_M(n) induced by s_n is stored explicitly, saved to
  disk, and printed with the triangulations in each block, so a stalled run
  can be inspected and resumed instead of restarted.

* Search is *targeted*: once blocks have formed, every triangulation in every
  component of a block gives an entry point into the same component upstairs
  (all facet subdivisions of all its members).  A block of k singleton
  components on a 44-facet surface therefore offers 44k independent starting
  points rather than one.

* Automorphism orders are recorded.  For a connected closed surface Aut acts
  freely on flags, so |Aut(T)| equals the number of flags whose canonical
  traversal reproduces the canonical string.  If the surviving blocks track
  automorphism data, that is the structural fact worth reporting.

* Failed attempts per block pair are counted, so effort is spent on the pairs
  that have been probed least.
"""

from __future__ import annotations

import json
import os
import random
import sys
import time
from collections import deque

from coalesce import (normalize, canonical_form, edge_to_triangles,
                      diagonal_flips, subdivide_facet, vertices_of, f_vector,
                      euler_characteristic, _vertex_colors)
from fibers import (DSU, degrees, degree3_vertices, contract_degree3,
                    build_base_index, project, degree_reduction_walk,
                    pack, unpack, digest, _flip_on_edge)


# --------------------------------------------------------------------------
# automorphism order from the flag count
# --------------------------------------------------------------------------

def _traverse(tri, nbr, start, t0):
    lab = {start[0]: 0, start[1]: 1, start[2]: 2}
    nxt = 3
    out = []
    seen = {t0}
    queue = deque([start])
    while queue:
        x, y, z = queue.popleft()
        key = tuple(sorted((x, y, z)))
        out.append((lab[x], lab[y], lab[z]))
        for (p, q) in ((x, y), (y, z), (z, x)):
            e = (p, q) if p < q else (q, p)
            nb = nbr[(key, e)]
            if nb in seen:
                continue
            w = next(v for v in nb if v not in (p, q))
            if w not in lab:
                lab[w] = nxt
                nxt += 1
            seen.add(nb)
            queue.append((q, p, w))
    return tuple(out)


def automorphism_order(tri) -> int:
    """
    |Aut(T)| for a connected closed surface triangulation.

    Aut acts freely on the 6|F| flags (an automorphism fixing an ordered
    triangle propagates to the identity across the connected facet-adjacency
    graph), and two flags give the same canonical traversal string exactly
    when they lie in the same Aut-orbit.  So the number of flags realising the
    minimal string is |Aut(T)|.
    """
    tri = normalize(tri)
    e2t = edge_to_triangles(tri)
    nbr = {}
    for e, ts in e2t.items():
        if len(ts) != 2:
            raise ValueError("not a closed surface")
        nbr[(ts[0], e)] = ts[1]
        nbr[(ts[1], e)] = ts[0]

    best = None
    count = 0
    for t0 in tri:
        a, b, c = t0
        for start in ((a, b, c), (a, c, b), (b, a, c),
                      (b, c, a), (c, a, b), (c, b, a)):
            s = _traverse(tri, nbr, start, t0)
            if best is None or s < best:
                best, count = s, 1
            elif s == best:
                count += 1
    return count


# --------------------------------------------------------------------------
# serialisation helpers
# --------------------------------------------------------------------------

def tri_to_json(tri):
    return [list(t) for t in tri]


def tri_from_json(obj):
    return normalize(tuple(tuple(t) for t in obj))


# --------------------------------------------------------------------------
# the state object
# --------------------------------------------------------------------------

class FiberState:
    """
    Records the partition of the components of F_M(n) induced by s_n,
    together with witnesses, attempt counts and automorphism data.
    """

    def __init__(self, comps, n, chi=None, orientable=None):
        self.n = n
        self.chi = chi
        self.orientable = orientable
        self.comps = [sorted(c) for c in comps]
        self.dsu = DSU(len(self.comps))
        self.witnesses = {}      # "i,j" -> summary dict
        self.attempts = {}       # "i,j" -> failed targeted probes
        self.aut = None          # list of |Aut| per component representative
        self.base_index = build_base_index(self.comps)

    # -- automorphism data -------------------------------------------------

    def compute_aut(self):
        self.aut = []
        for c in self.comps:
            self.aut.append(sorted({automorphism_order(t) for t in c}))
        return self.aut

    # -- partition ---------------------------------------------------------

    def blocks(self):
        """List of blocks, each a sorted list of component ids."""
        b = [sorted(x) for x in self.dsu.classes()]
        b.sort(key=lambda blk: (-len(blk), blk[0]))
        return b

    def n_blocks(self):
        return self.dsu.n_classes

    def merge(self, i, j, witness=None):
        if self.dsu.union(i, j):
            key = f"{min(i, j)},{max(i, j)}"
            if witness is not None:
                self.witnesses[key] = witness
            return True
        return False

    def note_failure(self, bi, bj):
        key = f"{min(bi, bj)},{max(bi, bj)}"
        self.attempts[key] = self.attempts.get(key, 0) + 1

    def pair_attempts(self, bi, bj):
        return self.attempts.get(f"{min(bi, bj)},{max(bi, bj)}", 0)

    # -- entry points ------------------------------------------------------

    def entry_points(self, block):
        """
        Every facet subdivision of every triangulation in every component of
        `block` lies in the SAME component of F_M(n+1), by the lift-fiber
        lemma.  Returns them all as independent starting points.
        """
        pts = []
        for i in block:
            for t in self.comps[i]:
                for facet in t:
                    pts.append(subdivide_facet(t, facet, self.n))
        return pts

    # -- reporting ---------------------------------------------------------

    def report(self, stream=sys.stdout, show_triangulations=True, max_show=None):
        w = stream.write
        blocks = self.blocks()
        w(f"\n{'=' * 70}\n")
        w(f"Class structure of s_{self.n} : Pi(n={self.n}) -> Pi(n={self.n + 1})\n")
        if self.chi is not None:
            w(f"  surface: chi = {self.chi}, "
              f"orientable = {self.orientable}\n")
        w(f"  components at n={self.n}: {len(self.comps)} "
          f"(sizes {[len(c) for c in self.comps]})\n")
        w(f"  distinct components at n={self.n + 1}: {len(blocks)}\n")
        w(f"  fiber sizes: {[len(b) for b in blocks]}\n")
        w(f"{'=' * 70}\n")

        for bi, blk in enumerate(blocks):
            auts = ""
            if self.aut:
                vals = sorted({a for i in blk for a in self.aut[i]})
                auts = f"  |Aut| in {vals}"
            w(f"\nBLOCK {bi}  ({len(blk)} component"
              f"{'s' if len(blk) != 1 else ''}, "
              f"{sum(len(self.comps[i]) for i in blk)} triangulations)"
              f"{auts}\n")
            w(f"  component ids: {blk}\n")
            if show_triangulations:
                shown = blk if max_show is None else blk[:max_show]
                for i in shown:
                    for t in self.comps[i]:
                        w(f"    [{i:3d}] {list(map(list, t))}\n")
                if max_show is not None and len(blk) > max_show:
                    w(f"    ... {len(blk) - max_show} more components\n")

        if len(blocks) > 1:
            w(f"\nUNRESOLVED PAIRS (probe count):\n")
            for a in range(len(blocks)):
                for b in range(a + 1, len(blocks)):
                    ra, rb = blocks[a][0], blocks[b][0]
                    w(f"  block {a} <-> block {b}: "
                      f"{self.pair_attempts(ra, rb)} probes\n")

        if self.witnesses:
            w(f"\n{len(self.witnesses)} merge witnesses recorded\n")
        w("\n")

    # -- persistence -------------------------------------------------------

    def save(self, path):
        obj = {
            "n": self.n,
            "chi": self.chi,
            "orientable": self.orientable,
            "comps": [[tri_to_json(t) for t in c] for c in self.comps],
            "parent": self.dsu.p,
            "n_classes": self.dsu.n_classes,
            "witnesses": self.witnesses,
            "attempts": self.attempts,
            "aut": self.aut,
        }
        tmp = path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(obj, fh)
        os.replace(tmp, path)

    @classmethod
    def load(cls, path):
        with open(path) as fh:
            obj = json.load(fh)
        comps = [[tri_from_json(t) for t in c] for c in obj["comps"]]
        st = cls(comps, obj["n"], obj.get("chi"), obj.get("orientable"))
        st.dsu.p = list(obj["parent"])
        st.dsu.n_classes = obj["n_classes"]
        st.witnesses = obj.get("witnesses", {})
        st.attempts = obj.get("attempts", {})
        st.aut = obj.get("aut")
        return st


# --------------------------------------------------------------------------
# targeted search between blocks
# --------------------------------------------------------------------------

def targeted_round(state, rng, walks_per_block=40, max_steps=300,
                   checkpoint=None, verbose=True):
    """
    One round of targeted search.

    For each block, launch walks from randomly chosen entry points (all facet
    subdivisions of all members) and union on every degree-3 projection.
    Blocks are visited least-probed first.  Returns the number of merges.
    """
    merges = 0
    blocks = state.blocks()
    # least-probed block first
    order = sorted(range(len(blocks)),
                   key=lambda b: sum(state.pair_attempts(blocks[b][0], blocks[c][0])
                                     for c in range(len(blocks)) if c != b))

    for bi in order:
        blocks = state.blocks()
        if len(blocks) == 1:
            break
        if bi >= len(blocks):
            continue
        blk = blocks[bi]
        pts = state.entry_points(blk)
        rng.shuffle(pts)
        home = blk[0]

        for k in range(min(walks_per_block, len(pts))):
            start = pts[k]
            hits, path = degree_reduction_walk(start, state.base_index, rng,
                                              max_steps=max_steps)
            for j in hits:
                if state.dsu.find(j) == state.dsu.find(home):
                    continue
                wit = {
                    "from": home,
                    "to": j,
                    "path_len": len(path),
                    "start": tri_to_json(start),
                    "end": tri_to_json(path[-1]),
                }
                if state.merge(home, j, wit):
                    merges += 1
                    if verbose:
                        print(f"    MERGE {home} <-> {j} "
                              f"(path length {len(path)})", file=sys.stderr)
                    if checkpoint:
                        state.save(checkpoint)
        # record the probe against every other block still standing
        for other in state.blocks():
            if state.dsu.find(other[0]) != state.dsu.find(home):
                state.note_failure(home, other[0])

    if checkpoint:
        state.save(checkpoint)
    return merges


def targeted_search(state, rounds=50, walks_per_block=40, max_steps=300,
                    seed=0, checkpoint=None, verbose=True):
    rng = random.Random(seed)
    stalled = 0
    for r in range(rounds):
        t0 = time.time()
        m = targeted_round(state, rng, walks_per_block=walks_per_block,
                           max_steps=max_steps, checkpoint=checkpoint,
                           verbose=verbose)
        if verbose:
            print(f"  round {r + 1}: {m} merges, {state.n_blocks()} blocks "
                  f"({time.time() - t0:.1f}s)", file=sys.stderr)
        if state.n_blocks() == 1:
            break
        stalled = stalled + 1 if m == 0 else 0
        if stalled >= 8 and verbose:
            print(f"  no merges in {stalled} rounds; the partition looks "
                  f"stable at {state.n_blocks()} blocks", file=sys.stderr)
    return state
