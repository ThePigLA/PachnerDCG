"""
Fibers of the stabilization map s_n : Pi_M(n) -> Pi_M(n+1), for closed surfaces.

The idea
--------
If U is an (n+1)-vertex triangulation with a vertex v of degree 3, then the
link of v is a triangle, so U is the stellar subdivision S^sigma of the
n-vertex triangulation S obtained by the 3--1 move at v.  Hence U lies in the
lift fiber L(S), and by the lift-fiber lemma

        [U]  =  s_n([S])   in  Pi_M(n+1).

So every degree-3 vertex met while walking in F_M(n+1) projects down to a
*known* component downstairs.  If a walk started from a subdivision of S_i
reaches some U that projects into class j, then

        s_n(C_i) = [U] = s_n(C_j),

and we may union i and j.  Two classes merge under stabilization exactly when
they are joined by a chain of such observations.

Why this is better than exhaustive BFS upstairs
-----------------------------------------------
Exhaustive multi-source BFS has to make two frontiers physically collide in
F_M(n+1), which for a high-genus surface means storing millions of
triangulations.  Degree-3 projection instead reduces every observation to a
lookup in a table with |Pi_M(n)| entries, so the search needs only
O(walk length) memory and each walk is independent.

Completeness
------------
Positive merges are cheap and certified: a merge is witnessed by an explicit
flip path plus one 3--1 move.  Proving that two classes do NOT merge still
requires exhausting their components upstairs, because the down-set of a
component upstairs is precisely its fiber.  `exhaustive_fibers` does that with
hash-only bookkeeping when it is feasible; `walk_fibers` is the tool for the
positive direction.
"""

from __future__ import annotations

import random
import sys
from collections import deque
from hashlib import blake2b

from coalesce import (normalize, canonical_form, edge_to_triangles,
                      diagonal_flips, subdivide_facet, vertices_of,
                      f_vector, euler_characteristic)


# --------------------------------------------------------------------------
# compact representation
# --------------------------------------------------------------------------

def pack(tri) -> bytes:
    """Sorted triangulation -> bytes (3 bytes per triangle).  Needs n <= 256."""
    out = bytearray()
    for a, b, c in normalize(tri):
        out += bytes((a, b, c))
    return bytes(out)


def unpack(b: bytes):
    return tuple((b[i], b[i + 1], b[i + 2]) for i in range(0, len(b), 3))


def digest(tri) -> int:
    """128-bit fingerprint of a canonical form.  Collisions are negligible."""
    return int.from_bytes(blake2b(pack(tri), digest_size=16).digest(), "big")


# --------------------------------------------------------------------------
# degrees and the 3--1 move
# --------------------------------------------------------------------------

def degrees(tri) -> dict:
    deg = {}
    for e in edge_to_triangles(tri):
        deg[e[0]] = deg.get(e[0], 0) + 1
        deg[e[1]] = deg.get(e[1], 0) + 1
    return deg


def degree3_vertices(tri):
    return [v for v, d in degrees(tri).items() if d == 3]


def contract_degree3(tri, v):
    """
    3--1 move at a degree-3 vertex v.  Returns the smaller triangulation, or
    None if the move is not simplicial (the link triangle already present).
    """
    incident = [t for t in tri if v in t]
    if len(incident) != 3:
        return None
    link = sorted({u for t in incident for u in t if u != v})
    if len(link) != 3:
        return None
    new_face = tuple(link)
    rest = [t for t in tri if v not in t]
    if new_face in set(rest):
        return None
    rest.append(new_face)
    return normalize(rest)


# --------------------------------------------------------------------------
# union-find over the downstairs classes
# --------------------------------------------------------------------------

class DSU:
    def __init__(self, n):
        self.p = list(range(n))
        self.n_classes = n

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return False
        self.p[ry] = rx
        self.n_classes -= 1
        return True

    def classes(self):
        buckets = {}
        for x in range(len(self.p)):
            buckets.setdefault(self.find(x), []).append(x)
        return list(buckets.values())


# --------------------------------------------------------------------------
# the base index: canonical n-vertex triangulation -> component id
# --------------------------------------------------------------------------

def build_base_index(components):
    """
    components: list of sets of canonical n-vertex triangulations
    (one set per component of F_M(n)).  Returns dict canon -> component id.
    """
    index = {}
    for i, comp in enumerate(components):
        for t in comp:
            index[canonical_form(t)] = i
    return index


def project(tri, base_index):
    """
    Every degree-3 vertex of `tri` gives a component id downstairs.
    Returns the set of ids found.
    """
    found = set()
    for v in degree3_vertices(tri):
        small = contract_degree3(tri, v)
        if small is None:
            continue
        cid = base_index.get(canonical_form(small))
        if cid is not None:
            found.add(cid)
    return found


# --------------------------------------------------------------------------
# degree-reduction walk: try to manufacture a degree-3 vertex
# --------------------------------------------------------------------------

def _flip_on_edge(tri, e, e2t=None):
    """Perform the 2--2 flip on edge e if legal; else None."""
    e2t = e2t or edge_to_triangles(tri)
    ts = e2t.get(e)
    if ts is None or len(ts) != 2:
        return None
    t1, t2 = ts
    v1 = next(v for v in t1 if v not in e)
    v2 = next(v for v in t2 if v not in e)
    new_e = (v1, v2) if v1 < v2 else (v2, v1)
    if new_e in e2t:
        return None
    drop = {t1, t2}
    rest = [t for t in tri if t not in drop]
    rest.append(tuple(sorted((e[0], v1, v2))))
    rest.append(tuple(sorted((e[1], v1, v2))))
    return normalize(rest)


def degree_reduction_walk(start, base_index, rng, max_steps=400,
                          restart_temp=0.15):
    """
    Randomized walk that pushes some vertex's degree down toward 3, reporting
    every downstairs class id it can project to along the way.

    A 2--2 flip on edge vw lowers deg(v) and deg(w) by one, so repeatedly
    flipping edges at a low-degree vertex drives its degree to 3.
    """
    cur = normalize(start)
    hits = set()
    hits |= project(cur, base_index)
    path = [cur]

    for _ in range(max_steps):
        deg = degrees(cur)
        e2t = edge_to_triangles(cur)
        # focus on a low-degree vertex, with some randomness
        order = sorted(deg, key=lambda v: (deg[v], rng.random()))
        target = order[0] if rng.random() > restart_temp else rng.choice(order[:4])

        cands = [e for e in e2t if target in e]
        rng.shuffle(cands)
        moved = None
        for e in cands:
            nxt = _flip_on_edge(cur, e, e2t)
            if nxt is not None:
                moved = nxt
                break
        if moved is None:
            # target is stuck; take any legal flip
            allf = diagonal_flips(cur)
            if not allf:
                break
            moved = rng.choice(allf)

        cur = moved
        path.append(cur)
        new = project(cur, base_index)
        if new - hits:
            hits |= new
    return hits, path


def walk_fibers(seeds, base_index, n_base_classes, trials=40, max_steps=400,
                seed=0, verbose=True):
    """
    Positive-direction fiber computation.

    seeds: list where seeds[i] is one facet subdivision of one representative
           of component i downstairs (by the lift-fiber lemma, any choice will
           do).
    Returns (DSU, witnesses) where witnesses[(i, j)] is a flip path proving
    that classes i and j stabilize to the same component.
    """
    rng = random.Random(seed)
    dsu = DSU(n_base_classes)
    witnesses = {}

    for rounds in range(trials):
        for i, s in enumerate(seeds):
            if dsu.n_classes == 1:
                break
            hits, path = degree_reduction_walk(s, base_index, rng,
                                               max_steps=max_steps)
            for j in hits:
                if dsu.union(i, j):
                    witnesses[(min(i, j), max(i, j))] = path
        if verbose:
            print(f"  round {rounds + 1}: {dsu.n_classes} classes remaining",
                  file=sys.stderr)
        if dsu.n_classes == 1:
            break
    return dsu, witnesses


# --------------------------------------------------------------------------
# exhaustive mode, memory-lean
# --------------------------------------------------------------------------

def exhaustive_fibers(seeds, base_index, n_base_classes, node_budget=5_000_000,
                      verbose=True):
    """
    Exhaust the components upstairs, storing only 128-bit fingerprints in the
    visited set and packed bytes in the frontier.  Expands the class with the
    smallest frontier first, which merges classes far sooner than round-robin
    BFS.

    Only this mode can prove a NEGATIVE (that two classes really do not
    merge), since the down-set of a component upstairs is exactly its fiber.
    """
    dsu = DSU(n_base_classes)
    seen = set()                       # fingerprints, all classes pooled
    fronts = {}                        # class root -> deque of packed forms
    owner = {}                         # fingerprint -> class root

    for i, s in enumerate(seeds):
        c = canonical_form(s)
        h = digest(c)
        seen.add(h)
        owner[h] = i
        fronts[i] = deque([pack(c)])
        for j in project(c, base_index):
            dsu.union(i, j)

    explored = 0
    while any(fronts.values()):
        # smallest non-empty frontier first
        root = min((r for r, q in fronts.items() if q), key=lambda r: len(fronts[r]))
        root = dsu.find(root) if dsu.find(root) in fronts else root
        q = fronts[root]
        cur = unpack(q.popleft())
        explored += 1

        for nb_raw in diagonal_flips(cur):
            nb = canonical_form(nb_raw)
            h = digest(nb)
            if h in seen:
                other = owner.get(h)
                if other is not None:
                    dsu.union(root, other)
                continue
            seen.add(h)
            owner[h] = root
            for j in project(nb, base_index):
                dsu.union(root, j)
            q.append(pack(nb))

        if dsu.n_classes == 1:
            return dsu, explored, len(seen)
        if verbose and explored % 5000 == 0:
            print(f"  explored {explored}, seen {len(seen)}, "
                  f"{dsu.n_classes} classes, frontier {sum(map(len, fronts.values()))}",
                  file=sys.stderr)
        if explored > node_budget:
            raise RuntimeError(f"node budget exceeded; {dsu.n_classes} classes remain")

    return dsu, explored, len(seen)
