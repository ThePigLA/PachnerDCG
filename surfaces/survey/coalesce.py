"""
Flip-graph coalescence checker for triangulated closed surfaces.

Question addressed: given the components of F_M(n), do their stabilizations
(1--3 facet subdivisions) lie in a single component of F_M(n+1)?

Key algorithmic points vs. the naive approach:

  * canonical_form uses a FLAG-BASED canonical labeling (6*|F| candidate
    labelings) instead of brute-forcing all n! vertex permutations.
    For n=10, |F|=24 this is 144 labelings instead of 3,628,800.

  * By the lift-fiber lemma, all facet subdivisions of a triangulation lie
    in one component upstairs, so ONE subdivision of ONE representative per
    downstairs component suffices as a seed.

  * Coalescence is checked by multi-source BFS + union-find with early
    termination, so the full level-(n+1) flip graph is never built unless
    the answer turns out to be negative.
"""

from __future__ import annotations

import ast
import re
import sys
from collections import deque
from itertools import combinations

Triangulation = tuple  # tuple of sorted 3-tuples, sorted


# --------------------------------------------------------------------------
# basic combinatorics
# --------------------------------------------------------------------------

def normalize(triangles) -> Triangulation:
    """Sorted tuple of sorted vertex triples."""
    return tuple(sorted(tuple(sorted(t)) for t in triangles))


def vertices_of(tri: Triangulation):
    return sorted({v for t in tri for v in t})


def edge_to_triangles(tri: Triangulation) -> dict:
    """Map each edge (sorted pair) to the list of incident triangles."""
    e2t = {}
    for t in tri:
        a, b, c = t
        for e in ((a, b), (a, c), (b, c)):
            e2t.setdefault(e, []).append(t)
    return e2t


def f_vector(tri: Triangulation):
    e2t = edge_to_triangles(tri)
    return (len(vertices_of(tri)), len(e2t), len(tri))


def euler_characteristic(tri: Triangulation) -> int:
    v, e, f = f_vector(tri)
    return v - e + f


def is_closed_surface(tri: Triangulation) -> bool:
    """Every edge in exactly two triangles and every vertex link a cycle."""
    e2t = edge_to_triangles(tri)
    if any(len(ts) != 2 for ts in e2t.values()):
        return False
    # vertex links must be single cycles
    link_edges = {}
    for t in tri:
        a, b, c = t
        link_edges.setdefault(a, []).append((b, c))
        link_edges.setdefault(b, []).append((a, c))
        link_edges.setdefault(c, []).append((a, b))
    for v, edges in link_edges.items():
        deg = {}
        for x, y in edges:
            deg[x] = deg.get(x, 0) + 1
            deg[y] = deg.get(y, 0) + 1
        if any(d != 2 for d in deg.values()):
            return False
        # connectivity of the link
        adj = {}
        for x, y in edges:
            adj.setdefault(x, []).append(y)
            adj.setdefault(y, []).append(x)
        start = next(iter(adj))
        seen, stack = {start}, [start]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        if len(seen) != len(adj):
            return False
    return True


def is_orientable(tri: Triangulation) -> bool:
    """Try to coherently orient the triangles; fail => nonorientable."""
    e2t = edge_to_triangles(tri)
    orient = {}          # triangle -> ordered triple
    start = tri[0]
    orient[start] = (start[0], start[1], start[2])
    queue = deque([start])
    while queue:
        t = queue.popleft()
        x, y, z = orient[t]
        for (p, q) in ((x, y), (y, z), (z, x)):
            e = (p, q) if p < q else (q, p)
            for nb in e2t[e]:
                if nb == t:
                    continue
                w = next(v for v in nb if v not in (p, q))
                # coherent orientation reverses the shared edge
                want = (q, p, w)
                if nb in orient:
                    cur = orient[nb]
                    rots = {cur, (cur[1], cur[2], cur[0]), (cur[2], cur[0], cur[1])}
                    if want not in rots:
                        return False
                else:
                    orient[nb] = want
                    queue.append(nb)
    return True


# --------------------------------------------------------------------------
# canonical form  (flag-based; exact for connected closed surfaces)
# --------------------------------------------------------------------------

_CANON_CACHE: dict = {}
# The cache is bounded: at scale, each flip result is a distinct raw complex,
# so the hit rate approaches zero while an unbounded cache would store one
# entry per flip attempted (tens of millions of entries on a large search).
_CANON_CACHE_MAX = 20_000


def _vertex_colors(tri, e2t, rounds=2):
    """Isomorphism-invariant vertex colors: degree, then refined by neighbors."""
    adj = {}
    for (u, v) in e2t:
        adj.setdefault(u, set()).add(v)
        adj.setdefault(v, set()).add(u)
    color = {v: len(adj[v]) for v in adj}
    for _ in range(rounds):
        sig = {v: (color[v], tuple(sorted(color[u] for u in adj[v]))) for v in adj}
        order = {s: i for i, s in enumerate(sorted(set(sig.values())))}
        color = {v: order[sig[v]] for v in adj}
    return color


def canonical_form(tri) -> Triangulation:
    """
    Canonical representative of the isomorphism class of a connected closed
    surface triangulation.

    Method.  A *flag* is a triangle together with an ordering of its three
    vertices.  Each flag determines a deterministic BFS over the
    facet-adjacency graph; labelling vertices in order of discovery turns the
    complex into a labelled one, and recording the triangles in BFS order
    gives a string.  The lexicographic minimum of that string over all flags
    depends only on the isomorphism class, so it is a canonical form.  This is
    exact because a connected closed surface triangulation is reconstructible
    from a single rooted traversal.

    Cost: O(|F|^2) worst case (6|F| flags, O(|F|) work each), versus O(n! |F|)
    for brute-force permutation search.  Two further prunings are applied:

      * only flags whose vertex-color triple is lexicographically minimal are
        considered (sound, because the coloring is isomorphism-invariant);
      * a traversal is abandoned as soon as its prefix exceeds the best
        string found so far.
    """
    tri = normalize(tri)
    hit = _CANON_CACHE.get(tri)
    if hit is not None:
        return hit

    e2t = edge_to_triangles(tri)
    nbr = {}
    for e, ts in e2t.items():
        if len(ts) != 2:
            raise ValueError("not a closed surface: edge %r in %d triangles"
                             % (e, len(ts)))
        nbr[(ts[0], e)] = ts[1]
        nbr[(ts[1], e)] = ts[0]

    color = _vertex_colors(tri, e2t)

    # invariant restriction of the flag set
    flags = []
    for t0 in tri:
        a, b, c = t0
        for start in ((a, b, c), (a, c, b), (b, a, c),
                      (b, c, a), (c, a, b), (c, b, a)):
            flags.append(((color[start[0]], color[start[1]], color[start[2]]),
                          t0, start))
    best_key = min(f[0] for f in flags)
    flags = [(t0, start) for key, t0, start in flags if key == best_key]

    nfac = len(tri)
    best = None
    for t0, start in flags:
        lab = {start[0]: 0, start[1]: 1, start[2]: 2}
        nxt = 3
        out = []
        seen = {t0}
        queue = deque([start])
        while queue:
            x, y, z = queue.popleft()
            key = tuple(sorted((x, y, z)))
            out.append((lab[x], lab[y], lab[z]))
            if best is not None and out[-1] != best[len(out) - 1]:
                if out[-1] > best[len(out) - 1]:
                    out = None
                    break
                best = None           # strictly better prefix; accept this one
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
        if out is None:
            continue
        cand = tuple(out)
        if len(cand) != nfac:
            raise ValueError("triangulation is not connected")
        if best is None or cand < best:
            best = cand

    result = tuple(sorted(tuple(sorted(t)) for t in best))
    if len(_CANON_CACHE) >= _CANON_CACHE_MAX:
        _CANON_CACHE.clear()
    _CANON_CACHE[tri] = result
    return result


# --------------------------------------------------------------------------
# moves
# --------------------------------------------------------------------------

def diagonal_flips(tri: Triangulation):
    """All legal simplicial 2--2 flips, returned as normalized triangulations."""
    e2t = edge_to_triangles(tri)
    existing = set(e2t)
    out = []
    for e, ts in e2t.items():
        if len(ts) != 2:
            continue
        t1, t2 = ts
        v1 = next(v for v in t1 if v not in e)
        v2 = next(v for v in t2 if v not in e)
        new_e = (v1, v2) if v1 < v2 else (v2, v1)
        if new_e in existing:
            continue                      # simplicial condition fails
        drop = {t1, t2}
        rest = [t for t in tri if t not in drop]
        rest.append(tuple(sorted((e[0], v1, v2))))
        rest.append(tuple(sorted((e[1], v1, v2))))
        out.append(normalize(rest))
    return out


def subdivide_facet(tri: Triangulation, facet, new_vertex):
    """1--3 stellar subdivision of one facet."""
    rest = [t for t in tri if t != facet]
    a, b, c = facet
    rest += [tuple(sorted((a, b, new_vertex))),
             tuple(sorted((a, c, new_vertex))),
             tuple(sorted((b, c, new_vertex)))]
    return normalize(rest)


# --------------------------------------------------------------------------
# component machinery
# --------------------------------------------------------------------------

class DSU:
    def __init__(self, items):
        self.p = {x: x for x in items}

    def find(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx != ry:
            self.p[ry] = rx
            return True
        return False

    def add(self, x):
        self.p.setdefault(x, x)


def components_at_level(seeds, verbose=False):
    """
    Full component decomposition of the flip graph generated by `seeds`.
    Returns a list of sets of canonical forms.
    """
    seeds = [canonical_form(s) for s in seeds]
    dsu = DSU(seeds)
    seen = set(seeds)
    queue = deque(seeds)
    reps = {}
    for s in seeds:
        reps[s] = s
    while queue:
        cur = queue.popleft()
        for nb_raw in diagonal_flips(cur):
            nb = canonical_form(nb_raw)
            dsu.add(nb)
            dsu.union(cur, nb)
            if nb not in seen:
                seen.add(nb)
                queue.append(nb)
        if verbose and len(seen) % 500 == 0:
            print(f"    ... {len(seen)} triangulations seen", file=sys.stderr)
    buckets = {}
    for x in seen:
        buckets.setdefault(dsu.find(x), set()).add(x)
    return list(buckets.values())


def do_seeds_coalesce(seeds, node_budget=2_000_000, verbose=True):
    """
    Multi-source BFS with union-find and early termination.

    Returns (coalesced: bool, nodes_explored: int, dsu, seen).
    Stops as soon as all seeds are in one class.
    """
    seeds = [canonical_form(s) for s in seeds]
    if len(set(seeds)) == 1:
        return True, 0, None, set(seeds)

    dsu = DSU(seeds)
    seen = set(seeds)
    queue = deque(seeds)
    explored = 0

    while queue:
        cur = queue.popleft()
        explored += 1
        for nb_raw in diagonal_flips(cur):
            nb = canonical_form(nb_raw)
            dsu.add(nb)
            dsu.union(cur, nb)
            if nb not in seen:
                seen.add(nb)
                queue.append(nb)
        roots = {dsu.find(s) for s in seeds}
        if len(roots) == 1:
            return True, explored, dsu, seen
        if verbose and explored % 500 == 0:
            print(f"    ... explored {explored}, seen {len(seen)}, "
                  f"{len(roots)} classes remaining", file=sys.stderr)
        if explored > node_budget:
            raise RuntimeError("node budget exceeded")

    return False, explored, dsu, seen


# --------------------------------------------------------------------------
# parsing Lutz-style files
# --------------------------------------------------------------------------

def parse_triangulations_file(path):
    """
    Extract every [[...],[...],...] block from a Lutz-style manifold file.
    Tolerant of the '#12=' prefixes, line breaks and stray whitespace.
    """
    with open(path) as fh:
        content = fh.read()
    out = []
    for m in re.finditer(r'\[\s*\[.*?\]\s*\]', content, re.DOTALL):
        block = re.sub(r'\s+', '', m.group(0))
        try:
            data = ast.literal_eval(block)
        except (ValueError, SyntaxError):
            continue
        if not data or not all(isinstance(t, (list, tuple)) and len(t) == 3 for t in data):
            continue
        out.append(normalize(data))
    return out


def select_surface(triangulations, chi, orientable):
    """Keep only closed surfaces with the given chi and orientability."""
    keep = []
    for t in triangulations:
        if not is_closed_surface(t):
            continue
        if euler_characteristic(t) != chi:
            continue
        if is_orientable(t) != orientable:
            continue
        keep.append(t)
    return keep


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def main(path, chi=-2, orientable=False):
    print(f"Parsing {path} ...")
    raw = parse_triangulations_file(path)
    print(f"  extracted {len(raw)} triangulation blocks")

    target = select_surface(raw, chi=chi, orientable=orientable)
    print(f"  {len(target)} are closed surfaces with chi={chi}, "
          f"orientable={orientable}")
    if not target:
        print("  nothing to do")
        return

    print("Canonicalizing ...")
    level_n = {canonical_form(t) for t in target}
    n = len(vertices_of(next(iter(level_n))))
    print(f"  {len(level_n)} isomorphism classes at n={n}")

    print(f"Computing components of F(n={n}) ...")
    comps = components_at_level(level_n)
    comps.sort(key=len, reverse=True)
    print(f"  {len(comps)} components, sizes {[len(c) for c in comps]}")

    # one subdivision of one representative per component is enough
    print(f"Stabilizing one representative per component to n={n+1} ...")
    seeds = []
    for c in comps:
        rep = min(c)
        seeds.append(subdivide_facet(rep, rep[0], new_vertex=n))
    print(f"  {len(seeds)} seeds")

    print(f"Searching F(n={n+1}) for coalescence ...")
    ok, explored, dsu, seen = do_seeds_coalesce(seeds)
    print(f"  explored {explored} triangulations, {len(seen)} discovered")
    if ok:
        print(f"\nCOALESCED: all {len(comps)} components merge at n={n+1}. "
              f"Absorption time = 1.")
    else:
        roots = {}
        for s in [canonical_form(x) for x in seeds]:
            roots.setdefault(dsu.find(s), []).append(s)
        print(f"\nSTILL BRANCHED: the {len(comps)} rays land in "
              f"{len(roots)} distinct components at n={n+1}.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        print("usage: python coalesce.py <lutz_file> [chi] [orientable:0|1]")
        sys.exit(1)
    chi = int(sys.argv[2]) if len(sys.argv) > 2 else -2
    ori = bool(int(sys.argv[3])) if len(sys.argv) > 3 else False
    main(sys.argv[1], chi=chi, orientable=ori)
