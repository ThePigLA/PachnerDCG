"""
Initial component trees for every closed surface occurring at small vertex
counts.

For each surface M and each level n for which complete data is available,
this records:

  * the components of F_M(n): how many, and their sizes;
  * the fibers of the stabilization map s_n : Pi_M(n) -> Pi_M(n+1);
  * the components at level n+1 that are NOT in the image of s_n
    ("newly born" branches, in the language of the component tree);
  * the largest one-step fiber;
  * the coalescence distance: over each fiber, the largest shortest flip
    distance in F_M(n+1) between the stabilized representatives of its
    members, maximised over fibers.

DATA REQUIREMENT
----------------
Component counts and fiber sizes need only the triangulations at level n.
Newly born branches need the COMPLETE enumeration at level n+1, because a
newborn component is by definition unreachable from below.  So this module
consumes complete level files (Lutz's manifold_lex_d2_n*_o*), not flip
closures of seeds.

WELL-DEFINEDNESS
----------------
s_n([T]) is computed from one facet subdivision of one representative.  By
the lift-fiber lemma every other choice lands in the same component, and
`verify_stabilization_welldefined` checks that empirically across all
components and all facets -- an independent test of the lemma on real data.
"""

from __future__ import annotations

import json
import sys
import time
from collections import deque, defaultdict

from coalesce import (normalize, canonical_form, diagonal_flips,
                      subdivide_facet, vertices_of, f_vector,
                      euler_characteristic, is_closed_surface, is_orientable,
                      parse_triangulations_file, edge_to_triangles)


# --------------------------------------------------------------------------
# small local helpers
#
# These are deliberately defined here rather than imported, so that the survey
# pipeline is exactly coalesce.py + survey.py + batch.py with no further
# dependencies.
# --------------------------------------------------------------------------

def degrees(tri) -> dict:
    """Vertex degrees (number of incident edges)."""
    deg = {}
    for e in edge_to_triangles(tri):
        deg[e[0]] = deg.get(e[0], 0) + 1
        deg[e[1]] = deg.get(e[1], 0) + 1
    return deg


def contract_degree3(tri, v):
    """
    The 3--1 move at a degree-3 vertex v: replace the three triangles at v by
    the triangle spanned by its link.  Returns None when that is not
    simplicial (the link triangle is already a face, which for a closed
    surface happens only for the 4-vertex sphere).
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
# surface naming
# --------------------------------------------------------------------------

def surface_name(chi, orientable):
    if orientable:
        g = (2 - chi) // 2
        return "S^2" if g == 0 else ("T^2" if g == 1 else f"Sigma_{g}")
    k = 2 - chi
    return {1: "RP^2", 2: "Klein"}.get(k, f"N_{k}")


def min_vertices(chi, orientable):
    """
    Minimal vertex number of a triangulation, from the Heawood/Ringel bound

        n >= ceil( (7 + sqrt(49 - 24 chi)) / 2 ).

    Jungerman and Ringel showed the bound is attained for every closed surface
    EXCEPT three: the orientable surface of genus 2, the Klein bottle, and N_3,
    each of which needs one vertex more than the bound.
    """
    import math
    n = math.ceil((7 + math.sqrt(49 - 24 * chi)) / 2)
    exceptions = {
        (-2, True): 10,     # Sigma_2: bound gives 9
        (0, False): 8,      # Klein bottle: bound gives 7
        (-1, False): 9,     # N_3: bound gives 8
    }
    return max(n, exceptions.get((chi, orientable), 0))


# --------------------------------------------------------------------------
# loading complete levels
# --------------------------------------------------------------------------

def load_levels(paths, verbose=True):
    """
    paths: iterable of Lutz-format files, any mix of vertex counts and
    orientability.  Returns {(chi, orientable): {n: sorted list of canonical
    triangulations}}.
    """
    data = defaultdict(lambda: defaultdict(set))
    for p in paths:
        raw = parse_triangulations_file(p)
        kept = 0
        for t in raw:
            if not is_closed_surface(t):
                continue
            chi = euler_characteristic(t)
            ori = is_orientable(t)
            n = len(vertices_of(t))
            data[(chi, ori)][n].add(canonical_form(t))
            kept += 1
        if verbose:
            print(f"  {p}: {len(raw)} blocks, {kept} closed surfaces",
                  file=sys.stderr)
    return {k: {n: sorted(v) for n, v in sorted(lv.items())}
            for k, lv in data.items()}


# --------------------------------------------------------------------------
# components of a complete level
# --------------------------------------------------------------------------

def components_of_level(level, progress_every=5000, label=""):
    """
    level: complete collection of canonical triangulations with n vertices.
    Returns (components, index) where components is a list of sorted lists and
    index maps each canonical form to its component number.

    Because the level is complete, no exploration is needed: flips stay inside
    the level, so this is plain connected-components on the induced subgraph.
    """
    level = list(level)
    present = set(level)
    index = {}
    comps = []
    done = 0
    t_start = time.time()
    for start in level:
        if start in index:
            continue
        cid = len(comps)
        comp = []
        queue = deque([start])
        index[start] = cid
        while queue:
            cur = queue.popleft()
            comp.append(cur)
            done += 1
            if progress_every and done % progress_every == 0:
                rate = done / max(time.time() - t_start, 1e-9)
                eta = (len(level) - done) / max(rate, 1e-9)
                print(f"    {label}{done}/{len(level)} "
                      f"(in component #{len(comps) + 1}, "
                      f"{rate:.0f}/s, ~{eta / 60:.1f} min left)",
                      file=sys.stderr)
            for nb_raw in diagonal_flips(cur):
                nb = canonical_form(nb_raw)
                if nb not in present:
                    raise ValueError(
                        "flip left the level: the input enumeration is "
                        "incomplete or inconsistent")
                if nb not in index:
                    index[nb] = cid
                    queue.append(nb)
        comps.append(sorted(comp))
    return comps, index


# --------------------------------------------------------------------------
# stabilization
# --------------------------------------------------------------------------

def stabilize(tri, n, facet=None):
    return canonical_form(subdivide_facet(tri, facet or tri[0], n))


def stabilization_map(comps_n, n, index_np1):
    """
    s_n on components.  Returns list `img` with img[i] = component id at
    level n+1 containing the stabilizations of component i.
    """
    img = []
    for comp in comps_n:
        rep = comp[0]
        up = stabilize(rep, n)
        if up not in index_np1:
            raise ValueError("stabilized triangulation missing from level "
                             f"{n+1}: enumeration incomplete")
        img.append(index_np1[up])
    return img


def verify_stabilization_welldefined(comps_n, n, index_np1, max_checks=None):
    """
    Empirical check of the lift-fiber lemma: every facet subdivision of every
    member of a component must land in the same component upstairs.
    Returns (checks_done, violations).
    """
    checks = 0
    bad = 0
    for comp in comps_n:
        target = None
        for t in comp:
            for facet in t:
                up = canonical_form(subdivide_facet(t, facet, n))
                cid = index_np1.get(up)
                if target is None:
                    target = cid
                elif cid != target:
                    bad += 1
                checks += 1
                if max_checks and checks >= max_checks:
                    return checks, bad
    return checks, bad


# --------------------------------------------------------------------------
# newborn detection without the level below
# --------------------------------------------------------------------------

def is_subdivision(tri):
    """
    True iff `tri` is a stellar subdivision of a smaller triangulation, i.e.
    it has a vertex of degree 3 whose 3--1 move is simplicial.  Returns the
    witnessing vertex, or None.

    For n >= 5 the link triangle of a degree-3 vertex cannot already be a
    face, so this is just "has a vertex of degree 3"; the contraction is
    attempted anyway so the n = 4 case is handled correctly too.
    """
    for v, d in degrees(tri).items():
        if d == 3 and contract_degree3(tri, v) is not None:
            return v
    return None


def component_is_newborn(comp):
    """
    A component C of F_M(n) is in the image of s_{n-1} iff some member is a
    stellar subdivision:

      (<=)  if U in C has a degree-3 vertex v, then U = S^sigma for the S
            obtained by the 3--1 move, so U lies in the lift fiber L(S) and
            C = s_{n-1}([S]).
      (=>)  if C = s_{n-1}(D) then C contains T^sigma for some T in D, and
            the subdividing vertex has degree 3.

    So C is NEWBORN iff no member has such a vertex -- a test that needs only
    level n, never the level below.  Returns (newborn, witness, min_degree).
    """
    min_deg = None
    for t in comp:
        d = degrees(t)
        lo = min(d.values())
        min_deg = lo if min_deg is None else min(min_deg, lo)
        v = is_subdivision(t)
        if v is not None:
            return False, (t, v), min_deg
    return True, None, min_deg


def newborn_analysis(comps, progress_every=200, label=""):
    """
    Classify every component of a level as newborn or not.  Short-circuits on
    the first subdivision found, so large components are usually settled
    immediately.
    """
    out = []
    for i, comp in enumerate(comps):
        if progress_every and i and i % progress_every == 0:
            print(f"    {label}newborn scan {i}/{len(comps)}",
                  file=sys.stderr)
        newborn, witness, min_deg = component_is_newborn(comp)
        out.append({
            "id": i,
            "size": len(comp),
            "newborn": newborn,
            "min_degree": min_deg,
        })
    return out


# --------------------------------------------------------------------------
# flip distances
# --------------------------------------------------------------------------

def bfs_distances(source, present, targets=None, progress=None):
    """
    Shortest flip distances from `source` within the level `present`.

    If `targets` is given, the search stops as soon as every target has been
    labelled -- for a fiber of k members we only ever need k-1 distances, and
    exploring the whole level to find them is pure waste.
    """
    dist = {source: 0}
    queue = deque([source])
    remaining = set(targets) - {source} if targets else None
    seen = 0
    while queue:
        cur = queue.popleft()
        seen += 1
        if progress and seen % 5000 == 0:
            progress(seen, len(dist), len(remaining) if remaining else 0)
        d = dist[cur] + 1
        for nb_raw in diagonal_flips(cur):
            nb = canonical_form(nb_raw)
            if nb in present and nb not in dist:
                dist[nb] = d
                queue.append(nb)
                if remaining is not None and nb in remaining:
                    remaining.discard(nb)
                    if not remaining:
                        return dist
    return dist


def coalescence_distance(fiber, comps_n, n, level_np1_set):
    """
    Largest shortest flip distance in F_M(n+1) between the stabilized
    representatives of the components in one fiber.  Returns (max_dist, pair)
    or (0, None) for a singleton fiber.
    """
    if len(fiber) < 2:
        return 0, None
    reps = [stabilize(comps_n[i][0], n) for i in fiber]
    best, pair = 0, None
    for a in range(len(reps)):
        need = set(reps[a + 1:])
        if not need:
            break
        prog = None
        if len(level_np1_set) > 20000:
            def prog(seen, labelled, left, _a=a, _t=len(reps) - 1):
                print(f"      distance BFS {_a + 1}/{_t}: {seen} expanded, "
                      f"{labelled} labelled, {left} target(s) left",
                      file=sys.stderr)
        dist = bfs_distances(reps[a], level_np1_set, targets=need,
                             progress=prog)
        for b in range(a + 1, len(reps)):
            d = dist.get(reps[b])
            if d is None:
                raise ValueError("fiber members are not connected upstairs; "
                                 "the fiber computation is inconsistent")
            if d > best:
                best, pair = d, (fiber[a], fiber[b])
    return best, pair


# --------------------------------------------------------------------------
# the survey
# --------------------------------------------------------------------------

def survey_surface(chi, orientable, levels, verbose=True, verify=True):
    """
    levels: {n: complete sorted list of canonical triangulations}
    Returns a list of per-level records.
    """
    name = surface_name(chi, orientable)
    ns = sorted(levels)
    cache = {}
    for n in ns:
        t0 = time.time()
        comps, index = components_of_level(levels[n],
                                           label=f"{name} n={n}: ")
        comps.sort(key=len, reverse=True)
        index = {t: i for i, c in enumerate(comps) for t in c}
        cache[n] = (comps, index)
        if verbose:
            print(f"  {name} n={n}: {len(levels[n])} triangulations, "
                  f"{len(comps)} components ({time.time()-t0:.1f}s)",
                  file=sys.stderr)

    records = []
    for n in ns:
        comps, index = cache[n]
        na = newborn_analysis(comps, label=f"{name} n={n}: ")
        nb_ids = [a["id"] for a in na if a["newborn"]]
        floor = min_vertices(chi, orientable)

        rec = {
            "surface": name,
            "chi": chi,
            "orientable": orientable,
            "n": n,
            "n_triangulations": len(levels[n]),
            "n_components": len(comps),
            "component_sizes": [len(c) for c in comps],
            # newborn detection from level n alone, via the degree-3 criterion
            "at_minimum_n": n == floor,
            "min_vertices": floor,
            "newborn_here": nb_ids,
            "n_newborn_here": len(nb_ids),
            "newborn_here_sizes": [len(comps[i]) for i in nb_ids],
            "min_degree_by_component": [a["min_degree"] for a in na],
        }

        # At the minimal vertex number no triangulation can have a degree-3
        # vertex (the 3--1 move would produce a smaller one), so EVERY
        # component must come out newborn.  Any other answer is a bug.
        if n == floor and len(nb_ids) != len(comps):
            print(f"  !! {name} n={n}: at minimum vertex count but only "
                  f"{len(nb_ids)}/{len(comps)} components are newborn",
                  file=sys.stderr)
            rec["consistency_warning"] = "minimal level with non-newborn component"

        if n + 1 in cache:
            comps1, index1 = cache[n + 1]
            level1 = set(levels[n + 1])
            img = stabilization_map(comps, n, index1)

            fibers = defaultdict(list)
            for i, target in enumerate(img):
                fibers[target].append(i)
            fibers = {k: sorted(v) for k, v in fibers.items()}

            newborn = [j for j in range(len(comps1)) if j not in fibers]

            fiber_sizes = sorted((len(v) for v in fibers.values()), reverse=True)
            worst_d, worst_pair, worst_fiber = 0, None, None
            for target, fib in fibers.items():
                d, pair = coalescence_distance(fib, comps, n, level1)
                if d > worst_d:
                    worst_d, worst_pair, worst_fiber = d, pair, target

            rec.update({
                "n_components_next": len(comps1),
                "component_sizes_next": [len(c) for c in comps1],
                "image_size": len(fibers),
                "fiber_sizes": fiber_sizes,
                "largest_fiber": fiber_sizes[0] if fiber_sizes else 0,
                "fibers": {str(k): v for k, v in fibers.items()},
                "newborn_components": newborn,
                "n_newborn": len(newborn),
                "newborn_sizes": [len(comps1[j]) for j in newborn],
                "coalescence_distance": worst_d,
                "coalescence_pair": list(worst_pair) if worst_pair else None,
                "coalescence_fiber": worst_fiber,
                "absorption_one_step": len(fibers) == 1,
            })

            # Independent cross-check: the components at level n+1 missed by
            # s_n must be exactly those with no degree-3 vertex.  The two
            # computations share no code path.
            na1 = newborn_analysis(comps1, label=f"{name} n={n+1}: ")
            by_degree = sorted(a["id"] for a in na1 if a["newborn"])
            rec["newborn_by_degree_criterion"] = by_degree
            rec["newborn_criteria_agree"] = (by_degree == sorted(newborn))
            if by_degree != sorted(newborn):
                print(f"  !! {name} n={n}: newborn sets disagree -- "
                      f"fibers say {sorted(newborn)}, degree-3 says "
                      f"{by_degree}", file=sys.stderr)

            if verify:
                checks, bad = verify_stabilization_welldefined(comps, n, index1)
                rec["welldefined_checks"] = checks
                rec["welldefined_violations"] = bad
                if bad:
                    print(f"  !! {name} n={n}: {bad} well-definedness "
                          f"violations", file=sys.stderr)

        records.append(rec)
    return records


def run_survey(paths, verbose=True, verify=True):
    data = load_levels(paths, verbose=verbose)
    out = []
    for (chi, ori) in sorted(data, key=lambda k: (-k[0], not k[1])):
        out.extend(survey_surface(chi, ori, data[(chi, ori)],
                                  verbose=verbose, verify=verify))
    return out


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def print_table(records, stream=sys.stdout):
    w = stream.write
    hdr = (f"{'surface':>9} {'chi':>4} {'n':>3} {'#tri':>7} {'comps':>6} "
           f"{'sizes':>22} {'nbrn*':>6} {'minD':>5} {'->img':>6} "
           f"{'fibers':>14} {'newborn':>8} {'maxfib':>7} {'dist':>5}")
    w(hdr + "\n" + "-" * len(hdr) + "\n")
    for r in records:
        sizes = str(r["component_sizes"])
        if len(sizes) > 22:
            sizes = sizes[:19] + "..."
        nb = str(r.get("n_newborn_here", "-"))
        if r.get("at_minimum_n"):
            nb += "!"
        mind = r.get("min_degree_by_component") or []
        mind = str(min(x for x in mind if x is not None)) if mind else "-"
        if "image_size" in r:
            fs = str(r["fiber_sizes"])
            if len(fs) > 14:
                fs = fs[:11] + "..."
            w(f"{r['surface']:>9} {r['chi']:>4} {r['n']:>3} "
              f"{r['n_triangulations']:>7} {r['n_components']:>6} "
              f"{sizes:>22} {nb:>6} {mind:>5} {r['image_size']:>6} "
              f"{fs:>14} {r['n_newborn']:>8} {r['largest_fiber']:>7} "
              f"{r['coalescence_distance']:>5}\n")
        else:
            w(f"{r['surface']:>9} {r['chi']:>4} {r['n']:>3} "
              f"{r['n_triangulations']:>7} {r['n_components']:>6} "
              f"{sizes:>22} {nb:>6} {mind:>5} {'-':>6} {'-':>14} "
              f"{'-':>8} {'-':>7} {'-':>5}\n")


def print_detail(records, stream=sys.stdout):
    w = stream.write
    for r in records:
        w(f"\n{'=' * 66}\n{r['surface']}  (chi = {r['chi']}, "
          f"{'orientable' if r['orientable'] else 'nonorientable'})  "
          f"level n = {r['n']}\n{'=' * 66}\n")
        w(f"  triangulations: {r['n_triangulations']}\n")
        w(f"  components of F(n): {r['n_components']}  "
          f"sizes {r['component_sizes']}\n")
        floor = r.get("min_vertices")
        w(f"  minimum vertex number for this surface: {floor}\n")
        if r.get("at_minimum_n"):
            w(f"  n is MINIMAL: no triangulation can have a degree-3 vertex, "
              f"so all {r['n_newborn_here']} component(s) are newborn "
              f"vacuously\n")
        else:
            w(f"  newborn at this level (no member is a subdivision): "
              f"{r['n_newborn_here']} of {r['n_components']}")
            if r["n_newborn_here"]:
                w(f"  ids {r['newborn_here']}, sizes "
                  f"{r['newborn_here_sizes']}")
            w("\n")
        md = [x for x in (r.get("min_degree_by_component") or [])
              if x is not None]
        if md:
            w(f"  minimum vertex degree per component: "
              f"{r['min_degree_by_component']}\n")
        if "image_size" not in r:
            w("  (no level n+1 data; stabilization not computed)\n")
            continue
        w(f"  components of F(n+1): {r['n_components_next']}  "
          f"sizes {r['component_sizes_next']}\n")
        w(f"  image of s_{r['n']}: {r['image_size']} component(s)\n")
        w(f"  fibers: {r['fiber_sizes']}   largest = {r['largest_fiber']}\n")
        for tgt, fib in sorted(r["fibers"].items(), key=lambda kv: -len(kv[1])):
            w(f"    component {tgt} at n+1  <-  components {fib} at n\n")
        w(f"  newly born at n+1: {r['n_newborn']}")
        if r["n_newborn"]:
            w(f"  (ids {r['newborn_components']}, "
              f"sizes {r['newborn_sizes']})")
        w("\n")
        w(f"  coalescence distance: {r['coalescence_distance']}")
        if r["coalescence_pair"]:
            w(f"  (between components {r['coalescence_pair']} of F({r['n']}))")
        w("\n")
        if "newborn_criteria_agree" in r:
            w(f"  newborn cross-check: fibers {sorted(r['newborn_components'])} "
              f"vs degree-3 {r['newborn_by_degree_criterion']} -- "
              f"{'AGREE' if r['newborn_criteria_agree'] else 'DISAGREE'}\n")
        if r.get("welldefined_checks"):
            w(f"  lift-fiber check: {r['welldefined_checks']} subdivisions, "
              f"{r['welldefined_violations']} violations\n")
        if r["absorption_one_step"]:
            w(f"  ALL components merge at n+1 (absorption time 1)\n")


def to_csv(records, path):
    import csv
    cols = ["surface", "chi", "orientable", "n", "n_triangulations",
            "n_components", "component_sizes", "min_vertices", "at_minimum_n",
            "n_newborn_here", "newborn_here_sizes", "min_degree_by_component",
            "n_components_next", "image_size", "fiber_sizes", "largest_fiber",
            "n_newborn", "newborn_sizes", "newborn_criteria_agree",
            "coalescence_distance", "absorption_one_step",
            "welldefined_checks", "welldefined_violations"]
    with open(path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        wr.writeheader()
        for r in records:
            wr.writerow({k: (json.dumps(v) if isinstance(v, list) else v)
                         for k, v in r.items() if k in cols})


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        print("usage: python survey.py <lutz_file> [<lutz_file> ...] "
              "[--json out.json] [--csv out.csv] [--detail]")
        sys.exit(1)
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = sys.argv[1:]
    files = []
    i = 0
    skip = False
    for j, a in enumerate(sys.argv[1:]):
        if skip:
            skip = False
            continue
        if a in ("--json", "--csv"):
            skip = True
        elif a == "--detail":
            pass
        else:
            files.append(a)

    recs = run_survey(files)
    print_table(recs)
    if "--detail" in opts:
        print_detail(recs)
    if "--json" in opts:
        p = opts[opts.index("--json") + 1]
        with open(p, "w") as fh:
            json.dump(recs, fh, indent=1)
        print(f"\nwrote {p}")
    if "--csv" in opts:
        p = opts[opts.index("--csv") + 1]
        to_csv(recs, p)
        print(f"wrote {p}")
