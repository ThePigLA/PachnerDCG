"""
Batch driver: run the component-tree survey over a directory of Lutz files.

Handles both of Lutz's naming conventions

    manifolds_lex_d2_n<N>_o<0|1>_g<G>.txt     one surface per file
    manifold_lex_d2_n<N>_o<0|1>.txt           all surfaces of that
                                              orientability at that level

and groups triangulations by their COMPUTED (chi, orientability), not by the
filename -- the filename is used only as a cross-check, and a mismatch is
reported.

Coverage
--------
Fibers of s_n need complete levels n AND n+1 for the same surface.  A pile of
files may cover many surfaces without containing a single consecutive pair,
so this prints a coverage matrix first and names the exact files that would
close each gap.  Per-level component structure is computed for every level
present, consecutive or not.

Usage
-----
    python batch.py <dir-or-files...> [options]

    --coverage-only   print the coverage matrix and exit
    --json PATH       write full records
    --csv PATH        write the summary table
    --detail          expand every row
    --cache PATH      cache parsed levels (default .levelcache.json)
    --no-verify       skip the lift-fiber well-definedness checks
    --max-n N         ignore levels above N
"""

from __future__ import annotations

import json
import os
import re
import sys
import time
from collections import defaultdict

from coalesce import (canonical_form, vertices_of, f_vector,
                      euler_characteristic, is_closed_surface, is_orientable,
                      parse_triangulations_file)
from survey import (surface_name, min_vertices, survey_surface, print_table,
                    print_detail, to_csv)


FNAME = re.compile(
    r"manifolds?_lex_d(?P<d>\d+)_n(?P<n>\d+)_o(?P<o>[01])(?:_g(?P<g>\d+))?",
    re.I)


def parse_filename(path):
    """Metadata declared by the filename, or None."""
    m = FNAME.search(os.path.basename(path))
    if not m:
        return None
    return {
        "d": int(m.group("d")),
        "n": int(m.group("n")),
        "orientable": m.group("o") == "1",
        "genus": int(m.group("g")) if m.group("g") else None,
    }


def genus_of(chi, orientable):
    return (2 - chi) // 2 if orientable else (2 - chi)


def expected_filename(chi, orientable, n):
    return (f"manifolds_lex_d2_n{n}_o{1 if orientable else 0}"
            f"_g{genus_of(chi, orientable)}.txt")


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------

# Bump whenever the on-disk cache layout changes; older files are then
# ignored and rebuilt rather than misread.
CACHE_VERSION = 2


def _valid_cache_entry(entries):
    """A cache entry must be {group_key: {"w": int, "d": hex string}}."""
    if not isinstance(entries, dict):
        return False
    for v in entries.values():
        if not (isinstance(v, dict) and isinstance(v.get("w"), int)
                and isinstance(v.get("d"), str)):
            return False
    return True


def _hexpack(tri):
    """Canonical triangulation -> hex string, 3 bytes (6 hex chars) per face."""
    return bytes(v for t in tri for v in t).hex()


def _hexunpack_all(blob, width):
    """Inverse of _hexpack over a concatenated blob; `width` = bytes per item."""
    step = 2 * width
    out = []
    for i in range(0, len(blob), step):
        b = bytes.fromhex(blob[i:i + step])
        out.append(tuple((b[j], b[j + 1], b[j + 2]) for j in range(0, len(b), 3)))
    return out


def collect_files(args):
    files = []
    for a in args:
        if os.path.isdir(a):
            for fn in sorted(os.listdir(a)):
                if fn.endswith(".txt") and FNAME.search(fn):
                    files.append(os.path.join(a, fn))
        elif os.path.isfile(a):
            files.append(a)
        else:
            print(f"  skipping {a}: not found", file=sys.stderr)
    return files


def load_all(files, cache_path=None, max_n=None, verbose=True):
    """
    Returns (levels, problems) where
      levels[(chi, orientable)][n] = sorted list of canonical triangulations
      problems = list of human-readable warnings
    """
    cache = {}
    if cache_path and os.path.exists(cache_path):
        try:
            with open(cache_path) as fh:
                blob = json.load(fh)
            if (isinstance(blob, dict)
                    and blob.get("cache_version") == CACHE_VERSION
                    and isinstance(blob.get("entries"), dict)):
                cache = blob["entries"]
            else:
                if verbose:
                    print(f"  cache {cache_path}: format from an older "
                          f"version, rebuilding", file=sys.stderr)
        except Exception as exc:
            if verbose:
                print(f"  cache {cache_path}: unreadable ({exc}), rebuilding",
                      file=sys.stderr)

    levels = defaultdict(lambda: defaultdict(set))
    problems = []
    new_cache = {}

    for path in files:
        meta = parse_filename(path)
        stat = os.stat(path)
        key = f"{os.path.basename(path)}:{stat.st_size}:{int(stat.st_mtime)}"

        entries = cache.get(key)
        if entries is not None and not _valid_cache_entry(entries):
            entries = None          # stale or malformed: re-parse
        if entries is not None:
            if verbose:
                print(f"  {os.path.basename(path)}: cached "
                      f"({len(entries)} groups)", file=sys.stderr)
        else:
            t0 = time.time()
            raw = parse_triangulations_file(path)
            groups = defaultdict(list)
            skipped = 0
            for idx, t in enumerate(raw):
                if verbose and idx and idx % 20000 == 0:
                    rate = idx / max(time.time() - t0, 1e-9)
                    print(f"      {idx}/{len(raw)} parsed ({rate:.0f}/s, "
                          f"~{(len(raw) - idx) / max(rate, 1e-9):.0f}s left)",
                          file=sys.stderr)
                if not is_closed_surface(t):
                    skipped += 1
                    continue
                chi = euler_characteristic(t)
                ori = is_orientable(t)
                n = len(vertices_of(t))
                groups[f"{chi}|{int(ori)}|{n}"].append(canonical_form(t))
            # store hex-packed (3 bytes per triangle) rather than nested
            # lists: ~4 KB per triangulation becomes ~170 B, which is the
            # difference between 2 GB and 0.1 GB on a half-million-complex
            # level.
            entries = {k: {"w": 3 * len(v[0]), "d": "".join(_hexpack(t) for t in v)}
                       for k, v in groups.items() if v}
            new_cache[key] = entries
            if verbose:
                n_surf = sum(len(v["d"]) // (2 * v["w"]) for v in entries.values())
                print(f"  {os.path.basename(path)}: {len(raw)} blocks, "
                      f"{n_surf} surfaces"
                      + (f", {skipped} skipped" if skipped else "")
                      + f" ({time.time() - t0:.1f}s)", file=sys.stderr)
            if skipped:
                problems.append(f"{os.path.basename(path)}: {skipped} blocks "
                                f"were not closed surfaces")

        for k, blob in entries.items():
            chi_s, ori_s, n_s = k.split("|")
            chi, ori, n = int(chi_s), bool(int(ori_s)), int(n_s)
            if max_n and n > max_n:
                continue
            for t in _hexunpack_all(blob["d"], blob["w"]):
                levels[(chi, ori)][n].add(t)

        # cross-check the filename against what the file actually contains
        if meta:
            found = {(int(k.split("|")[0]), bool(int(k.split("|")[1])),
                      int(k.split("|")[2])) for k in entries}
            for chi, ori, n in found:
                if meta["n"] != n:
                    problems.append(
                        f"{os.path.basename(path)}: declares n={meta['n']} "
                        f"but contains {n}-vertex triangulations")
                if meta["orientable"] != ori:
                    problems.append(
                        f"{os.path.basename(path)}: declares "
                        f"o={int(meta['orientable'])} but contains "
                        f"{'orientable' if ori else 'nonorientable'} surfaces")
                if meta["genus"] is not None and meta["genus"] != genus_of(chi, ori):
                    problems.append(
                        f"{os.path.basename(path)}: declares g={meta['genus']} "
                        f"but contains genus {genus_of(chi, ori)} "
                        f"({surface_name(chi, ori)})")

    if cache_path and new_cache:
        cache.update(new_cache)
        tmp = cache_path + ".tmp"
        with open(tmp, "w") as fh:
            json.dump({"cache_version": CACHE_VERSION, "entries": cache}, fh)
        os.replace(tmp, cache_path)

    out = {k: {n: sorted(v) for n, v in sorted(lv.items())}
           for k, lv in levels.items()}
    return out, problems


# --------------------------------------------------------------------------
# coverage
# --------------------------------------------------------------------------

def coverage_report(levels, stream=sys.stdout):
    """Print what is present, what is computable, and what is missing."""
    w = stream.write
    if not levels:
        w("no surfaces found\n")
        return [], []

    all_n = sorted({n for lv in levels.values() for n in lv})
    keys = sorted(levels, key=lambda k: (not k[1], -k[0]))

    w("\nCOVERAGE  (levels present per surface)\n")
    head = f"{'surface':>9} {'chi':>4} " + "".join(f"{n:>7}" for n in all_n)
    w(head + "\n" + "-" * len(head) + "\n")
    for k in keys:
        chi, ori = k
        row = f"{surface_name(chi, ori):>9} {chi:>4} "
        for n in all_n:
            row += f"{len(levels[k][n]):>7}" if n in levels[k] else f"{'.':>7}"
        w(row + "\n")

    computable, missing = [], []
    for k in keys:
        chi, ori = k
        ns = sorted(levels[k])
        for n in ns:
            if n + 1 in levels[k]:
                computable.append((k, n))
            else:
                missing.append((k, n, expected_filename(chi, ori, n + 1)))

    w(f"\nFIBERS COMPUTABLE: {len(computable)}\n")
    for (chi, ori), n in computable:
        w(f"  {surface_name(chi, ori):>9}  s_{n} : n={n} -> n={n+1}\n")
    if not computable:
        w("  (none -- no surface has two consecutive levels)\n")

    w(f"\nGAPS: {len(missing)} level(s) have no successor\n")
    for (chi, ori), n, fname in missing:
        w(f"  {surface_name(chi, ori):>9}  n={n} present, n={n+1} missing"
          f"   -> {fname}\n")

    w("\nAlso worth having, to extend each chain downward:\n")
    for k in keys:
        chi, ori = k
        lo = min(levels[k])
        floor = min_vertices(chi, ori)
        if lo > floor:
            wanted = [expected_filename(chi, ori, m) for m in range(floor, lo)]
            w(f"  {surface_name(chi, ori):>9}  has n>={lo}, minimum is "
              f"n={floor}: {', '.join(wanted)}\n")
    w("\n")
    return computable, missing


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def main(argv):
    opts = set(a for a in argv if a.startswith("--"))
    def val(flag, default=None):
        return argv[argv.index(flag) + 1] if flag in argv else default

    skip = set()
    for flag in ("--json", "--csv", "--cache", "--max-n"):
        if flag in argv:
            skip.add(argv.index(flag) + 1)
    args = [a for i, a in enumerate(argv)
            if not a.startswith("--") and i not in skip]

    files = collect_files(args)
    if not files:
        print(__doc__)
        return 1
    print(f"Loading {len(files)} file(s) ...", file=sys.stderr)

    max_n = int(val("--max-n")) if "--max-n" in argv else None
    cache = val("--cache", ".levelcache.json")
    levels, problems = load_all(files, cache_path=cache, max_n=max_n)

    if problems:
        print("\nWARNINGS")
        for p in sorted(set(problems)):
            print(f"  ! {p}")

    computable, missing = coverage_report(levels)

    if "--coverage-only" in opts:
        return 0
    if not computable:
        print("Nothing to stabilize. Computing per-level component structure "
              "only.\n")

    records = []
    for k in sorted(levels, key=lambda k: (not k[1], -k[0])):
        chi, ori = k
        records.extend(survey_surface(chi, ori, levels[k],
                                      verify="--no-verify" not in opts))

    print()
    print_table(records)
    if "--detail" in opts:
        print_detail(records)

    if "--json" in argv:
        with open(val("--json"), "w") as fh:
            json.dump(records, fh, indent=1)
        print(f"\nwrote {val('--json')}")
    if "--csv" in argv:
        to_csv(records, val("--csv"))
        print(f"wrote {val('--csv')}")

    tot = sum(r.get("welldefined_checks", 0) for r in records)
    bad = sum(r.get("welldefined_violations", 0) for r in records)
    if tot:
        print(f"\nlift-fiber lemma: {tot} subdivisions checked, "
              f"{bad} violations")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
