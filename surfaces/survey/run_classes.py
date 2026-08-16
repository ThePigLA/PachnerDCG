"""
Compute and inspect the class structure of s_n, with checkpointing.

Usage
-----
    python run_classes.py <file> <chi> <orientable:0|1> [options]

Options
-------
    --checkpoint PATH   state file (default: state_chi<chi>_n<n>.json)
    --resume            load the checkpoint instead of rebuilding
    --report-only       print the current class structure and exit
    --rounds N          targeted search rounds (default 50)
    --walks N           walks launched per block per round (default 40)
    --steps N           max flips per walk (default 300)
    --seed N            RNG seed (default 0)
    --no-aut            skip automorphism computation
    --brief             do not print the triangulations themselves

Examples
--------
    # first run, nonorientable genus 4 at n=9
    python run_classes.py manifold_lex_d2_n9_o0.txt -2 0

    # orientable genus 6 at n=12, checkpointed
    python run_classes.py manifold_lex_d2_n12_o1.txt -10 1 \
        --checkpoint sigma6.json --rounds 200 --walks 120

    # inspect what survived, without searching
    python run_classes.py manifold_lex_d2_n12_o1.txt -10 1 \
        --checkpoint sigma6.json --resume --report-only

    # keep hammering the survivors from a saved state
    python run_classes.py manifold_lex_d2_n12_o1.txt -10 1 \
        --checkpoint sigma6.json --resume --rounds 500 --walks 300 --seed 7
"""

import os
import sys
import time

from coalesce import (parse_triangulations_file, select_surface, canonical_form,
                      components_at_level, vertices_of, f_vector)
from fiberstate import FiberState, targeted_search


def parse_args(argv):
    pos, opt = [], {}
    i = 0
    flags = {"--resume", "--report-only", "--no-aut", "--brief"}
    while i < len(argv):
        a = argv[i]
        if a in flags:
            opt[a] = True
            i += 1
        elif a.startswith("--"):
            opt[a] = argv[i + 1]
            i += 2
        else:
            pos.append(a)
            i += 1
    return pos, opt


def build_state(path, chi, orientable):
    print(f"Parsing {path} ...")
    raw = parse_triangulations_file(path)
    print(f"  {len(raw)} blocks extracted")

    target = select_surface(raw, chi=chi, orientable=orientable)
    print(f"  {len(target)} closed surfaces with chi={chi}, "
          f"orientable={orientable}")
    if not target:
        sys.exit("nothing matched the filter")
    print(f"  f-vector {f_vector(target[0])}")

    level = {canonical_form(t) for t in target}
    n = len(vertices_of(next(iter(level))))
    print(f"  {len(level)} isomorphism classes at n={n}")

    print(f"Components of F(n={n}) ...")
    t0 = time.time()
    comps = components_at_level(level)
    comps.sort(key=len, reverse=True)
    print(f"  {len(comps)} components, sizes {[len(c) for c in comps]} "
          f"({time.time() - t0:.1f}s)")

    return FiberState(comps, n, chi, orientable)


def main(argv):
    pos, opt = parse_args(argv)
    if len(pos) < 3:
        print(__doc__)
        return 1
    path, chi, orientable = pos[0], int(pos[1]), bool(int(pos[2]))

    ckpt = opt.get("--checkpoint")
    if opt.get("--resume"):
        if not ckpt or not os.path.exists(ckpt):
            sys.exit("--resume needs an existing --checkpoint file")
        print(f"Resuming from {ckpt} ...")
        state = FiberState.load(ckpt)
        print(f"  {len(state.comps)} components, {state.n_blocks()} blocks, "
              f"{len(state.witnesses)} witnesses")
    else:
        state = build_state(path, chi, orientable)
        if ckpt is None:
            ckpt = f"state_chi{chi}_n{state.n}.json"
        state.save(ckpt)
        print(f"  checkpoint -> {ckpt}")

    if state.aut is None and not opt.get("--no-aut"):
        print("Computing automorphism orders ...")
        t0 = time.time()
        state.compute_aut()
        print(f"  done ({time.time() - t0:.1f}s)")
        state.save(ckpt)

    if not opt.get("--report-only"):
        rounds = int(opt.get("--rounds", 50))
        walks = int(opt.get("--walks", 40))
        steps = int(opt.get("--steps", 300))
        seed = int(opt.get("--seed", 0))
        print(f"\nTargeted search: {rounds} rounds, {walks} walks/block, "
              f"{steps} steps, seed {seed}")
        print(f"  checkpointing every merge -> {ckpt}")
        t0 = time.time()
        targeted_search(state, rounds=rounds, walks_per_block=walks,
                        max_steps=steps, seed=seed, checkpoint=ckpt)
        print(f"  search finished ({time.time() - t0:.1f}s)")
        state.save(ckpt)

    state.report(show_triangulations=not opt.get("--brief"))

    if state.n_blocks() == 1:
        print(f"All {len(state.comps)} components coalesce at "
              f"n={state.n + 1}.  Absorption time = 1.")
    else:
        print(f"{state.n_blocks()} classes survive one stabilization.")
        print("Walk mode cannot prove non-merging: the down-set of a")
        print("component upstairs IS its fiber, so a negative needs")
        print("exhaustion.  Treat this as 'not yet merged'.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
