#!/usr/bin/env python3
"""Strict, dimension-generic Pachner certificate helpers.

A d-dimensional bistellar move lives on d+2 vertices U.  If alpha and beta
are disjoint simplices with alpha union beta = U, the move is

    alpha * boundary(beta)  ->  boundary(alpha) * beta.

At the facet-list level this removes {U - {b}: b in beta} and adds
{U - {a}: a in alpha}.  Legality requires:

* the removed side is exactly the live star of alpha; and
* beta is not already a face of the current complex.

The implementation below is deliberately independent of any search engine.
It is suitable for replaying both 3-dimensional and 4-dimensional
certificates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, FrozenSet, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

Facet = FrozenSet[int]
Face = FrozenSet[int]


def sorted_facet(F: Iterable[int]) -> Tuple[int, ...]:
    return tuple(sorted(int(v) for v in F))


def normalize_facets(facets: Iterable[Iterable[int]], dimension: Optional[int] = None) -> Set[Facet]:
    out = {frozenset(int(v) for v in F) for F in facets}
    if not out:
        raise ValueError("facet list is empty")
    sizes = {len(F) for F in out}
    if len(sizes) != 1:
        raise ValueError(f"mixed facet sizes: {sorted(sizes)}")
    facet_size = next(iter(sizes))
    if dimension is not None and facet_size != dimension + 1:
        raise ValueError(
            f"expected facets of size {dimension + 1}, found {facet_size}"
        )
    return out


def facets_jsonable(facets: Iterable[Facet]) -> List[List[int]]:
    return [list(sorted(F)) for F in sorted(facets, key=lambda F: tuple(sorted(F)))]


def canonical_complex_bytes(facets: Iterable[Facet]) -> bytes:
    return json.dumps(
        facets_jsonable(facets),
        separators=(",", ":"),
        sort_keys=False,
    ).encode("utf-8")


def canonical_complex_sha256(facets: Iterable[Facet]) -> str:
    return hashlib.sha256(canonical_complex_bytes(facets)).hexdigest()


def load_facets(path: Path, dimension: Optional[int] = None) -> Set[Facet]:
    data = json.loads(path.read_text())
    if isinstance(data, dict):
        if "facets" in data:
            data = data["facets"]
        elif "final_facets" in data:
            data = data["final_facets"]
        else:
            raise ValueError(f"cannot locate facet list in {path}")
    return normalize_facets(data, dimension)


def save_facets(path: Path, facets: Iterable[Facet]) -> None:
    path.write_text(json.dumps(facets_jsonable(facets), indent=2) + "\n")


def vertices(facets: Iterable[Facet]) -> Set[int]:
    ans: Set[int] = set()
    for F in facets:
        ans.update(F)
    return ans


def face_index(facets: Iterable[Facet], max_size: Optional[int] = None) -> Dict[Face, Set[Facet]]:
    fs = list(facets)
    facet_size = len(fs[0])
    if max_size is None:
        max_size = facet_size
    idx: Dict[Face, Set[Facet]] = defaultdict(set)
    for F in fs:
        items = sorted(F)
        for size in range(1, min(max_size, facet_size) + 1):
            for face in combinations(items, size):
                idx[frozenset(face)].add(F)
    return idx


def live_star(facets: Set[Facet], face: Iterable[int]) -> Set[Facet]:
    A = frozenset(face)
    if not A:
        return set(facets)
    return {F for F in facets if A <= F}


def has_face(facets: Set[Facet], face: Iterable[int]) -> bool:
    A = frozenset(face)
    return any(A <= F for F in facets)


def f_vector(facets: Iterable[Facet]) -> List[int]:
    fs = list(facets)
    facet_size = len(fs[0])
    return [
        len({frozenset(face) for F in fs for face in combinations(sorted(F), size)})
        for size in range(1, facet_size + 1)
    ]


def euler_characteristic(facets: Iterable[Facet]) -> int:
    fv = f_vector(facets)
    return sum((1 if i % 2 == 0 else -1) * n for i, n in enumerate(fv))


def check_closed_pseudomanifold(facets: Set[Facet]) -> Tuple[bool, Dict[Face, int]]:
    facet_size = len(next(iter(facets)))
    counts: Counter[Face] = Counter()
    for F in facets:
        for ridge in combinations(sorted(F), facet_size - 1):
            counts[frozenset(ridge)] += 1
    bad = {R: c for R, c in counts.items() if c != 2}
    return not bad, bad


def check_connected(facets: Set[Facet]) -> bool:
    if not facets:
        return False
    facet_size = len(next(iter(facets)))
    ridge_to_facets: Dict[Face, List[Facet]] = defaultdict(list)
    for F in facets:
        for ridge in combinations(sorted(F), facet_size - 1):
            ridge_to_facets[frozenset(ridge)].append(F)
    adj: Dict[Facet, Set[Facet]] = defaultdict(set)
    for incident in ridge_to_facets.values():
        for A in incident:
            for B in incident:
                if A != B:
                    adj[A].add(B)
    root = next(iter(facets))
    seen = {root}
    q = deque([root])
    while q:
        F = q.popleft()
        for G in adj[F]:
            if G not in seen:
                seen.add(G)
                q.append(G)
    return len(seen) == len(facets)


def _ridge_induced_sign(F: Tuple[int, ...], omitted: int, facet_sign: int) -> int:
    return facet_sign * (-1 if omitted % 2 else 1)


def check_orientable(facets: Set[Facet]) -> bool:
    """Check coherent orientability of a closed pure simplicial complex."""
    tuples = [tuple(sorted(F)) for F in facets]
    facet_size = len(tuples[0])
    ridge_to_entries: Dict[Tuple[int, ...], List[Tuple[int, int]]] = defaultdict(list)
    for i, F in enumerate(tuples):
        for omitted in range(facet_size):
            ridge = F[:omitted] + F[omitted + 1 :]
            ridge_to_entries[ridge].append((i, omitted))
    if any(len(entries) != 2 for entries in ridge_to_entries.values()):
        return False

    signs: Dict[int, int] = {}
    for start in range(len(tuples)):
        if start in signs:
            continue
        signs[start] = 1
        q = deque([start])
        while q:
            i = q.popleft()
            F = tuples[i]
            for omitted in range(facet_size):
                ridge = F[:omitted] + F[omitted + 1 :]
                entries = ridge_to_entries[ridge]
                (a, oa), (b, ob) = entries
                j, oj = (b, ob) if a == i else (a, oa)
                # induced orientations on a common ridge must be opposite
                required = -_ridge_induced_sign(F, omitted, signs[i])
                # solve required = sign_j * (-1)^oj
                sign_j = required * (-1 if oj % 2 else 1)
                if j in signs and signs[j] != sign_j:
                    return False
                if j not in signs:
                    signs[j] = sign_j
                    q.append(j)
    return True


@dataclass(frozen=True)
class PachnerMove:
    dimension: int
    removed: Tuple[Facet, ...]
    added: Tuple[Facet, ...]
    alpha: Face
    beta: Face

    @property
    def support(self) -> Face:
        return self.alpha | self.beta

    @property
    def move_type(self) -> str:
        return f"{len(self.removed)}-{len(self.added)}"

    def to_json(self, **extra: Any) -> Dict[str, Any]:
        out: Dict[str, Any] = {
            "type": self.move_type,
            "removed": facets_jsonable(self.removed),
            "added": facets_jsonable(self.added),
            "support": sorted(self.support),
            "alpha": sorted(self.alpha),
            "beta": sorted(self.beta),
        }
        out.update(extra)
        return out


def move_from_alpha_beta(dimension: int, alpha: Iterable[int], beta: Iterable[int]) -> PachnerMove:
    A = frozenset(int(v) for v in alpha)
    B = frozenset(int(v) for v in beta)
    if A & B:
        raise ValueError("alpha and beta must be disjoint")
    U = A | B
    if len(U) != dimension + 2:
        raise ValueError(
            f"dimension {dimension} move needs {dimension + 2} support vertices; "
            f"found {len(U)}"
        )
    removed = tuple(sorted((U - {b} for b in B), key=lambda F: tuple(sorted(F))))
    added = tuple(sorted((U - {a} for a in A), key=lambda F: tuple(sorted(F))))
    return PachnerMove(dimension, removed, added, A, B)


def infer_move(
    dimension: int,
    removed: Iterable[Iterable[int]],
    added: Iterable[Iterable[int]],
) -> PachnerMove:
    R = tuple(normalize_facets(removed, dimension))
    Aadd = tuple(normalize_facets(added, dimension))
    all_facets = R + Aadd
    U: Set[int] = set()
    for F in all_facets:
        U.update(F)
    if len(U) != dimension + 2:
        raise ValueError(
            f"move support has {len(U)} vertices, expected {dimension + 2}"
        )
    alpha = set.intersection(*(set(F) for F in R)) if R else set()
    beta = set(U) - alpha
    move = move_from_alpha_beta(dimension, alpha, beta)
    if set(move.removed) != set(R) or set(move.added) != set(Aadd):
        raise ValueError("removed/added facets are not complementary Pachner sides")
    return move


def move_legality(facets: Set[Facet], move: PachnerMove) -> Tuple[bool, str]:
    if not set(move.removed) <= facets:
        return False, "removed-facets-not-live"
    expected_star = set(move.removed)
    actual_star = live_star(facets, move.alpha)
    if actual_star != expected_star:
        return False, "removed-side-not-exact-live-star"
    if has_face(facets, move.beta):
        return False, "opposite-face-has-live-cofaces"
    if set(move.added) & facets:
        return False, "added-facet-already-live"
    return True, "ok"


def apply_move(facets: Set[Facet], move: PachnerMove, strict: bool = True) -> Set[Facet]:
    if strict:
        ok, reason = move_legality(facets, move)
        if not ok:
            raise ValueError(f"illegal {move.move_type} move: {reason}")
    return (set(facets) - set(move.removed)) | set(move.added)


def certificate_document(
    dimension: int,
    input_facets: Set[Facet],
    moves: Sequence[Mapping[str, Any]],
    final_facets: Set[Facet],
    metadata: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    return {
        "format": "pachner-certificate-v1",
        "dimension": dimension,
        "input_complex_sha256": canonical_complex_sha256(input_facets),
        "num_moves": len(moves),
        "moves": list(moves),
        "final_complex_sha256": canonical_complex_sha256(final_facets),
        "final_facets": facets_jsonable(final_facets),
        "metadata": dict(metadata or {}),
    }


def verify_certificate(
    input_facets: Set[Facet],
    cert: Mapping[str, Any],
    paranoid: bool = False,
) -> Tuple[Set[Facet], Counter[str]]:
    dimension = int(cert["dimension"])
    expected_hash = cert.get("input_complex_sha256")
    actual_hash = canonical_complex_sha256(input_facets)
    if expected_hash and expected_hash != actual_hash:
        raise ValueError(
            f"input complex hash mismatch: cert={expected_hash} actual={actual_hash}"
        )
    facets = normalize_facets(input_facets, dimension)
    counts: Counter[str] = Counter()
    for i, row in enumerate(cert.get("moves", [])):
        try:
            move = infer_move(dimension, row["removed"], row["added"])
        except Exception as exc:
            raise ValueError(f"move {i}: malformed: {exc}") from exc
        declared = row.get("type")
        if declared and declared != move.move_type:
            raise ValueError(
                f"move {i}: declared type {declared}, inferred {move.move_type}"
            )
        ok, reason = move_legality(facets, move)
        if not ok:
            raise ValueError(f"move {i}: illegal {move.move_type}: {reason}")
        facets = apply_move(facets, move, strict=False)
        counts[move.move_type] += 1
        if paranoid:
            pm, bad = check_closed_pseudomanifold(facets)
            if not pm:
                raise ValueError(f"move {i}: pseudomanifold failure at {len(bad)} ridges")
            if not check_connected(facets):
                raise ValueError(f"move {i}: complex disconnected")
            if not check_orientable(facets):
                raise ValueError(f"move {i}: complex nonorientable")
    expected_final = cert.get("final_facets")
    if expected_final is not None:
        expected = normalize_facets(expected_final, dimension)
        if facets != expected:
            raise ValueError("final complex does not match certificate")
    expected_final_hash = cert.get("final_complex_sha256")
    actual_final_hash = canonical_complex_sha256(facets)
    if expected_final_hash and expected_final_hash != actual_final_hash:
        raise ValueError("final complex hash does not match certificate")
    return facets, counts


def summarize_complex(facets: Set[Facet]) -> Dict[str, Any]:
    d = len(next(iter(facets))) - 1
    pm, bad = check_closed_pseudomanifold(facets)
    return {
        "dimension": d,
        "vertices": len(vertices(facets)),
        "facets": len(facets),
        "f_vector": f_vector(facets),
        "euler_characteristic": euler_characteristic(facets),
        "closed_pseudomanifold": pm,
        "bad_ridges": len(bad),
        "connected": check_connected(facets),
        "orientable": check_orientable(facets),
        "canonical_sha256": canonical_complex_sha256(facets),
    }


def cli() -> None:
    ap = argparse.ArgumentParser(description="Verify a strict Pachner certificate")
    ap.add_argument("input", type=Path)
    ap.add_argument("certificate", type=Path)
    ap.add_argument("--paranoid", action="store_true")
    args = ap.parse_args()

    cert = json.loads(args.certificate.read_text())
    dimension = int(cert["dimension"])
    input_facets = load_facets(args.input, dimension)
    print("input:", json.dumps(summarize_complex(input_facets), sort_keys=True))
    final, counts = verify_certificate(input_facets, cert, paranoid=args.paranoid)
    print(f"verified {sum(counts.values())} moves: {dict(sorted(counts.items()))}")
    print("final:", json.dumps(summarize_complex(final), sort_keys=True))
    print("VERDICT: certificate VALID")


if __name__ == "__main__":
    cli()
