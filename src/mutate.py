"""Mutation-library generators.

Produces candidate variants by combining (a) a fixed set of literature-backed
mutations and (b) a curated pool of position-specific substitutions learned
from the brightness dataset.
"""
from __future__ import annotations
import random
import itertools
from typing import Iterable, Sequence
from .constants import AA20, CHROMOPHORE_POS_0IDX


def apply_mutations(parent: str, muts: Iterable[tuple]) -> str:
    """Apply a list of (pos_1idx, new_aa) mutations to a parent sequence."""
    s = list(parent)
    for pos, aa in muts:
        i = pos - 1
        if 0 <= i < len(s):
            s[i] = aa
    return "".join(s)


def parse_mut(token: str) -> tuple:
    """Parse 'F64L' or 'A206K' into (pos_1idx, new_aa). Returns (pos, aa)."""
    return int(token[1:-1]), token[-1]


def parse_mut_set(spec: str) -> list:
    """Parse 'S30R:Y39N:F64L' into [(30,'R'),(39,'N'),(64,'L')]."""
    return [parse_mut(t) for t in spec.split(":") if t]


def all_combinations(parent: str, mut_pool: Sequence[tuple], max_combo: int = 4) -> list:
    """All ≤max_combo subsets of the mutation pool, applied to parent."""
    out = []
    for k in range(0, max_combo + 1):
        for combo in itertools.combinations(mut_pool, k):
            # Skip combos that touch the same position twice with different AAs
            seen_pos = {}
            ok = True
            for pos, aa in combo:
                if pos in seen_pos and seen_pos[pos] != aa:
                    ok = False
                    break
                seen_pos[pos] = aa
            if ok:
                out.append(apply_mutations(parent, combo))
    return out


def random_combinations(parent: str, mut_pool: Sequence[tuple], n: int,
                        size_range: tuple = (3, 10), seed: int = 0) -> list:
    """Sample n random subsets of the pool with sizes in [lo, hi]."""
    rng = random.Random(seed)
    pool = list(mut_pool)
    out = set()
    lo, hi = size_range
    attempts = 0
    while len(out) < n and attempts < n * 50:
        attempts += 1
        k = rng.randint(lo, min(hi, len(pool)))
        combo = rng.sample(pool, k)
        # de-dup by position
        by_pos = {}
        for pos, aa in combo:
            by_pos[pos] = aa
        seq = apply_mutations(parent, list(by_pos.items()))
        out.add(seq)
    return list(out)


def safe_position(pos_1idx: int) -> bool:
    """Refuse to mutate the chromophore-forming residues (positions 65-67)."""
    return (pos_1idx - 1) not in CHROMOPHORE_POS_0IDX


def filter_safe(muts: Iterable[tuple]) -> list:
    return [(p, a) for p, a in muts if safe_position(p)]
