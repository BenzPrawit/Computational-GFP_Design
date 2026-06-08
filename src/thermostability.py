"""Thermostability prior for sfGFP-scaffold designs.

Encodes literature-derived mutations and structural heuristics that correlate
with improved folding stability and resistance to heat denaturation.

References (in repo: referencepaper/):
  * Pédelacq et al. 2006 — Superfolder GFP (sfGFP)
  * Don Paul et al. 2013 — TGP, an extremely stable non-aggregating FP
  * Zacharias et al. 2002 — A206K monomerizing mutation
  * Costantini et al. 2015 — robust folder mutations on sfGFP scaffold
  * Hostettler et al. 2017 — mGreenLantern stability
  * StayGold (Hirano 2022), mBaoJin — bright dimeric/monomeric variants

POSITION NUMBERING NOTE
-----------------------
All positions in this module are 1-indexed canonical avGFP/sfGFP numbering
(1..238). After the 2026-05-27 competition correction (extra residue removed
from the distributed sfGFP reference), competition numbering == canonical
numbering, so `canon_to_comp` is now an identity function. It is retained
for callers that still invoke it.
"""
from __future__ import annotations
from dataclasses import dataclass
from .constants import SFGFP, CHROMOPHORE_POS_0IDX


def canon_to_comp(canonical_pos: int) -> int:
    """Identity mapping (kept for backwards compatibility).

    Before 2026-05-27 the competition sfGFP carried an erroneous extra residue
    at position 171, requiring a +1 shift for canonical positions >= 172.
    After the official correction the competition sfGFP is the canonical
    238-aa sequence, so no shift is needed.
    """
    return canonical_pos


@dataclass
class ThermoMutation:
    pos: int          # 1-indexed canonical avGFP/sfGFP position (1..238)
    aa: str           # target amino acid
    weight: float     # confidence weight in (0, 1]
    source: str
    note: str = ""


# ---------------------------------------------------------------------------
# Curated thermostabilizing single mutations on the sfGFP backbone.
# Positions are 1-indexed canonical avGFP/sfGFP numbering (1..238).
# ---------------------------------------------------------------------------
def _C(canon_pos: int, aa: str, weight: float, source: str, note: str = "") -> ThermoMutation:
    """Constructor that takes a canonical position and stores it directly."""
    parent_aa = SFGFP[canon_pos - 1]
    full_note = f"{parent_aa}{canon_pos}{aa} (canonical). {note}"
    return ThermoMutation(canon_pos, aa, weight, source, full_note)


CANDIDATES: list = [
    # Monomerization → less aggregation under heat → better residual brightness
    _C(206, "K", 0.95, "Zacharias 2002",   "A206K monomerizing mutation"),
    # Folding-robustness mutations
    _C(167, "T", 0.85, "mClover/sfGFP lineage", "I167T improves folding kinetics"),
    _C(149, "K", 0.65, "consensus charge", "Surface charge consensus near barrel"),
    _C(147, "T", 0.55, "consensus",        "Surface stabilizing"),
    _C(177, "F", 0.55, "TGP-inspired",     "Hydrophobic core packing"),
    _C(223, "I", 0.60, "TGP-inspired",     "Barrel hydrophobic packing"),
    _C(224, "L", 0.55, "TGP-inspired",     "Barrel hydrophobic packing"),
    _C(46,  "L", 0.50, "robust folder",    "F46L stabilizes loop"),
    _C(80,  "R", 0.50, "robust folder",    "Q80R surface charge"),
    _C(231, "L", 0.45, "robust folder",    "T231L barrel packing"),
    _C(78,  "K", 0.40, "consensus charge", "Q78K surface ion pair"),
    _C(176, "K", 0.40, "consensus charge", "N176K surface lysine"),
    _C(154, "L", 0.45, "TGP-like",         "Improved core packing"),
    # Additional consensus mutations gleaned from comparative GFP studies
    _C(43,  "L", 0.35, "consensus",        "F43L loop stabilization"),
    _C(141, "L", 0.30, "consensus",        "I141L barrel packing"),
    _C(195, "I", 0.30, "consensus",        "V195I core packing"),
    # Avoid Cys-based disulfide bonds because the cell-free E. coli system is reducing
    # and disulfides can mis-fold; we explicitly do not include any C mutations here.
]


def candidate_pool() -> list:
    """Return [(pos_comp, aa, weight)] usable by the mutation engine."""
    return [(c.pos, c.aa, c.weight) for c in CANDIDATES
            if (c.pos - 1) not in CHROMOPHORE_POS_0IDX]


def thermostability_score(seq: str, parent: str = SFGFP) -> float:
    """Sum weights of supported thermostabilizing mutations present in `seq`,
    minus a small penalty for unsupported deviations from the parent.

    This is a *prior* used to bias the cascade-funnel search, not a predictor.
    """
    if len(seq) != len(parent):
        return 0.0
    by_pos = {c.pos: c for c in CANDIDATES}
    s = 0.0
    for i, (a, b) in enumerate(zip(parent, seq), 1):
        if a == b:
            continue
        c = by_pos.get(i)
        if c and c.aa == b:
            s += c.weight
        else:
            s -= 0.05  # unsupported drift
    return s


def explain(seq: str, parent: str = SFGFP) -> list:
    """Return human-readable list of supported mutations present in seq."""
    by_pos = {c.pos: c for c in CANDIDATES}
    out = []
    for i, (a, b) in enumerate(zip(parent, seq), 1):
        if a != b:
            c = by_pos.get(i)
            tag = f"{a}{i}{b}"
            if c and c.aa == b:
                out.append(f"{tag} (+{c.weight:.2f}, {c.source})")
            else:
                out.append(f"{tag} (unsupported)")
    return out


if __name__ == "__main__":
    from .mutate import apply_mutations
    print(f"Candidates: {len(CANDIDATES)} (all in competition numbering)")
    print(f"Parent length: {len(SFGFP)}")
    # Apply A206K + I167T + N149K (canonical) → comp positions
    test = apply_mutations(SFGFP, [(c.pos, c.aa) for c in CANDIDATES[:3]])
    print(f"Score (top-3 mutations applied): {thermostability_score(test):.3f}")
    for line in explain(test):
        print(" ", line)
