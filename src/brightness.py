"""Data-driven brightness prior built from GFP_data.xlsx (141K mutations).

We learn a position-by-AA effect table from the avGFP brightness dataset.
Mutation strings are parsed (e.g. "A109D:N145D:I187V") and for each (pos, aa)
substitution we aggregate the mean log10-brightness across all sequences that
contain it. Subtracting the WT baseline gives a Δlog10-brightness effect.

This is the "independent-site approximation": we ignore epistasis. It's not
perfect, but it's an unbiased prior that we use only as a *ranker*, not as
ground truth. The cell-free wet-lab measurement is the true judge.

Position numbering in the source data is canonical avGFP (1..238). After the
2026-05-27 competition correction, competition sfGFP is also canonical 238 aa,
so `thermostability.canon_to_comp` is now an identity function.
"""
from __future__ import annotations
import json
from pathlib import Path
from collections import defaultdict
import pandas as pd
from .constants import SFGFP, AVGFP
from .thermostability import canon_to_comp

DATA_PATH = Path(__file__).resolve().parents[1] / "data" / "GFP_data.xlsx"
CACHE_PATH = Path(__file__).resolve().parent / "_brightness_cache.json"
WT_LABEL = "WT"


def _parse_mut_string(s: str) -> list:
    """Parse 'A109D:N145D' into [(109,'A','D'), (145,'N','D')]. Returns [] for WT."""
    if not isinstance(s, str) or s.upper() == WT_LABEL or not s.strip():
        return []
    out = []
    for tok in s.split(":"):
        tok = tok.strip()
        if len(tok) < 3:
            continue
        try:
            ref, pos, alt = tok[0], int(tok[1:-1]), tok[-1]
            out.append((pos, ref, alt))
        except (ValueError, IndexError):
            continue
    return out


def build_effect_table(family: str = "avGFP", min_support: int = 3,
                       func_threshold: float = 2.5) -> dict:
    """Build a per-mutation prior table that is robust to the data's floor effect.

    For each (pos, alt) substitution we record:
      * `effect`        — Δlog10-brightness using the mean of the TOP 25% of
                          brightness values among sequences containing the
                          mutation (best-case effect, less polluted by
                          catastrophic combinations).
      * `func_rate`     — fraction of sequences containing the mutation that
                          retain functional brightness (log10 ≥ func_threshold,
                          i.e. ~316 linear ≥ ~50 % of avGFP WT). High = the
                          mutation is permissive / non-destructive.
      * `support`       — sample count.

    Rationale: most rows in the dataset combine many mutations and the
    aggregate sometimes lands at the brightness floor (~1.30 log10) regardless
    of any individual mutation's true contribution. Mean-of-all is dragged
    down by these. The top-quartile mean and functional-retention rate are
    much better single-mutation priors.
    """
    df = pd.read_excel(DATA_PATH, sheet_name="brightness")
    df = df[df["GFP type"] == family].copy()
    if df.empty:
        raise RuntimeError(f"No rows for family {family}")
    wt_brightness = float(df[df["aaMutations"].astype(str).str.upper() == WT_LABEL]["Brightness"].mean())

    bucket = defaultdict(list)
    for _, row in df.iterrows():
        muts = _parse_mut_string(str(row["aaMutations"]))
        if not muts:
            continue
        b = float(row["Brightness"])
        for pos, ref, alt in muts:
            bucket[(pos, alt)].append(b)

    table = {}
    for key, vals in bucket.items():
        if len(vals) < min_support:
            continue
        vals_sorted = sorted(vals, reverse=True)
        top_q = vals_sorted[: max(1, len(vals_sorted) // 4)]
        top_mean = sum(top_q) / len(top_q)
        func_rate = sum(1 for v in vals if v >= func_threshold) / len(vals)
        table[key] = {
            "effect": top_mean - wt_brightness,    # Δlog10-brightness (best case)
            "func_rate": func_rate,                # 0..1
            "support": len(vals),
        }
    return {"wt_brightness": wt_brightness, "family": family,
            "n_entries": int(len(table)), "table": table}


def save_table(table_obj: dict, path: Path = CACHE_PATH) -> None:
    # JSON keys must be strings
    serial = {f"{p}|{a}": v for (p, a), v in table_obj["table"].items()}
    payload = {**table_obj, "table": serial}
    path.write_text(json.dumps(payload, indent=2))


def load_table(path: Path = CACHE_PATH) -> dict:
    payload = json.loads(path.read_text())
    payload["table"] = {tuple([int(k.split("|")[0]), k.split("|")[1]]): v
                         for k, v in payload["table"].items()}
    return payload


class BrightnessScorer:
    """Score a sfGFP-numbered sequence by summing per-substitution Δlog10-brightness
    effects looked up from the avGFP dataset. Substitutions not in the table are
    treated as zero-effect (neutral prior) so we don't over-penalize novelty."""

    def __init__(self, table_obj: dict | None = None, family: str = "avGFP"):
        if table_obj is None:
            if CACHE_PATH.exists():
                table_obj = load_table()
            else:
                table_obj = build_effect_table(family=family)
                save_table(table_obj)
        self.wt_brightness = table_obj["wt_brightness"]
        self.table = table_obj["table"]
        self.family = table_obj["family"]

    def score(self, seq: str, parent: str = SFGFP) -> float:
        """Sum the Δlog10-brightness effects of all substitutions in `seq`
        relative to the parent. Positions are 1-indexed canonical avGFP numbering.

        After the 2026-05-27 competition correction, the sfGFP reference is
        the canonical 238-aa sequence, so no position-shift is needed. The
        legacy 239-aa code path is retained for backwards compatibility if a
        caller still passes the old erroneous parent.
        """
        if len(seq) != len(parent):
            return 0.0
        total = 0.0
        is_239_legacy = (len(parent) == 239)  # only true if caller passes the OLD erroneous parent
        for pos_1idx, (a, b) in enumerate(zip(parent, seq), 1):
            if a == b:
                continue
            if is_239_legacy:
                if pos_1idx <= 171:
                    canon_pos = pos_1idx
                elif pos_1idx == 172:
                    continue
                else:
                    canon_pos = pos_1idx - 1
            else:
                canon_pos = pos_1idx  # canonical 238-aa parent → no shift
            entry = self.table.get((canon_pos, b))
            if entry is not None:
                total += float(entry["effect"])
        return total

    def per_site_options(self, top_k: int = 3, min_func_rate: float = 0.5,
                         min_support: int = 5) -> dict:
        """Return {comp_pos: [(aa, effect), ...]} of substitutions that
        (a) appear in functional sequences ≥ `min_func_rate` of the time and
        (b) have positive top-quartile effect, ranked by effect.

        Only the top-k options per site are returned.
        """
        per_pos = defaultdict(list)
        for (canon_pos, aa), v in self.table.items():
            if v.get("support", 0) < min_support:
                continue
            if v.get("func_rate", 0) < min_func_rate:
                continue
            per_pos[canon_pos].append((aa, v["effect"], v["support"], v["func_rate"]))
        out = {}
        for canon_pos, items in per_pos.items():
            items.sort(key=lambda t: t[1], reverse=True)
            comp_pos = canon_to_comp(canon_pos)
            out[comp_pos] = [(aa, eff) for aa, eff, _, _ in items[:top_k] if eff > 0]
        return out


if __name__ == "__main__":
    table = build_effect_table()
    print(f"avGFP WT log10 brightness: {table['wt_brightness']:.3f}")
    print(f"Effect table size: {table['n_entries']:,} (pos, aa) entries")
    save_table(table)
    print(f"Cached → {CACHE_PATH}")

    scorer = BrightnessScorer()
    # Score sfGFP itself relative to itself: should be 0
    print(f"sfGFP self-score: {scorer.score(SFGFP):.3f}")

    # Score the previous Top-1 from 2025 to sanity-check
    top2025 = (
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVLCFSRYPDHMKQH"
        "DFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNPHNVYIMADKQKNGIK"
        "AYFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLDTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    )
    print(f"2025 top seq length: {len(top2025)} (vs sfGFP {len(SFGFP)})")
    if len(top2025) == len(SFGFP):
        print(f"2025-top vs sfGFP brightness prior: {scorer.score(top2025):.3f}")

    options = scorer.per_site_options(top_k=2)
    print(f"\nSites with positive-effect mutations: {len(options)}")
    # Print a few sample positions
    for comp_pos in sorted(options)[:8]:
        print(f"  pos {comp_pos}: {options[comp_pos]}")
