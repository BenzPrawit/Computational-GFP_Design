"""Fast O(1) exclusion-list membership check via SHA-256 hashing.

The Exclusion_List.csv has ~135K sequences (~30 MB). We hash each entry to a
fixed-size set so the membership check is constant-time and memory-stable.
"""
from __future__ import annotations
import csv
import hashlib
from pathlib import Path
from typing import Iterable


def _h(seq: str) -> str:
    return hashlib.sha256(seq.strip().upper().encode("ascii")).hexdigest()


class ExclusionIndex:
    def __init__(self, hashes: set):
        self._hashes = hashes

    @classmethod
    def from_csv(cls, path: str | Path) -> "ExclusionIndex":
        hashes: set = set()
        path = Path(path)
        with path.open(newline="") as f:
            reader = csv.reader(f)
            header = next(reader, None)  # 'Sequence'
            for row in reader:
                if not row:
                    continue
                seq = row[0].strip().upper()
                if seq:
                    hashes.add(_h(seq))
        return cls(hashes)

    def __len__(self) -> int:
        return len(self._hashes)

    def contains(self, seq: str) -> bool:
        return _h(seq) in self._hashes

    def filter_novel(self, seqs: Iterable[str]) -> list:
        return [s for s in seqs if not self.contains(s)]


if __name__ == "__main__":
    # Smoke test
    idx = ExclusionIndex.from_csv(
        Path(__file__).resolve().parents[1] / "data" / "Exclusion_List.csv"
    )
    print(f"Loaded {len(idx):,} excluded sequences")
    sample = "APAMKIECRITGTLNGVEFELVGGGEGTPEQGRMTNKMKSTKGALTFSPYLLSAVMGYGFYHFGTYPSGYENPFLHAINNGGYTNTRIEKYEDGGVLHVSFSYRYEAGRVIGDFKVVGTGFPEDSVIFTDKIIRSNATVEHLHPMGDNVLVGSFARTFSLRDGGYYSFVVDSHMHFKSARHPSILQNGGPMFAFRRVEELHSNTELGIVEYQHAFKTPIAFA"
    print("Known excluded sample present?", idx.contains(sample))
    print("Random novel sequence present?", idx.contains("MNOVEL" * 40))
