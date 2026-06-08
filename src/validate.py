"""Hard-rule validation for competition submissions."""
from __future__ import annotations
from dataclasses import dataclass
from .constants import AA_SET, MIN_LEN, MAX_LEN


@dataclass
class ValidationResult:
    ok: bool
    reasons: list

    def __bool__(self) -> bool:
        return self.ok


def validate_sequence(seq: str) -> ValidationResult:
    """Apply every hard rule from the competition spec.

    Rules:
      * length 220–250 aa
      * starts with M
      * only the 20 standard uppercase amino acids
      * no stop codons (*) or punctuation
    """
    reasons = []
    if not isinstance(seq, str):
        return ValidationResult(False, ["sequence is not a string"])

    if len(seq) < MIN_LEN or len(seq) > MAX_LEN:
        reasons.append(f"length {len(seq)} not in [{MIN_LEN},{MAX_LEN}]")

    if not seq or seq[0] != "M":
        reasons.append("does not start with M")

    bad = sorted({c for c in seq if c not in AA_SET})
    if bad:
        reasons.append(f"illegal characters present: {bad}")

    if "*" in seq:
        reasons.append("contains stop codon (*)")

    return ValidationResult(len(reasons) == 0, reasons)


def validate_team_submission(records: list, team_name: str, max_seqs: int = 6) -> ValidationResult:
    """Validate the full team submission (≤6 sequences, distinct, all valid)."""
    reasons = []
    if len(records) == 0:
        reasons.append("no sequences submitted")
    if len(records) > max_seqs:
        reasons.append(f"too many sequences: {len(records)} > {max_seqs}")

    seqs_seen = set()
    for i, rec in enumerate(records, 1):
        seq = rec.get("Sequence", "")
        if seq in seqs_seen:
            reasons.append(f"sequence {i} is a duplicate of an earlier submission")
        seqs_seen.add(seq)

        v = validate_sequence(seq)
        if not v.ok:
            reasons.append(f"sequence {i} ({rec.get('Seq_ID')}) failed: {v.reasons}")

        if rec.get("Team_Name") != team_name:
            reasons.append(f"sequence {i} has wrong Team_Name")

    return ValidationResult(len(reasons) == 0, reasons)


if __name__ == "__main__":
    from .constants import SFGFP
    print(validate_sequence(SFGFP))
    print(validate_sequence("MAAA*"))
    print(validate_sequence("XSKGEELFT"))
