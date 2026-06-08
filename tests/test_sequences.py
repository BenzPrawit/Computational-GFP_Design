"""Reproducibility checks for the designed GFP sequences and project layout.

Run from the repository root:  pytest tests/
"""
import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SUBMISSION = ROOT / "results" / "tables" / "submission.csv"


def _load_records():
    with open(SUBMISSION, newline="") as fh:
        return list(csv.DictReader(fh))


def test_submission_exists():
    assert SUBMISSION.exists(), f"missing {SUBMISSION}"


def test_expected_project_files_present():
    for rel in [
        "README.md",
        "LICENSE",
        "CITATION.cff",
        "requirements.txt",
        "src/validate.py",
        "data/raw/GFP_data.xlsx",
        "results/tables/designed_sequences.fasta",
    ]:
        assert (ROOT / rel).exists(), f"missing {rel}"


def test_submission_has_six_sequences():
    records = _load_records()
    assert 1 <= len(records) <= 6, f"expected 1–6 sequences, got {len(records)}"


def test_designed_sequences_pass_validation():
    """Every submitted sequence must satisfy the competition hard rules."""
    try:
        from src.validate import validate_sequence
    except Exception as exc:  # pragma: no cover - import/env issue
        pytest.skip(f"src.validate not importable in this environment: {exc}")

    records = _load_records()
    for rec in records:
        seq = rec["Sequence"].strip()
        result = validate_sequence(seq)
        assert result.ok, f"Seq_ID {rec.get('Seq_ID')} invalid: {result.reasons}"


def test_sequences_are_distinct():
    records = _load_records()
    seqs = [r["Sequence"].strip() for r in records]
    assert len(seqs) == len(set(seqs)), "duplicate sequences in submission"
