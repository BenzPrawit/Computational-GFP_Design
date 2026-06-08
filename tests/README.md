# Tests

Lightweight, dependency-free reproducibility checks that can run in CI without
the heavy ML/MD stack.

```bash
pip install pytest
pytest tests/
```

- `test_sequences.py` — validates the designed sequences in
  `results/tables/submission.csv` against the competition rules implemented in
  `src/validate.py` (length, amino-acid alphabet, M-start), and checks that the
  expected project files exist.

Extend with model/pipeline unit tests as the project matures (e.g. mutation
generation in `src/mutate.py`, exclusion-list logic in `src/exclusion.py`).
