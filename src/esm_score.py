"""Optional ESM-2 pseudo-log-likelihood scorer (CPU-friendly).

Uses the smallest ESM-2 model (t6_8M, ~8M params) by default — runs at a few
hundred ms per 240-aa sequence on a modern CPU. For each candidate sequence we
compute the model's log-probability of every position summed over the sequence
(a single forward pass). Higher = more "natural-looking" to the protein
language model, which has been shown to correlate with experimental fitness
for many proteins (Meier et al. 2021, Lin et al. 2023).

Importing this module never fails. If `torch` or `fair-esm` aren't installed,
`available()` returns False and the scorer returns 0.0 for every sequence
(safe no-op so the pipeline keeps running).

Install (optional):
    pip install torch fair-esm
"""
from __future__ import annotations

_HAVE_TORCH = False
_HAVE_ESM = False
try:
    import torch
    _HAVE_TORCH = True
except Exception:
    torch = None  # type: ignore

try:
    import esm  # fair-esm
    _HAVE_ESM = True
except Exception:
    esm = None  # type: ignore


def available() -> bool:
    return _HAVE_TORCH and _HAVE_ESM


class ESMScorer:
    """Single-forward-pass log-likelihood scorer using ESM-2 t6_8M."""

    def __init__(self, model_name: str = "esm2_t6_8M_UR50D"):
        if not available():
            self.model = None
            self.alphabet = None
            self.batch_converter = None
            return
        self.model, self.alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        self.model.eval()
        self.batch_converter = self.alphabet.get_batch_converter()

    def score(self, seq: str) -> float:
        """Return mean log-probability of each residue under a single forward pass."""
        if self.model is None:
            return 0.0
        data = [("seq", seq)]
        _, _, tokens = self.batch_converter(data)
        with torch.no_grad():
            out = self.model(tokens, repr_layers=[], return_contacts=False)
        logits = out["logits"][0]                    # [L+2, vocab]
        log_probs = torch.log_softmax(logits, dim=-1)
        # Skip BOS/EOS positions (index 0 and -1)
        per_pos = []
        for i, aa in enumerate(seq, start=1):
            tok_id = self.alphabet.get_idx(aa)
            per_pos.append(float(log_probs[i, tok_id]))
        return sum(per_pos) / len(per_pos)


class _NullScorer:
    """No-op fallback when ESM isn't installed."""
    def score(self, seq: str) -> float:
        return 0.0


def get_scorer():
    return ESMScorer() if available() else _NullScorer()


if __name__ == "__main__":
    print("ESM available:", available())
    s = get_scorer()
    from .constants import SFGFP
    print(f"sfGFP score: {s.score(SFGFP):.4f}")
