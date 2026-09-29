"""Consensus ties must not depend on hash randomization.

Paper fidelity angle: the repo pins RNG seeds everywhere for reproducibility
("Sampling is seeded (deterministic)"); a repeat string that flips with
PYTHONHASHSEED breaks that guarantee end to end. Raw-fidelity audit 2026-09-29
finding 2 (wf_8e67c019-baf): max(set(col), key=col.count) breaks modal-count
ties by set iteration order.
"""
import os
import subprocess
import sys

COLS = ["AAAA", "ACAC", "AAAA", "ACAC"]


def _consensus_under(seed: str) -> str:
    code = (
        "from artharness.arrays import _consensus_block;"
        f"print(_consensus_block({COLS!r})[2])"
    )
    env = dict(os.environ, PYTHONHASHSEED=seed)
    out = subprocess.run([sys.executable, "-c", code], env=env,
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()


def test_consensus_letter_is_hashseed_independent():
    a = _consensus_under("0")
    b = _consensus_under("3")
    c = _consensus_under("5")
    assert a == b == c, (
        f"consensus differs across PYTHONHASHSEED: {a!r} vs {b!r} vs {c!r} — "
        "modal-count ties must resolve in fixed base order, not set order"
    )


def test_chain_score_consensus_letter_is_hashseed_independent():
    # _chain_score's consensus join has the same max(set(col), ...) pattern
    code = (
        "import random\n"
        "from artharness.arrays import _chain_score\n"
        "win = ('ACACGTACGTAC' * 40) + 'ACGT'\n"
        "chain = [0, 12, 24, 36, 48]\n"
        "print(_chain_score(win, chain, 8)[1])\n"
    )
    outs = []
    for seed in ("0", "3", "5"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        r = subprocess.run([sys.executable, "-c", code], env=env,
                           capture_output=True, text=True, check=True)
        outs.append(r.stdout.strip())
    assert len(set(outs)) == 1, (
        f"_chain_score consensus differs across PYTHONHASHSEED: {outs!r}"
    )
