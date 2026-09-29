"""Campaign configuration — every constant traces to the paper (see docs/paper-notes.md).

"NOT-IN-PAPER" marks values we had to choose because the paper does not state them.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class CampaignConfig:
    # Harness scale (paper Methods "Autonomous research harness", p.28)
    max_concurrent_sessions: int = 58  # paper: "permitted up to 58 concurrent sessions"
    # NOT-IN-PAPER: sessions actually dispatched in parallel, capped at
    # max_concurrent_sessions. Default 1 keeps offline tests and dry runs
    # deterministic (task ids follow creation order); a live campaign sets it up to
    # max_concurrent_sessions (paper: 58).
    dispatch_concurrency: int = 1
    sandbox_cpus: int = 60  # paper: sandbox with 60 CPU cores
    sandbox_mem_gib: int = 192  # paper: 192 GiB of memory, no GPU

    # Stall handling (paper Methods "Campaign accounting", p.30):
    # one task stalled "after ten revisions" and one "after ten failed completion checks"
    max_revisions: int = 10
    max_gate_failures: int = 10

    # Report tournament (paper Methods "Ranking of reports", p.30)
    judge_weights: dict[str, float] = field(
        default_factory=lambda: {
            "impact": 0.35,
            "novelty": 0.30,
            "soundness": 0.25,
            "actionability": 0.10,
        }
    )
    soundness_auto_lose: int = 2  # "a report with a soundness score of 2 or less lost automatically"

    # Benchmark (paper Methods "Fixed-input benchmark", p.38)
    benchmark_attempts_per_cell: int = 100  # 7 models x 5 levels x 100 = 3,500 attempts
    benchmark_max_output_tokens: int = 1_000_000  # "maximum output budget of 1 million tokens"
    benchmark_levels: tuple[str, ...] = ("L1", "L2", "L3", "L4", "L5")

    # Budgets. NOT-IN-PAPER: the paper states no token/price caps; these are safety rails
    # so a runaway campaign dies instead of billing silently.
    max_tasks_total: int = 500
    max_tokens_total: int = 300_000_000  # paper campaign used 215.6M

    # Session backend gate. Live LLM backends refuse to run unless this env var is set
    # by the owner (README "Safety gates").
    live_backend_env: str = "ARTHARNESS_ALLOW_LIVE"

    # Brief provenance (brief/research_brief.md): every campaign receipt must
    # record that the brief is a reconstruction, never Anthropic's text. Written
    # into the campaign record by run_stage_chain so it cannot be forgotten.
    brief_provenance: str = "brief: reconstruction-v0 (6 anchors), not verbatim"
