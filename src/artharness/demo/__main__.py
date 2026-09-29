"""python -m artharness.demo: run the whole ART workflow at toy scale and explain it."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .campaign import DemoOptions, run_demo


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="python -m artharness.demo", description=__doc__)
    ap.add_argument("--agent", choices=["rules", "ollama"], default="rules",
                    help="rules: deterministic policy, offline. ollama: a small local model.")
    ap.add_argument("--model", default="qwen3:0.6b", help="Ollama model name")
    ap.add_argument("--ollama-url", default="http://localhost:11434")
    ap.add_argument("--out", type=Path, default=Path("demo-run"))
    ap.add_argument("--force", action="store_true", help="replace an existing --out")
    ap.add_argument("--no-stage-revision", action="store_true",
                    help="rules agent: skip the staged sloppy first pass of stage 4")
    ap.add_argument("--concurrency", type=int, default=1,
                    help="sessions dispatched in parallel (paper: up to 58)")
    ap.add_argument("--quiet", action="store_true", help="no 'why' commentary")
    a = ap.parse_args(argv)
    try:
        res = run_demo(DemoOptions(
            out=a.out, agent=a.agent, model=a.model, ollama_url=a.ollama_url,
            stage_revision=not a.no_stage_revision, concurrency=a.concurrency,
            explain=not a.quiet, force=a.force))
    except (FileExistsError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0 if res.error is None else 1


if __name__ == "__main__":
    raise SystemExit(main())
