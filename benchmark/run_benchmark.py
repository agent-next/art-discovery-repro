"""Benchmark orchestrator (paper Methods "Fixed-input benchmark" p.38).

Attempts = models x levels x N (paper: 7 models x 5 levels x 100 = 3,500
attempts). Each attempt runs in an isolated run dir with a maximum output
budget of 1,000,000 tokens and produces report.md + submission.json.

Grading: a judge model scores the 10 rubric claims -> score 0-10; the headline
metric "recognized repeat array" is the repeat-array claim asserted. The paper's
judge is an LLM; the default judge here is the offline rubric proxy from
rubric.py (paper rule, NOT-IN-PAPER mechanism).

Usage: ``python -m benchmark.run_benchmark --help``
"""

from __future__ import annotations

import argparse
import importlib
import json
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

if __package__ in (None, ""):  # executed as a script: python3 benchmark/run_benchmark.py
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark.levels import (
    LEVELS,
    Environment,
    build_environment,
    disallowed_for,
    materialize,
    write_synthetic_inputs,
)
from benchmark.rubric import Finding, Rubric, Submission, default_rubric, grade_submission

# paper Methods "Fixed-input benchmark" p.38
ATTEMPTS_PER_CELL = 100  # 7 models x 5 levels x 100 = 3,500 attempts
MAX_OUTPUT_TOKENS = 1_000_000  # "maximum output budget of 1 million tokens"

Backend = Callable[["AttemptSpec"], dict[str, Any]]
Judge = Callable[[Submission, Rubric], Any]


@dataclass
class AttemptSpec:
    model: str
    level: str
    attempt: int  # 1-based index within the model x level cell
    run_dir: Path  # isolated output dir for this attempt
    prompt: str
    env: Environment
    max_output_tokens: int = MAX_OUTPUT_TOKENS
    seed: int = 0


def default_backend(model: str) -> Backend:
    """Live backend via ClaudeCodeBackend (owner-gated).

    Imports lazily so ``--help`` and offline use work without ``artharness``
    installed. Raises RuntimeError unless ARTHARNESS_ALLOW_LIVE=1 — the same
    gate as the campaign runner (README "Safety gates").
    """
    from artharness.config import CampaignConfig
    from artharness.runner.base import ClaudeCodeBackend, SessionSpec

    backend = ClaudeCodeBackend(model, CampaignConfig())  # raises unless gated env set

    def _run(spec: AttemptSpec) -> dict[str, Any]:
        backend.run(SessionSpec(
            role="worker",
            task_id=f"benchmark/{spec.model}/{spec.level}/{spec.attempt}",
            system_prompt="You are an autonomous research agent in the fixed-input benchmark.",
            user_prompt=spec.prompt,
            workdir=spec.run_dir,
            max_output_tokens=spec.max_output_tokens,
            disallowed_tools=disallowed_for(spec.env),  # enforce the L ladder (C1)
        ))
        # The session is instructed to write report.md + submission.json itself.
        # The RAW TEXT is returned: parsing happens in run_attempt so a
        # model-written invalid JSON is graded as a model zero, never raised
        # in here to masquerade as a backend (infra) error.
        report_path = spec.run_dir / "report.md"
        submission_path = spec.run_dir / "submission.json"
        return {
            "report": report_path.read_text() if report_path.exists() else "",
            "submission_text": submission_path.read_text()
            if submission_path.exists() else None,
        }

    return _run


def resolve_backend(dotted_path: str) -> Backend:
    """Load ``module.attr`` for --backend (callable AttemptSpec -> dict)."""
    module_name, _, attr = dotted_path.rpartition(".")
    if not module_name:
        raise ValueError(f"--backend must be a dotted path, got {dotted_path!r}")
    return getattr(importlib.import_module(module_name), attr)


def submission_from_dict(data: dict[str, Any],
                         report_path: str | None = None) -> Submission:
    findings = [
        Finding(
            claim=f["claim"],
            evidence=f.get("evidence", ""),
            confidence=f.get("confidence", 1.0),
            asserted=f.get("asserted"),
        )
        for f in data.get("findings", [])
    ]
    return Submission(findings=findings, report_path=report_path)


def run_attempt(spec: AttemptSpec, backend: Backend, judge: Judge,
                rubric: Rubric) -> dict[str, Any]:
    spec.run_dir.mkdir(parents=True, exist_ok=True)
    materialize(spec.env, spec.run_dir)  # C1: L3+ inputs must be IN the run dir
    base = {
        "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "model": spec.model,
        "level": spec.level,
        "attempt": spec.attempt,
        "run_dir": str(spec.run_dir),
        "max_output_tokens": spec.max_output_tokens,
    }
    try:
        out = backend(spec)
    except Exception as exc:  # infra failure must not grade as a model zero (C)
        base.update({"status": "backend_error", "error": str(exc)[:500],
                     "score": None, "recognized_repeat_array": None,
                     "asserted": None})
        return base
    report_path = spec.run_dir / "report.md"
    submission_path = spec.run_dir / "submission.json"
    report_path.write_text(out.get("report", ""))
    invalid_reason: str | None = None
    if "submission_text" in out:
        # session-written file: parse HERE. Invalid JSON is a MODEL failure
        # (audited 2026-09-28: it used to raise inside the backend and land as
        # backend_error/score=None — excluded from grading, i.e. treated
        # better than writing nothing, which scores 0).
        text = out["submission_text"]
        submission_data = {"findings": []}
        if text is not None:
            try:
                parsed = json.loads(text)
                if not isinstance(parsed, dict):
                    raise ValueError("submission.json is not a JSON object")
                submission_data = parsed
            except (json.JSONDecodeError, ValueError) as exc:
                submission_data = {"findings": []}
                invalid_reason = str(exc)[:200]
    else:
        submission_data = out.get("submission", {"findings": []})
    submission_path.write_text(json.dumps(submission_data, indent=2) + "\n")
    submission = submission_from_dict(submission_data, report_path=str(report_path))
    grade = judge(submission, rubric)
    base.update({
        "status": "invalid_submission" if invalid_reason else "ok",
        "parse_error": invalid_reason,
        "score": grade.score,
        "recognized_repeat_array": grade.recognized_repeat_array,
        "asserted": grade.asserted,
        "report_path": str(report_path),
        "submission_path": str(submission_path),
    })
    return base


def format_matrix(models: list[str], levels: list[str], attempts: int) -> str:
    lines = [f"benchmark matrix: {len(models)} models x {len(levels)} levels "
             f"x {attempts} attempts = {len(models) * len(levels) * attempts} total"]
    for m in models:
        for lv in levels:
            lines.append(f"  {m} {lv} x{attempts}")
    return "\n".join(lines)


def run_benchmark(models: list[str], levels: list[str], attempts: int,
                  outdir: Path, seed: int = 0, backend: Backend | None = None,
                  judge: Judge | None = None,
                  inputs_dir: Path | None = None,
                  allow_synthetic: bool = False) -> list[dict[str, Any]]:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rubric = default_rubric()
    judge = judge or grade_submission
    inputs_root = Path(inputs_dir) if inputs_dir else outdir / "inputs"
    envs = {}
    for level in levels:
        ldir = inputs_root / level
        if not (ldir / "loci").exists():
            if not allow_synthetic:
                # S1 finding C: silent fabrication graded synthetic fixtures as
                # if they were the paper's 96 loci. Opt in explicitly.
                raise FileNotFoundError(
                    f"no inputs at {ldir}/loci — pass --inputs with the real "
                    "benchmark inputs, or --synthetic to fabricate fixtures")
            write_synthetic_inputs(ldir)
        envs[level] = build_environment(level, ldir)
    records = []
    results_path = outdir / "results.jsonl"
    done: set[tuple[str, str, int]] = set()
    if results_path.exists():  # resume: (model, level, attempt) cells already on disk
        for line in results_path.read_text().splitlines():
            try:
                r = json.loads(line)
                done.add((r["model"], r["level"], r["attempt"]))
            except (json.JSONDecodeError, KeyError):
                continue
    with results_path.open("a") as fh:
        for model in models:
            be = backend if backend is not None else default_backend(model)
            for level in levels:
                env = envs[level]
                for i in range(1, attempts + 1):
                    if (model, level, i) in done:
                        continue  # append mode used to duplicate cells on re-run
                    spec = AttemptSpec(
                        model=model, level=level, attempt=i,
                        run_dir=outdir / "runs" / model / level / f"attempt_{i:03d}",
                        prompt=env.prompt, env=env,
                        max_output_tokens=MAX_OUTPUT_TOKENS, seed=seed,
                    )
                    rec = run_attempt(spec, be, judge, rubric)
                    fh.write(json.dumps(rec) + "\n")
                    records.append(rec)
    return records


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="run_benchmark",
        description="Fixed-input benchmark: models x levels x attempts "
                    "(paper Methods 'Fixed-input benchmark' p.38).")
    p.add_argument("--models", required=True,
                   help="comma-separated model names")
    p.add_argument("--levels", default=",".join(LEVELS),
                   help="comma-separated levels (default: L1..L5)")
    p.add_argument("--attempts", type=int, default=ATTEMPTS_PER_CELL,
                   help="attempts per model x level cell (default: 100, paper p.38)")
    p.add_argument("--outdir", default="benchmark-runs",
                   help="output dir for runs/ and results.jsonl")
    p.add_argument("--seed", type=int, default=0)  # NOT-IN-PAPER: paper states no seed
    p.add_argument("--backend", default=None,
                   help="dotted path to a callable AttemptSpec -> "
                        "{report, submission} (default: live ClaudeCodeBackend, "
                        "requires ARTHARNESS_ALLOW_LIVE=1)")
    p.add_argument("--inputs", default=None,
                   help="dir with real benchmark inputs (default: <outdir>/inputs)")
    p.add_argument("--synthetic", action="store_true",
                   help="explicitly fabricate synthetic fixtures when --inputs is "
                        "absent (S1: silent fabrication graded fake data)")
    p.add_argument("--dry-run", action="store_true",
                   help="print the attempt matrix without running")
    args = p.parse_args(argv)

    models = [m.strip() for m in args.models.split(",") if m.strip()]
    levels = [lv.strip() for lv in args.levels.split(",") if lv.strip()]
    bad = [lv for lv in levels if lv not in LEVELS]
    if bad:
        p.error(f"unknown levels {bad}; expected subset of {LEVELS}")

    if args.dry_run:
        print(format_matrix(models, levels, args.attempts))
        return 0

    backend = resolve_backend(args.backend) if args.backend else None
    records = run_benchmark(
        models=models, levels=levels, attempts=args.attempts,
        outdir=Path(args.outdir), seed=args.seed, backend=backend,
        inputs_dir=Path(args.inputs) if args.inputs else None,
        allow_synthetic=args.synthetic,
    )
    print(f"wrote {len(records)} attempt records to "
          f"{Path(args.outdir) / 'results.jsonl'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
