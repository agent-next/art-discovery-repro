"""Live commentary: what is happening, and why the paper's harness does it.

Wraps the session backend, the stage gates and the triage policy so every harness
mechanism announces itself once with its paper citation. Purely observational: the
wrapped objects behave exactly as before.
"""

from __future__ import annotations

import json
import threading
from collections.abc import Callable
from pathlib import Path

from ..runner.base import BackendOutput, SessionSpec

WHY = {
    "worker": "worker = the agent that proposes a plan and executes it, then submits a "
              "written summary with its files (paper p.28).",
    "supervisor": "supervisor = a second agent that reviews plan, summary and files and "
                  "either accepts or returns them for revision; 49 of the paper's 119 "
                  "tasks were revised at least once.",
    "curator": "curator = after each accepted task it enters the findings into a shared "
               "knowledge base that later workers and supervisors read (paper p.28).",
    "editor": "editor = reviews a report before it is filed; the paper used 52 editor "
              "sessions for 19 reports.",
    "gate": "gate = a SCRIPTED completion check that closes a stage. No task of a later "
            "stage opens before it passes; the check is code, not a model's opinion "
            "(paper p.28).",
    "triage": "triage = follow-up tasks proposed by agents queue here; the harness "
              "releases each or rejects it with a written reason (10 rejections in the "
              "paper).",
    "revise": "revision = the supervisor's 'no'. The worker must pass the completion "
              "check again; after 10 revisions a task stalls in the paper (3 here).",
}


class Narrator:
    def __init__(self, say: Callable[[str], None] = print, explain: bool = True):
        self.say, self.explain = say, explain
        self._told: set[str] = set()
        self._lock = threading.Lock()

    def note(self, key: str) -> None:
        if self.explain and key not in self._told:
            self._told.add(key)
            self.say(f"      why: {WHY[key]}")

    def emit(self, line: str, *why: str) -> None:
        with self._lock:
            self.say(line)
            for key in why:
                self.note(key)


class NarratedBackend:
    def __init__(self, inner, narrator: Narrator, records: Path):
        self.inner, self.n, self.records = inner, narrator, records
        self._opened: set[str] = set()

    def _meta(self, tid: str | None) -> dict:
        p = self.records / (tid or "") / "meta.json"
        return json.loads(p.read_text()) if tid and p.exists() else {}

    def run(self, spec: SessionSpec) -> BackendOutput:
        meta = self._meta(spec.task_id)
        tid = spec.task_id or "-"
        if spec.role == "worker" and tid not in self._opened:
            self._opened.add(tid)
            brief = (spec.workdir / "brief.md").read_text().splitlines()[0][:96]
            origin = meta.get("origin", "?") + (f" of {meta['parent']}" if meta.get("parent")
                                                else "")
            self.n.emit(f"  task {tid} [{meta.get('stage', '?')}] ({origin}): {brief}")
        out = self.inner.run(spec)
        r = out.result
        tok = f"{r.input_tokens_uncached}->{r.output_tokens} tok" if (
            r.input_tokens_uncached or r.output_tokens) else "no model"
        detail = ""
        if spec.role == "worker":
            log = spec.workdir / "artifacts" / "tool_log.jsonl"
            calls = [json.loads(ln) for ln in log.read_text().splitlines()] if log.exists() \
                else []
            detail = ("tools: " + "; ".join(
                f"{c['tool']} -> {c['headline'] or c['output'][:50]}" for c in calls)
                if calls else "ran NO tool")
            if out.proposed_followups:
                detail += f"; proposes {len(out.proposed_followups)} follow-up(s)"
        elif spec.role == "supervisor":
            detail = f"VERDICT: {out.verdict}"
            if out.verdict != "accept":
                first = (out.verdict_notes or "no reason given").strip().splitlines()[0]
                detail += f" -- {first[:110]}"
        elif spec.role == "editor":
            detail = f"FILE: {out.verdict}"
        elif spec.role == "curator":
            detail = "entered findings into the knowledge base"
        self.n.emit(f"    {spec.role:<10} {tid} {r.duration_s:5.1f}s {tok:<14} {detail}",
                    spec.role, *(["revise"] if spec.role == "supervisor"
                                 and out.verdict != "accept" else []))
        return out
