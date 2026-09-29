"""Two small "agents" that drive the same harness roles.

- :class:`RuleAgent`: no model. A fixed policy that reads the brief, calls the
  tools, writes plan/summary, reviews, curates. It exists so the whole workflow can be
  shown deterministically and offline, and so tests have an oracle.
- :class:`OllamaAgent`: a real small language model (default ``qwen3:0.6b``) behind an
  Ollama server, playing worker / supervisor / curator / editor through the same
  prompts the harness sends a big model. It calls tools with ``TOOL:`` lines, writes
  its own summary, and its verdicts are subject to the same scripted checks.

Both share :func:`evidence_problems`, the demo's version of the paper's supervisor
duty to "check quantitative claims against the artifacts": the summary must quote
each tool's ``HEADLINE`` line verbatim, and any multi-digit number it states must
appear in the tool output. A weak model that invents numbers gets a "revise" from the
script, not from another guess.
"""

from __future__ import annotations

import json
import os
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from ..accounting import SessionResult
from ..runner.base import BackendOutput, SessionSpec
from .tools import TOOLS, Sandbox, ToolError, ToolResult, run_tool

_NUM = re.compile(r"(?<![\w.])\d{2,}(?:\.\d+)?(?![\w.])")
# small models decorate lines with bullets and backticks; accept that, nothing more
_TOOL_LINE = re.compile(r"^[\s>*`-]*TOOL:\s*`?([A-Za-z_]+)\b([^`\n]*)`?\s*$", re.MULTILINE)
_INVENTED_RESULT = re.compile(r"^[\s>*`-]*\[?RESULT\b.*", re.MULTILINE | re.DOTALL)
_TASK = re.compile(r"# Task (t\d+) \(([^)]*)\)")
_LOCUS = re.compile(r"\blocus\s+(L\d+)\b", re.IGNORECASE)
_TOOLS_LINE = re.compile(r"^Tools:\s*(.+)$", re.MULTILINE)


def brief_tools(brief: str) -> list[str]:
    """Tools a brief names: an explicit `Tools:` line, else any tool it mentions."""
    m = _TOOLS_LINE.search(brief)
    if m:
        return [t.strip() for t in m.group(1).split(",") if t.strip() in TOOLS]
    return [t for t in TOOLS if t in brief] or (
        ["deep_dive"] if re.search(r"deep dive", brief, re.IGNORECASE) else [])


# -- tool log + evidence check --------------------------------------------------------

def log_path(workdir: Path) -> Path:
    return workdir / "artifacts" / "tool_log.jsonl"


def read_log(workdir: Path) -> list[dict]:
    p = log_path(workdir)
    if not p.exists():
        return []
    return [json.loads(ln) for ln in p.read_text().splitlines() if ln.strip()]


def exec_tool(sb: Sandbox, workdir: Path, name: str, args: str) -> tuple[str, ToolResult | None]:
    """Run a tool on behalf of a worker and append it to the task's tool log."""
    try:
        res = run_tool(sb, name, args)
        entry = {"tool": name, "args": args, "ok": True, "headline": res.headline,
                 "output": res.render()}
        shown = res.render()
    except ToolError as exc:
        res = None
        entry = {"tool": name, "args": args, "ok": False, "headline": "",
                 "output": f"ERROR: {exc}"}
        shown = entry["output"]
    with log_path(workdir).open("a") as fh:
        fh.write(json.dumps(entry) + "\n")
    return shown, res


def evidence_problems(summary: str, log: list[dict]) -> list[str]:
    """Scripted claim check. Empty list = the summary is backed by the tool log."""
    ok_calls = [e for e in log if e["ok"]]
    if not ok_calls:
        return ["no tool call succeeded in this task, so no claim in the summary has "
                "evidence; run the tools named in the brief"]
    problems = []
    for e in ok_calls:
        line = f"HEADLINE: {e['headline']}"
        if line not in summary:
            problems.append(f"the summary must quote the tool result verbatim; "
                            f"missing line: {line}")
    body = "\n".join(ln for ln in summary.splitlines() if not ln.startswith("HEADLINE:"))
    evidence = "\n".join(e["output"] for e in ok_calls)
    invented = sorted({n for n in _NUM.findall(body) if n not in evidence})
    if invented:
        problems.append("numbers stated without support in any tool output: "
                        + ", ".join(invented))
    return problems


def _result(spec: SessionSpec, t0: float, tin: int = 0, tout: int = 0) -> SessionResult:
    return SessionResult(role=spec.role, task_id=spec.task_id,
                         duration_s=round(time.monotonic() - t0, 3),
                         input_tokens_uncached=tin, output_tokens=tout, cache_write_tokens=0)


def _task_meta(spec: SessionSpec) -> tuple[str, str]:
    m = _TASK.search(spec.user_prompt)
    return (m.group(1), m.group(2)) if m else (spec.task_id or "", "")


def _attempt(spec: SessionSpec) -> int:
    return len(list(spec.workdir.glob("verdict-r*.md")))


def _brief_of(spec: SessionSpec) -> str:
    p = spec.workdir / "brief.md"
    return p.read_text() if p.exists() else spec.user_prompt


class DemoAgent:
    """Shared plumbing: sandbox wiring, role dispatch, files a worker must leave."""

    def __init__(self, genomes: Path, shared: Path):
        self.genomes, self.shared = Path(genomes), Path(shared)

    def sandbox(self, workdir: Path) -> Sandbox:
        (workdir / "artifacts").mkdir(parents=True, exist_ok=True)
        return Sandbox(self.genomes, self.shared, workdir / "artifacts")

    def run(self, spec: SessionSpec) -> BackendOutput:
        return getattr(self, f"_{spec.role}")(spec)

    def _worker(self, spec):  # pragma: no cover - abstract
        raise NotImplementedError

    def _supervisor(self, spec):  # pragma: no cover - abstract
        raise NotImplementedError

    def _curator(self, spec):  # pragma: no cover - abstract
        raise NotImplementedError

    def _editor(self, spec):  # pragma: no cover - abstract
        raise NotImplementedError

    # scripted review shared by both agents ----------------------------------------
    def evidence_review(self, spec: SessionSpec) -> list[str]:
        summary_p = spec.workdir / "summary.md"
        if not summary_p.exists():
            return ["no summary.md was submitted"]
        return evidence_problems(summary_p.read_text(), read_log(spec.workdir))


# -- rule agent (no model) ---------------------------------------------------------------

class RuleAgent(DemoAgent):
    """Deterministic policy. ``stage_revision`` makes the first pass of the stage-4
    task sloppy (summary from memory, no tool run) so the revise loop is on screen; it
    is a staged demonstration, labelled as such in the task record."""

    def __init__(self, genomes: Path, shared: Path, stage_revision: bool = True):
        super().__init__(genomes, shared)
        self.stage_revision = stage_revision

    def _worker(self, spec: SessionSpec) -> BackendOutput:
        t0 = time.monotonic()
        tid, stage = _task_meta(spec)
        brief = _brief_of(spec)
        sb = self.sandbox(spec.workdir)
        first_pass = not list(spec.workdir.glob("verdict-r*.md"))
        if self.stage_revision and stage.startswith("4_") and first_pass \
                and brief.startswith("Stage 4"):
            (spec.workdir / "plan.md").write_text("# plan\nScan upstream windows.\n")
            (spec.workdir / "summary.md").write_text(
                "# summary\nAbout 4 arrays were found among the loci; the census is "
                "complete. (Staged sloppy pass: written from memory, no tool was run, "
                "to show the supervisor's check.)\n")
            return BackendOutput(result=_result(spec, t0))
        calls: list[tuple[str, str]] = []
        for tool in brief_tools(brief):
            m = _LOCUS.search(brief)
            calls.append((tool, f"locus={m.group(1)}" if tool == "deep_dive" and m else ""))
        (spec.workdir / "artifacts" / "tool_log.jsonl").unlink(missing_ok=True)
        heads, sugg = [], []
        for tool, args in calls:
            _, res = exec_tool(sb, spec.workdir, tool, args)
            if res:
                heads.append(f"HEADLINE: {res.headline}")
                sugg += res.suggestions
        (spec.workdir / "plan.md").write_text(
            "# plan\n" + "".join(f"{i}. run `{t} {a}`.\n" for i, (t, a) in enumerate(calls, 1))
            + "Then quote every HEADLINE and pass suggestions on as follow-ups.\n")
        (spec.workdir / "summary.md").write_text(
            f"# summary ({tid})\nRan {', '.join(t for t, _ in calls)}; tables are in "
            "artifacts/ and the shared sandbox.\n" + "\n".join(heads) + "\n")
        return BackendOutput(result=_result(spec, t0), proposed_followups=sugg)

    def _supervisor(self, spec: SessionSpec) -> BackendOutput:
        t0 = time.monotonic()
        problems = self.evidence_review(spec)
        if problems:
            return BackendOutput(result=_result(spec, t0), verdict="revise",
                                 verdict_notes="\n".join(f"- {p}" for p in problems))
        proposed = []
        log = read_log(spec.workdir)
        for e in log:
            if e["tool"] == "scan_arrays" and "not_assessed=0" not in e["headline"]:
                loci = [ln.split("\t")[0] for ln in e["output"].splitlines()
                        if "\tnot_assessed" in ln]
                proposed += [f"Deep dive locus {lid}: the scan could not assess it, so "
                             "look for an array anyway." for lid in loci]
        return BackendOutput(result=_result(spec, t0), verdict="accept",
                             verdict_notes="every HEADLINE is quoted and no unsupported "
                                           "number appears",
                             proposed_followups=proposed)

    def _curator(self, spec: SessionSpec) -> BackendOutput:
        t0 = time.monotonic()
        lines = [ln for ln in (spec.workdir / "summary.md").read_text().splitlines()
                 if ln.startswith("HEADLINE:")]
        text = "Findings (rule agent copies the verified headlines):\n" + "\n".join(
            f"- {ln[len('HEADLINE: '):]}" for ln in lines)
        return BackendOutput(result=_result(spec, t0), text=text)

    def _editor(self, spec: SessionSpec) -> BackendOutput:
        t0 = time.monotonic()
        problems = report_problems(spec.user_prompt, self.shared)
        return BackendOutput(
            result=_result(spec, t0), verdict="no" if problems else "yes",
            verdict_notes="\n".join(f"- {p}" for p in problems) or "claims match the tables")


def report_problems(report: str, shared: Path) -> list[str]:
    """Scripted editor check: every locus a report names is a known RT candidate and
    the array count it states matches arrays.tsv."""
    arrays_p, cands_p = shared / "arrays.tsv", shared / "rt_candidates.tsv"
    if not arrays_p.exists() or not cands_p.exists():
        return ["arrays.tsv or rt_candidates.tsv is missing; the report has nothing "
                "to stand on"]
    rows = [ln.split("\t") for ln in arrays_p.read_text().splitlines()[1:] if ln]
    status = {r[0]: r[1] for r in rows}
    known = {ln.split("\t")[0] for ln in cands_p.read_text().splitlines()[1:] if ln}
    problems = []
    for lid in sorted(set(re.findall(r"\bL\d{2}\b", report))):
        if lid not in known:
            problems.append(f"{lid} appears in the report but is not an RT candidate")
    m = re.search(r"(\d+) array-positive", report)
    n = sum(1 for s in status.values() if s == "array")
    if not m or int(m.group(1)) != n:
        problems.append(f"the report must state '{n} array-positive' loci as arrays.tsv does")
    return problems


# -- small-LLM agent ----------------------------------------------------------------------------

PROTOCOL = """\
How to work (plain text only):
1. First reply: two short lines starting with PLAN:, then one line per tool call,
   written exactly like `TOOL: scan_arrays` or `TOOL: deep_dive locus=L01`.
2. The system runs your tools and shows RESULT blocks.
3. Second reply: `SUMMARY:`, then copy every `HEADLINE:` line from the RESULT blocks
   unchanged, then ONE plain sentence with no digits. Do not add any other number.
4. Only if a RESULT block contains `SUGGEST_FOLLOWUP:` lines, add one
   `PROPOSE_FOLLOWUP: <same text>` line for each. Otherwise write no such line.
5. Use only the tools the brief names; anything else belongs in a follow-up.
Example of a correct second reply (made-up tool):
SUMMARY:
HEADLINE: widgets=3 sprockets=1
Every widget was found and a single sprocket was set aside.
Tools:
"""


@dataclass
class ChatReply:
    text: str
    tokens_in: int
    tokens_out: int


class OllamaChat:
    """Minimal Ollama /api/chat client (stdlib only)."""

    def __init__(self, model: str, base_url: str, timeout: float = 300.0, seed: int = 7):
        self.model, self.base_url, self.timeout, self.seed = model, base_url.rstrip("/"), \
            timeout, seed
        self._think: bool | None = False

    def chat(self, messages: list[dict], attempt: int = 0) -> ChatReply:
        # temperature 0 makes a retry repeat itself; the attempt number moves the seed
        body = {"model": self.model, "messages": messages, "stream": False,
                "options": {"temperature": 0, "seed": self.seed + attempt, "num_ctx": 6144,
                            "num_predict": 700}}
        if self._think is not None:
            body["think"] = self._think
        req = urllib.request.Request(f"{self.base_url}/api/chat",
                                     data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as r:
                data = json.loads(r.read())
        except urllib.error.HTTPError as exc:
            if exc.code == 400 and self._think is not None:
                self._think = None  # model has no thinking switch
                return self.chat(messages, attempt)
            raise
        text = data["message"]["content"]
        text = re.sub(r"<think>.*?</think>", "", text, flags=re.DOTALL).strip()
        return ChatReply(text, int(data.get("prompt_eval_count", 0)),
                         int(data.get("eval_count", 0)))


def check_ollama(model: str, base_url: str) -> None:
    """Fail early and clearly when the server or model is missing. A non-local server
    is a paid/remote model call, so it needs the same owner gate as the Claude backend."""
    host = re.sub(r"^https?://", "", base_url).split(":")[0].split("/")[0]
    if host not in ("localhost", "127.0.0.1", "::1") \
            and os.environ.get("ARTHARNESS_ALLOW_LIVE") != "1":
        raise RuntimeError("a non-local model server needs ARTHARNESS_ALLOW_LIVE=1 "
                           "(owner gate; see README Safety gates)")
    try:
        with urllib.request.urlopen(f"{base_url.rstrip('/')}/api/tags", timeout=5) as r:
            names = {m["name"] for m in json.loads(r.read()).get("models", [])}
    except (urllib.error.URLError, OSError) as exc:
        raise RuntimeError(f"cannot reach the model server at {base_url}: {exc}. "
                           "Start it (`ollama serve`) or use --agent rules.") from exc
    if model not in names and f"{model}:latest" not in names:
        raise RuntimeError(f"model {model!r} is not installed; `ollama pull {model}` "
                           f"(installed: {', '.join(sorted(names)) or 'none'})")


class OllamaAgent(DemoAgent):
    def __init__(self, genomes: Path, shared: Path, chat: OllamaChat, max_tool_calls: int = 3):
        super().__init__(genomes, shared)
        self.llm = chat
        self.max_tool_calls = max_tool_calls

    def _sys(self, spec: SessionSpec, with_protocol: bool = False) -> str:
        tools = "\n".join(f"- {desc}" for _, desc in TOOLS.values())
        return spec.system_prompt + ("\n" + PROTOCOL + tools if with_protocol else "")

    def _worker(self, spec: SessionSpec) -> BackendOutput:
        t0 = time.monotonic()
        sb = self.sandbox(spec.workdir)
        brief = _brief_of(spec)
        named = brief_tools(brief)
        locus = _LOCUS.search(brief)
        call = named[0] + (f" locus={locus.group(1)}" if named and named[0] == "deep_dive"
                           and locus else "") if named else ""
        hint = (f"\nThe brief names: {', '.join(named)}. Your first reply is `PLAN:` "
                f"lines, then the line `TOOL: {call}`." if named else "")
        msgs = [{"role": "system", "content": self._sys(spec, True)},
                {"role": "user", "content": spec.user_prompt + hint}]
        tin = tout = 0
        prior = [e for e in read_log(spec.workdir) if e["ok"]]
        if prior:
            # A revision keeps the tool results it already has; only the write-up is
            # in question, so the model is asked to fix its summary, not to rerun.
            blocks = [f"RESULT of `{e['tool']} {e['args']}`:\n{e['output'][:1800]}"
                      for e in prior]
        else:
            r1 = self.llm.chat(msgs, _attempt(spec))
            tin, tout = r1.tokens_in, r1.tokens_out
            # a model that writes its own RESULT blocks is inventing them: cut there
            r1.text = _INVENTED_RESULT.sub("", r1.text).strip()
            plan = "\n".join(ln for ln in r1.text.splitlines()
                             if not _TOOL_LINE.match(ln))
            (spec.workdir / "plan.md").write_text(f"# plan (model)\n{plan.strip()}\n")
            seen, blocks = set(), []
            for name, args in _TOOL_LINE.findall(r1.text):
                if (name, args.strip()) in seen or len(seen) >= self.max_tool_calls:
                    continue
                seen.add((name, args.strip()))
                if named and name not in named:
                    shown = (f"ERROR: {name} is outside this task's scope (the brief names "
                             f"{', '.join(named)}); propose a follow-up instead")
                    with log_path(spec.workdir).open("a") as fh:
                        fh.write(json.dumps({"tool": name, "args": args.strip(), "ok": False,
                                             "headline": "", "output": shown}) + "\n")
                else:
                    shown, _ = exec_tool(sb, spec.workdir, name, args.strip())
                blocks.append(f"RESULT of `{name} {args.strip()}`:\n{shown[:1800]}")
            if not blocks:
                (spec.workdir / "summary.md").write_text(
                    "# summary (model)\n" + r1.text.strip() + "\n")
                return BackendOutput(result=_result(spec, t0, tin, tout))
            msgs.append({"role": "assistant", "content": r1.text})
        msgs.append({"role": "user", "content": "\n\n".join(blocks)
                     + "\n\nNow write SUMMARY: (copy each HEADLINE line unchanged) and any "
                       "PROPOSE_FOLLOWUP lines."})
        r2 = self.llm.chat(msgs, _attempt(spec))
        tin, tout = tin + r2.tokens_in, tout + r2.tokens_out
        summary = re.sub(r"^\s*SUMMARY:\s*", "", r2.text, flags=re.MULTILINE)
        summary = "\n".join(ln for ln in summary.splitlines()
                            if not ln.strip().startswith("PROPOSE_FOLLOWUP:"))
        (spec.workdir / "summary.md").write_text(f"# summary (model)\n{summary.strip()}\n")
        follow = [re.sub(r"^(SUGGEST_FOLLOWUP:|PROPOSE_FOLLOWUP:|[-*\s])+", "", f).strip()
                  for f in re.findall(r"PROPOSE_FOLLOWUP:\s*(.+)", r2.text)]
        follow = [f for f in follow if len(f) > 12]
        return BackendOutput(result=_result(spec, t0, tin, tout), proposed_followups=follow)

    def _supervisor(self, spec: SessionSpec) -> BackendOutput:
        t0 = time.monotonic()
        problems = self.evidence_review(spec)
        checks = ("SCRIPTED CHECKS: PASS" if not problems else
                  "SCRIPTED CHECKS: FAIL\n" + "\n".join(f"- {p}" for p in problems))
        r = self.llm.chat([
            {"role": "system", "content": self._sys(spec)
             + "\nEnd with exactly one line `VERDICT: accept` or `VERDICT: revise`. "
               "If the scripted checks PASS and you see no concrete error, accept."},
            {"role": "user", "content": spec.user_prompt + "\n\n" + checks}],
            _attempt(spec))
        m = list(re.finditer(r"VERDICT:\s*(accept|revise)\b", r.text, re.IGNORECASE))
        model_verdict = m[-1].group(1).lower() if m else None
        # the script is a necessary condition; the model's verdict is the rest
        verdict = "revise" if problems else model_verdict
        notes = ("\n".join(f"- {p}" for p in problems) + "\n" if problems else "") \
            + r.text.strip()[:1200]
        proposed = re.findall(r"PROPOSED_TASK_BRIEF:\s*(.+?)(?=\n\n|\Z)", r.text, re.DOTALL)
        return BackendOutput(result=_result(spec, t0, r.tokens_in, r.tokens_out),
                             verdict=verdict, verdict_notes=notes[:2000] or None,
                             proposed_followups=[p.strip() for p in proposed])

    def _curator(self, spec: SessionSpec) -> BackendOutput:
        t0 = time.monotonic()
        r = self.llm.chat([
            {"role": "system", "content": self._sys(spec)
             + "\nWrite 2-3 factual sentences. Keep every number exactly as in the summary."},
            {"role": "user", "content": spec.user_prompt}], _attempt(spec))
        return BackendOutput(result=_result(spec, t0, r.tokens_in, r.tokens_out), text=r.text)

    def _editor(self, spec: SessionSpec) -> BackendOutput:
        t0 = time.monotonic()
        problems = report_problems(spec.user_prompt, self.shared)
        r = self.llm.chat([
            {"role": "system", "content": self._sys(spec)},
            {"role": "user", "content": spec.user_prompt + (
                "\n\nSCRIPTED CHECKS: FAIL\n" + "\n".join(f"- {p}" for p in problems)
                if problems else "\n\nSCRIPTED CHECKS: PASS")}])
        m = list(re.finditer(r"FILE:\s*(yes|no)\b", r.text, re.IGNORECASE))
        model_verdict = m[-1].group(1).lower() if m else None
        verdict = "no" if problems else model_verdict
        return BackendOutput(result=_result(spec, t0, r.tokens_in, r.tokens_out),
                             verdict=verdict, verdict_notes=r.text.strip()[:1500])
