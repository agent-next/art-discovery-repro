"""replay.html: a 2D animation of the run, built from its recorded events.

The page is one self-contained file (no network, no libraries). Each animation step is
one event the harness actually emitted (task opened, session finished, gate checked,
follow-up triaged), so the picture cannot drift from what happened.
"""

from __future__ import annotations

import json
from pathlib import Path

TEMPLATE = Path(__file__).with_name("replay_template.html")


def write_replay(out: Path, events: list[dict], meta: dict) -> Path:
    data = json.dumps({"events": events, "meta": meta}, ensure_ascii=False)
    # "</" inside a JSON script block would end it early
    data = data.replace("</", "<\\/")
    page = TEMPLATE.read_text().replace("__DATA__", data)
    (out / "events.json").write_text(json.dumps(events, indent=1) + "\n")
    path = out / "replay.html"
    path.write_text(page)
    return path
