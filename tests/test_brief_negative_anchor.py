"""Anchor 5 (paper Methods p.31, verbatim): "Neither this task brief nor the
research brief mentioned repeats or arrays." The reconstructed brief's body
must preserve that negative anchor — any mention of repeats/arrays leaks the
discovery the harness is supposed to find on its own. The provenance warning
block at the top is meta-commentary about the anchor and is exempt."""
import pathlib
import re

BRIEF = pathlib.Path(__file__).resolve().parents[1] / "brief" / "research_brief.md"
# stems, not whole words: "repetitive"/"repeated"/"tandemly" must not slip
# through (devin review F6); no legitimate brief term starts with these
LEAK = re.compile(r"\b(repeat|repetit|tandem|array)", re.IGNORECASE)


def _brief_body() -> str:
    text = BRIEF.read_text()
    # Drop the provenance warning block (blockquote after the title) — it names
    # the anchor itself, which is documentation, not brief content.
    return re.sub(r"^>.*?\n(?:>.*\n)*", "", text, count=1, flags=re.MULTILINE)


def test_brief_body_never_mentions_repeats_or_arrays():
    leaks = [(i + 1, line.strip()) for i, line in enumerate(_brief_body().splitlines())
             if LEAK.search(line)]
    assert not leaks, (
        "anchor 5 (negative) violated — brief body mentions repeats/arrays: "
        f"{leaks}; the discovery must not be spoiled (paper p.31: 'Neither "
        "this task brief nor the research brief mentioned repeats or arrays.')"
    )
