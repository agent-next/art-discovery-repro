# 20260929-demo: recorded runs of `python -m artharness.demo`

Commands (from the repo root at the commit that added the demo):

```
python -m artharness.demo --agent rules  --out /tmp/dm/final-rules  --force
python -m artharness.demo --agent ollama --out /tmp/dm/final-ollama --force   # qwen3:0.6b, local Ollama
```

| file | what |
| --- | --- |
| `console-rules.txt`, `WALKTHROUGH-rules.md` | rule policy: 9 tasks, 1 staged revision, 1 triage rejection, report filed, 4/4 arrays, 0 false calls |
| `console-ollama.txt`, `WALKTHROUGH-ollama.md` | qwen3:0.6b: 6 tasks, 5 revised, 4 triage rejections, 0 stalled, report filed, 4/4 arrays, 0 false calls |

The world is 9 synthetic contigs; nothing here supports a biological claim. Model output
varies with seed and prompt; the small-model run is a recording, not a guarantee.
