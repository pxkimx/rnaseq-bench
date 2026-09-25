"""An explicit lfc=0 (or any other falsy-but-valid threshold) must be respected, not silently
replaced by the built-in default.

bulk_pipeline.py picked its LFC and FDR thresholds with `float(params.get("x") or default)`. In
Python, `0 or default` evaluates to `default` because 0 is falsy — so a user who deliberately asked
for "no minimum fold-change cut-off" (lfc=0, a real, meaningful setting: call everything significant
regardless of effect size) silently got the 1.0 default instead, with no warning, and the report's
Parameters section then claimed a threshold that was never actually applied to the gene list.

Offline (RB_RESOURCE_DIR + no network needed); run: python tests/test_bulk_zero_threshold.py
"""
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from server import bulk_pipeline, demo  # noqa: E402
from server.common import Job, run_safely  # noqa: E402


def check(cond, what):
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        sys.exit(1)


with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    fs, _ = demo.airway(root / "input")
    job = Job(root, "bulk", "zero-lfc-test")
    run_safely(job, bulk_pipeline.run, fs, {"lfc": 0, "gsea": False, "activity": False, "enrichment": False})
    st = json.loads((root / "status.json").read_text())
    check(st["state"] == "done", f"bulk run with lfc=0 completes ({st.get('error', '')})")
    res = json.loads((root / "result.json").read_text())
    check(res["params"]["lfc"] == 0, f"lfc=0 is respected, not replaced by the 1.0 default (got {res['params']['lfc']})")

# sanity check: a non-default, non-zero value still passes through as before (guards against
# over-fixing this into something that ignores the default entirely)
with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    fs, _ = demo.airway(root / "input")
    job = Job(root, "bulk", "custom-alpha-test")
    run_safely(job, bulk_pipeline.run, fs, {"alpha": 0.2, "gsea": False, "activity": False, "enrichment": False})
    res = json.loads((root / "result.json").read_text())
    check(res["params"]["alpha"] == 0.2, f"a non-default alpha is respected (got {res['params']['alpha']})")

print("all bulk zero-threshold checks passed")
