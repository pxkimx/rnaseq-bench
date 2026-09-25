"""An explicit min_genes=0 (or max_genes=0 / max_mt=0) QC threshold must be respected, not silently
replaced by the auto-computed one.

sc_pipeline.py picked its QC cut-offs with `int(P["min_genes"] or auto_min)` etc. In Python, `0 or
auto_min` evaluates to `auto_min` because 0 is falsy — so a user who deliberately asked for "no
minimum genes-per-cell floor" (min_genes=0, a real setting for data that is already filtered, or for
single-nucleus data where very low counts are still valid) silently got the MAD-based auto threshold
instead, with no warning, and cells the user explicitly asked to keep were removed.

Offline (RB_RESOURCE_DIR + no network needed); run: python tests/test_sc_zero_threshold.py
"""
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import numpy as np  # noqa: E402

from server import sc_pipeline, demo  # noqa: E402
from server.common import Job, run_safely  # noqa: E402


def check(cond, what):
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        sys.exit(1)


with tempfile.TemporaryDirectory() as td:
    root = Path(td)
    rng = np.random.default_rng(0)
    a = demo.simulate_pbmc()
    a = a[np.sort(rng.choice(a.n_obs, 700, replace=False))].copy()
    (root / "input").mkdir(parents=True)
    a.write_h5ad(root / "input" / "sim.h5ad")
    job = Job(root, "sc", "zero-min-genes-test")
    run_safely(job, sc_pipeline.run, [root / "input" / "sim.h5ad"], {"min_genes": 0, "doublets": False, "n_hvg": 400})
    st = json.loads((root / "status.json").read_text())
    check(st["state"] == "done", f"sc run with min_genes=0 completes ({st.get('error', '')})")
    res = json.loads((root / "result.json").read_text())
    check(res["params"]["min_genes"] == 0, f"min_genes=0 is respected, not replaced by the auto threshold (got {res['params']['min_genes']})")

print("all sc zero-threshold checks passed")
