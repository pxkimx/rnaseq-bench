"""Quick end-to-end check of this installation:  python -m server.selftest
Runs a tiny simulated single-cell and bulk analysis (about a minute) and reports pass/fail."""
from __future__ import annotations

import json
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    from . import bulk_pipeline, demo, sc_pipeline
    from .common import Job, env_summary, run_safely

    print(env_summary())
    ok = True
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        # --- single-cell: 600 cells, 6 samples, planted IFN response
        a = demo.simulate_pbmc()
        rng = np.random.default_rng(0)
        a = a[np.sort(rng.choice(a.n_obs, 900, replace=False))].copy()
        d = root / "sc"; (d / "input").mkdir(parents=True)
        a.write_h5ad(d / "input" / "sim.h5ad")
        job = Job(d, "sc", "selftest-sc"); t = time.time()
        run_safely(job, sc_pipeline.run, [d / "input" / "sim.h5ad"], {"doublets": False, "n_hvg": 800})
        st = json.loads((d / "status.json").read_text())
        print(f"single-cell: {st['state']} in {time.time()-t:.0f}s", st.get("error", ""))
        if st["state"] != "done":
            ok = False; print(st.get("trace", "")[-2500:])
        else:
            r = json.loads((d / "result.json").read_text())
            print("  clusters:", r["summary"]["clusters"], "| sections:", ", ".join(s["id"] for s in r["sections"]))
        # --- bulk: example airway-like table
        fs, name = demo.airway(root / "bulk" / "input")
        d = root / "bulk"; job = Job(d, "bulk", "selftest-bulk"); t = time.time()
        run_safely(job, bulk_pipeline.run, fs, {})
        st = json.loads((d / "status.json").read_text())
        print(f"bulk: {st['state']} in {time.time()-t:.0f}s", st.get("error", ""))
        if st["state"] != "done":
            ok = False; print(st.get("trace", "")[-2500:])
        else:
            res = pd.read_csv(d / "deseq2_results.csv", index_col=0)
            print("  significant genes:", int((res.padj < 0.05).sum()))
    print("SELFTEST", "PASSED" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
