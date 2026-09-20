"""Shared job bookkeeping and figure helpers."""
from __future__ import annotations

import json
import time
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
    "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
    "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#393b79",
]
UP, DOWN, NS = "#c9304f", "#2a64b0", "#c7ccd1"


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    })


class Job:
    def __init__(self, root: Path, kind: str, name: str):
        self.root = root
        self.fig_dir = root / "figures"
        self.fig_dir.mkdir(parents=True, exist_ok=True)
        self.kind = kind
        self.result = {
            "kind": kind, "name": name, "tiles": [], "flags": [], "sections": [],
            "methods": [], "versions": {}, "downloads": [],
        }
        self._section = None
        self.set_status("running", 0, "Starting")

    # ---- progress ----
    def set_status(self, state, pct, msg, error=None, trace=None):
        s = {"state": state, "pct": pct, "message": msg, "time": time.time()}
        if error:
            s["error"] = error
        if trace:
            s["trace"] = trace
        (self.root / "status.json").write_text(json.dumps(s))

    def step(self, pct, msg):
        print(f"[{self.kind}] {pct:>3}% {msg}", flush=True)
        self.set_status("running", pct, msg)

    # ---- content ----
    def tile(self, value, label):
        self.result["tiles"].append({"value": value, "label": label})

    def flag(self, level, text):
        """level: ok | warn | info. text may contain <b>…</b>."""
        self.result["flags"].append({"level": level, "text": text})

    def section(self, sid, kicker, title, lede=""):
        self._section = {"id": sid, "kicker": kicker, "title": title, "lede": lede, "items": []}
        self.result["sections"].append(self._section)

    def figure(self, fid, title, fig=None, how="", yours="", wide=False, sub=""):
        fig = fig or plt.gcf()
        png = self.fig_dir / f"{fid}.png"
        fig.savefig(png, dpi=200, facecolor="white")
        plt.close("all")
        self._section["items"].append({
            "type": "figure", "id": fid, "title": title, "sub": sub,
            "src": f"figures/{fid}.png", "how": how, "yours": yours, "wide": wide,
        })

    def table(self, tid, title, columns, rows, note="", csv=None):
        item = {"type": "table", "id": tid, "title": title, "columns": columns,
                "rows": rows, "note": note}
        if csv:
            item["csv"] = csv
        self._section["items"].append(item)

    def widget(self, wtype, **kw):
        self._section["items"].append({"type": wtype, **kw})

    def method(self, head, text):
        self.result["methods"].append([head, text])

    def download(self, label, href):
        self.result["downloads"].append({"label": label, "href": href})

    def save(self):
        self.result["sections"] = [s for s in self.result["sections"] if s.get("items")]   # drop sections that ended up empty
        (self.root / "result.json").write_text(json.dumps(self.result, default=_json_default))


def _json_default(o):
    try:
        import numpy as np
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
    except Exception:
        pass
    return str(o)


def run_safely(job: Job, fn, *args, **kw):
    try:
        fn(job, *args, **kw)
        job.save()
        from .report import build_pdf
        job.step(98, "Writing PDF report")
        build_pdf(job.root)
        job.set_status("done", 100, "Analysis complete")
    except UserFacingError as e:
        job.set_status("error", 100, "Stopped", error=str(e))
    except Exception as e:  # pragma: no cover
        tb = traceback.format_exc()
        print(tb, flush=True)
        job.set_status("error", 100, "Stopped", error=f"{type(e).__name__}: {e}", trace=tb[-6000:] + "\n" + env_summary())


def log_exc(where: str) -> str:
    """Print the current exception's traceback to the server log and return a short 'Type: msg' string."""
    import sys
    et, ev, _ = sys.exc_info()
    print(f"[{where}] non-fatal error:\n" + traceback.format_exc(), flush=True)
    return f"{et.__name__ if et else 'Error'}: {str(ev)[:160]}"


def env_summary() -> str:
    import platform
    from importlib.metadata import version
    pk = ["numpy", "scipy", "pandas", "anndata", "scanpy", "matplotlib", "seaborn", "scikit-learn", "umap-learn",
          "leidenalg", "igraph", "harmonypy", "pydeseq2", "adjustText", "scikit-misc", "numba", "pynndescent"]
    vs = []
    for k in pk:
        try:
            vs.append(f"{k}={version(k)}")
        except Exception:  # noqa: BLE001
            vs.append(f"{k}=missing")
    return f"[env] python {platform.python_version()} {platform.machine()} | " + " ".join(vs)


class UserFacingError(Exception):
    pass


def fmt_n(v):
    v = float(v)
    a = abs(v)
    if a >= 1e6:
        return f"{v/1e6:.1f}M"
    if a >= 1e4:
        return f"{v/1e3:.0f}k"
    if a >= 1e3:
        return f"{v/1e3:.1f}k"
    return f"{v:.0f}" if a >= 10 or v == int(v) else f"{v:.2g}"


def commas(v):
    return f"{int(round(float(v))):,}"


def pct(a, b):
    return f"{100*a/b:.1f}%" if b else "0%"
