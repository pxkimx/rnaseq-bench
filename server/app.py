"""RNAseq Bench web server.  Run:  uvicorn server.app:app --port 8765"""
from __future__ import annotations

import json
import os
import threading
import uuid
from functools import lru_cache
from pathlib import Path

import numpy as np
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse, PlainTextResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import agent, bulk_pipeline, codegen, demo, dossier, geo, sc_pipeline
from .common import Cancelled, Job, run_safely
from .io_utils import guess_kind, unpack_archives

ROOT = Path(__file__).resolve().parent.parent
JOBS = Path(os.environ.get("TL_HOME", ROOT)) / "jobs"
JOBS.mkdir(parents=True, exist_ok=True)
VERSION = (ROOT / "VERSION").read_text().strip() if (ROOT / "VERSION").exists() else "dev"
app = FastAPI(title="RNAseq Bench", version=VERSION)
print(f"RNAseq Bench {VERSION} — jobs in {JOBS} — python {__import__('platform').python_version()}", flush=True)
_lock = threading.Semaphore(1)  # one heavy analysis at a time


def _start(kind: str, files: list[Path], name: str, params: dict, job_dir: Path):
    job = Job(job_dir, kind, name)
    (job_dir / "request.json").write_text(json.dumps({"kind": kind, "name": name, "params": params, "files": [str(f) for f in files]}))

    def work():
        with _lock:
            run_safely(job, sc_pipeline.run if kind == "sc" else bulk_pipeline.run, files, params)

    threading.Thread(target=work, daemon=True).start()
    return {"job": job_dir.name, "kind": kind}


@app.post("/api/analyze")
async def analyze(files: list[UploadFile] = File(...), kind: str = Form("auto"), params: str = Form("{}")):
    jid = uuid.uuid4().hex[:10]
    d = JOBS / jid
    (d / "input").mkdir(parents=True)
    saved = []
    for f in files:
        p = d / "input" / Path(f.filename).name
        with open(p, "wb") as out:
            while chunk := await f.read(1 << 20):
                out.write(chunk)
        saved.append(p)
    saved = unpack_archives(saved, d / "input")
    k = guess_kind(saved) if kind == "auto" else kind
    name = " + ".join(Path(f.filename).name for f in files)
    return _start(k, saved, name, json.loads(params or "{}"), d)


def _rerun(jid: str, params: dict):
    req = json.loads((JOBS / jid / "request.json").read_text())
    nid = uuid.uuid4().hex[:10]
    d = JOBS / nid
    d.mkdir(parents=True)
    return _start(req["kind"], [Path(p) for p in req["files"]], req["name"], {**req["params"], **params}, d)


@app.post("/api/rerun/{jid}")
async def rerun(jid: str, params: str = Form("{}")):
    return _rerun(jid, json.loads(params))


# ---------------------------------------------------------------- GEO import
class GeoList(BaseModel):
    accession: str


class GeoImport(BaseModel):
    accession: str
    urls: list[str]
    kind: str = "auto"
    params: dict = {}


@app.post("/api/geo/list")
def geo_list(body: GeoList):
    try:
        return geo.list_geo_files(body.accession)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(400, f"{type(e).__name__}: {e}")


@app.post("/api/geo/import")
def geo_import(body: GeoImport):
    jid = uuid.uuid4().hex[:10]
    d = JOBS / jid
    (d / "input").mkdir(parents=True)
    job = Job(d, "bulk", body.accession)
    job.set_status("running", 1, "Downloading from GEO…")

    def progress(i, n, name):
        job.check_cancelled()          # stop between files, so a large download can be abandoned
        job.set_status("running", 1 + int(8 * i / max(n, 1)), f"Downloading {name}")

    def work():
        try:
            files = geo.download_geo_files(body.urls, d / "input", progress=progress)
            files = unpack_archives(files, d / "input")
            k = guess_kind(files) if body.kind == "auto" else body.kind
            (d / "request.json").write_text(json.dumps({"kind": k, "name": body.accession.upper(), "params": body.params, "files": [str(f) for f in files]}))
            job2 = Job(d, k, body.accession.upper())
            with _lock:
                run_safely(job2, sc_pipeline.run if k == "sc" else bulk_pipeline.run, files, body.params)
        except Cancelled:
            job.set_status("cancelled", 100, "Stopped", error="You stopped this download.")
        except Exception as e:  # noqa: BLE001
            job.set_status("error", 100, "Stopped", error=f"{type(e).__name__}: {e}")

    threading.Thread(target=work, daemon=True).start()
    return {"job": jid, "kind": body.kind}


# ---------------------------------------------------------------- code & settings
@app.get("/api/jobs/{jid}/script", response_class=PlainTextResponse)
def job_script(jid: str):
    try:
        return codegen.script_for_job(JOBS / jid)
    except Exception as e:  # noqa: BLE001
        raise HTTPException(404, f"{type(e).__name__}: {e}")


@app.get("/api/source")
def source():
    return codegen.pipeline_sources()


class Settings(BaseModel):
    api_key: str | None = None
    model: str | None = None
    ncbi_email: str | None = None


@app.get("/api/settings")
def get_settings():
    cfg = geo.load_settings()
    key = cfg.get("api_key") or os.environ.get("ANTHROPIC_API_KEY") or ""
    return {"version": VERSION, "has_key": bool(key), "key_hint": (key[:7] + "…" + key[-4:]) if key else "", "model": cfg.get("model") or agent.DEFAULT_MODEL,
            "ncbi_email": cfg.get("ncbi_email") or ""}


@app.post("/api/settings")
def set_settings(body: Settings):
    d = {k: v for k, v in body.model_dump().items() if v}
    geo.save_settings(d)
    return get_settings()


@app.post("/api/quit")
def quit_server():
    """Stop the server (sidebar 'Quit'). Runs jobs are lost, finished ones stay on disk."""
    import threading as _t
    _t.Timer(0.5, lambda: os._exit(0)).start()
    return {"ok": True}


@app.get("/api/models")
def list_models():
    """Models the saved key can use (also verifies the key). Falls back to a static list."""
    return agent.available_models()


# ---------------------------------------------------------------- gene dossier
@app.get("/api/dossier/{gene}")
def gene_dossier(gene: str, species: str = "auto", context: str = "", refresh: bool = False):
    """NCBI Gene, UniProt, KEGG, AlphaFold DB and PubMed for one gene. Sources that cannot be reached come back
    as notes, so this only fails for a malformed request."""
    try:
        return dossier.build(gene, species, context, refresh)
    except ValueError as e:
        raise HTTPException(400, str(e))


# ---------------------------------------------------------------- agent
class ChatIn(BaseModel):
    conv: str
    message: str
    job: str | None = None


@app.post("/api/agent/chat")
def agent_chat(body: ChatIn):
    return StreamingResponse(agent.chat(body.conv, body.message, body.job, _rerun), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})


@app.post("/api/agent/reset")
def agent_reset(body: ChatIn):
    agent.reset(body.conv)
    return {"ok": True}


@app.post("/api/demo/{which}")
def run_demo(which: str):
    jid = uuid.uuid4().hex[:10]
    d = JOBS / jid
    if which == "sc":
        f, name = demo.pbmc3k(d / "input")
        return _start("sc", [f], name, {}, d)
    if which == "bulk":
        fs, name = demo.airway(d / "input")
        return _start("bulk", fs, name, {}, d)
    raise HTTPException(404)


@app.get("/api/jobs")
def list_jobs(limit: int = 60):
    """Every analysis on disk, newest first.

    The UI used to remember past runs in the browser's own storage, so anything started from another
    browser, from the assistant or from a re-run you had not opened yourself was simply unreachable —
    there was no way to get back to it without knowing its id.
    """
    rows = []
    if JOBS.exists():
        for d in sorted(JOBS.iterdir(), key=lambda x: x.stat().st_mtime, reverse=True):
            if not d.is_dir() or not (d / "request.json").exists():
                continue
            try:
                req = json.loads((d / "request.json").read_text())
                st = json.loads((d / "status.json").read_text()) if (d / "status.json").exists() else {}
                row = {"job": d.name, "kind": req.get("kind"), "name": req.get("name") or d.name,
                       "state": st.get("state"), "t": d.stat().st_mtime * 1000}
                if (d / "result.json").exists():
                    res = json.loads((d / "result.json").read_text())
                    row["summary"] = res.get("summary")
                elif st.get("state") in ("error", "cancelled"):
                    row["error"] = (st.get("error") or "")[:160]
                rows.append(row)
            except Exception:  # noqa: BLE001 - a half-written job must not break the list
                continue
            if len(rows) >= limit:
                break
    return rows


@app.post("/api/jobs/{jid}/cancel")
def cancel(jid: str):
    """Ask a running analysis to stop at its next step.

    A flag file rather than killing the worker: the pipeline checks it between steps, so it stops with
    its files consistent instead of leaving a half-written h5ad or PDF behind. A download in progress
    ends at the next step boundary too.
    """
    d = JOBS / jid
    if not (d / "status.json").exists():
        raise HTTPException(404, "Unknown job")
    st = json.loads((d / "status.json").read_text())
    if st.get("state") != "running":
        return {"ok": False, "state": st.get("state")}
    (d / "cancel").write_text("1")
    return {"ok": True}


@app.get("/api/jobs/{jid}/status")
def status(jid: str):
    p = JOBS / jid / "status.json"
    if not p.exists():
        raise HTTPException(404, "Unknown job")
    return JSONResponse(json.loads(p.read_text()))


@app.get("/api/jobs/{jid}/result")
def result(jid: str):
    p = JOBS / jid / "result.json"
    if not p.exists():
        raise HTTPException(404, "Result not ready")
    r = json.loads(p.read_text())
    r["job"] = jid
    r["report_name"] = report_name(jid)
    return JSONResponse(r)


@lru_cache(maxsize=2)
def _adata(jid: str):
    import scanpy as sc
    return sc.read_h5ad(JOBS / jid / "analyzed.h5ad")


class Selection(BaseModel):
    points: list[int]          # positions in the embedding arrays, i.e. what the user drew a loop around


@app.post("/api/jobs/{jid}/select_de")
def select_de(jid: str, body: Selection):
    """Rank the genes that separate a hand-drawn selection of cells from the rest.

    This is a ranking, not a test. The cells were selected by eye from a map built out of the same
    expression values, so any p-value here is optimistic by construction (the usual double-dipping
    problem), and cells from one sample are not independent replicates of each other. It is for finding
    what a population you can see is made of; a claim about a condition still needs the pseudobulk test.
    """
    import pandas as pd
    import scanpy as sc

    a = _adata(jid)
    emb = json.loads((JOBS / jid / "embedding.json").read_text())
    idx = np.asarray(emb["index"], dtype=int)
    pts = np.asarray(sorted(set(int(i) for i in body.points if 0 <= int(i) < len(idx))), dtype=int)
    if len(pts) < 20:
        raise HTTPException(400, "Select at least 20 cells — a smaller group gives nothing reliable.")
    mask = np.zeros(a.n_obs, dtype=bool)
    mask[idx[pts]] = True
    if int((~mask).sum()) < 20:
        raise HTTPException(400, "Almost every cell is selected, so there is nothing left to compare with.")

    lab = "cell_type_simple" if "cell_type_simple" in a.obs else "leiden"
    comp = (pd.Series(a.obs[lab].astype(str).values[mask]).value_counts().head(8) / int(mask.sum()) * 100)
    composition = [[k, round(float(v), 1)] for k, v in comp.items()]

    a.obs["_sel"] = pd.Categorical(np.where(mask, "selected", "rest"), categories=["rest", "selected"])
    sc.tl.rank_genes_groups(a, "_sel", groups=["selected"], reference="rest", method="wilcoxon")
    r = sc.get.rank_genes_groups_df(a, group="selected").dropna(subset=["logfoldchanges"])

    X = a[:, r["names"].tolist()].X
    dense = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    pct_in = (dense[mask] > 0).mean(0) * 100
    pct_out = (dense[~mask] > 0).mean(0) * 100
    r = r.assign(pct_in=pct_in, pct_out=pct_out)

    def rows(d):
        return [{"gene": x["names"], "lfc": round(float(x["logfoldchanges"]), 2),
                 "padj": float(x["pvals_adj"]), "pct_in": round(float(x["pct_in"]), 1),
                 "pct_out": round(float(x["pct_out"]), 1)} for _, x in d.iterrows()]

    up = r[r.logfoldchanges > 0].nlargest(25, "scores")
    down = r[r.logfoldchanges < 0].nsmallest(25, "scores")
    return {"n_selected": int(mask.sum()), "n_rest": int((~mask).sum()),
            "label_key": lab, "composition": composition,
            "up": rows(up), "down": rows(down)}


@app.get("/api/jobs/{jid}/gene/{gene}")
def gene(jid: str, gene: str):
    a = _adata(jid)
    emb = json.loads((JOBS / jid / "embedding.json").read_text())
    names = {g.upper(): g for g in a.var_names}
    g = names.get(gene.upper())
    if g is None:
        raise HTTPException(404, f"{gene} is not in this dataset")
    col = a[:, g].X
    v = np.asarray(col.toarray() if hasattr(col, "toarray") else col).ravel()[emb["index"]]
    cl = a.obs["leiden"].cat.codes.values
    means = [float(v[np.asarray(emb["cluster"]) == k].mean()) if (np.asarray(emb["cluster"]) == k).any() else 0 for k in range(len(emb["clusters"]))]
    return {"gene": g, "values": np.round(v, 3).tolist(), "cluster_means": means, "pct_expressing": float((v > 0).mean())}


def report_name(jid: str) -> str:
    """date_toolname_report — e.g. 2026-09-23_RNAseqBench-SingleCell_report.pdf. The date is the day the analysis
    ran (its request.json is written when it starts), so a report downloaded again later keeps the date of its contents."""
    import time as _time
    d = JOBS / jid
    try:
        req = d / "request.json"
        kind = json.loads(req.read_text()).get("kind")
        st = req.stat()
        t = getattr(st, "st_birthtime", st.st_mtime)   # creation time: the assistant can rewrite request.json later
    except Exception:  # noqa: BLE001
        kind, t = None, _time.time()
    tool = {"sc": "SingleCell", "bulk": "Bulk"}.get(kind, "Analysis")
    return f"{_time.strftime('%Y-%m-%d', _time.localtime(t))}_RNAseqBench-{tool}_report.pdf"


@app.get("/api/jobs/{jid}/files/{path:path}")
def files(jid: str, path: str):
    p = (JOBS / jid / path).resolve()
    if not str(p).startswith(str((JOBS / jid).resolve())) or not p.exists():
        raise HTTPException(404)
    dl = p.suffix in (".pdf", ".h5ad", ".csv")
    name = report_name(jid) if path == "report.pdf" else p.name
    return FileResponse(p, filename=name if dl else None)


app.mount("/", StaticFiles(directory=ROOT / "web", html=True), name="web")
