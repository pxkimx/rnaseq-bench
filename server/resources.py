"""Reference resources (gene-set libraries, pathway / TF networks, CellTypist models).

Everything is downloaded once and cached in ~/.rnaseq-bench so later runs work offline.
Every accessor returns None (never raises) when the resource cannot be obtained, and the
pipelines then skip that analysis with an explanatory note.

RB_RESOURCE_DIR: optional folder with pre-made files (used by the test-suite to run offline):
    genesets/<library>.json          {"term": ["GENE", ...]}
    progeny_<organism>.csv           columns source, target, weight
    collectri_<organism>.csv         columns source, target, weight
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pandas as pd

from .geo import CACHE

RES = Path(os.environ.get("RB_RESOURCE_DIR") or (CACHE / "resources"))
GS_DIR = RES / "genesets"

# Enrichr library names, per organism (the mouse KEGG library is species-specific; Hallmark/GO are human
# symbols and are used for mouse after upper-casing gene names — standard practice, imperfect for the
# handful of genes without a 1:1 orthologue).
LIBRARIES = {
    "human": ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"],
    "mouse": ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "KEGG_2019_Mouse", "Reactome_2022"],
}
SHORT = {"MSigDB_Hallmark_2020": "Hallmark", "GO_Biological_Process_2023": "GO BP", "KEGG_2021_Human": "KEGG",
         "KEGG_2019_Mouse": "KEGG", "Reactome_2022": "Reactome"}


def _log(msg):
    print(f"[resources] {msg}", flush=True)


def gene_sets(library: str, organism: str = "human") -> dict[str, list[str]] | None:
    """One Enrichr library as {term: [genes]}, cached as JSON."""
    GS_DIR.mkdir(parents=True, exist_ok=True)
    f = GS_DIR / f"{library}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except Exception:  # noqa: BLE001
            f.unlink(missing_ok=True)
    try:
        import gseapy as gp
        lib = gp.get_library(name=library, organism="Human" if organism == "human" else "Mouse")
        lib = {k: sorted(set(map(str, v))) for k, v in lib.items()}
        f.write_text(json.dumps(lib))
        _log(f"downloaded {library} ({len(lib)} sets)")
        return lib
    except Exception as e:  # noqa: BLE001
        _log(f"{library} unavailable: {type(e).__name__}: {str(e)[:120]}")
        return None


def all_gene_sets(organism: str) -> dict[str, dict[str, list[str]]]:
    """{short name: library} for every library that could be obtained."""
    out = {}
    for lib in LIBRARIES.get(organism, LIBRARIES["human"]):
        gs = gene_sets(lib, organism)
        if gs:
            out[SHORT.get(lib, lib)] = gs
    return out


def _net(kind: str, organism: str) -> pd.DataFrame | None:
    RES.mkdir(parents=True, exist_ok=True)
    f = RES / f"{kind}_{organism}.csv"
    if f.exists():
        try:
            return pd.read_csv(f)
        except Exception:  # noqa: BLE001
            f.unlink(missing_ok=True)
    try:
        import decoupler as dc
        if kind == "progeny":
            net = dc.op.progeny(organism=organism, top=500)
        elif kind == "collectri":
            net = dc.op.collectri(organism=organism)
        elif kind == "hallmark":
            net = dc.op.hallmark(organism=organism)
        else:
            return None
        net = net[["source", "target"] + (["weight"] if "weight" in net.columns else [])].copy()
        if "weight" not in net.columns:
            net["weight"] = 1.0
        net.to_csv(f, index=False)
        _log(f"downloaded {kind} ({organism}, {len(net)} edges)")
        return net
    except Exception as e:  # noqa: BLE001
        _log(f"{kind} ({organism}) unavailable: {type(e).__name__}: {str(e)[:120]}")
        return None


def progeny(organism: str = "human") -> pd.DataFrame | None:
    """PROGENy: 14 signalling pathways -> responsive genes with weights (top 500 per pathway)."""
    return _net("progeny", organism)


def collectri(organism: str = "human") -> pd.DataFrame | None:
    """CollecTRI: transcription factor -> target regulons with sign."""
    return _net("collectri", organism)


def celltypist_model(name: str) -> str | None:
    """Path to a CellTypist model, downloading it on first use."""
    try:
        import celltypist
        from celltypist import models
        models.models_path = str(RES / "celltypist")
        Path(models.models_path).mkdir(parents=True, exist_ok=True)
        p = Path(models.models_path) / name
        if not p.exists():
            models.download_models(model=name, force_update=False)
        return str(p) if p.exists() else None
    except Exception as e:  # noqa: BLE001
        _log(f"CellTypist model {name} unavailable: {type(e).__name__}: {str(e)[:120]}")
        return None


def status() -> dict:
    """What is cached (for the Settings dialog)."""
    have = {p.stem: int(p.stat().st_size // 1024) for p in GS_DIR.glob("*.json")} if GS_DIR.exists() else {}
    nets = [p.name for p in RES.glob("*.csv")] if RES.exists() else []
    ct = [p.name for p in (RES / "celltypist").glob("*.pkl")] if (RES / "celltypist").exists() else []
    return {"dir": str(RES), "gene_set_libraries": have, "networks": nets, "celltypist_models": ct}
