"""RNAseq Bench Web — the single-cell pipeline, running inside the browser.

Why this file exists at all: scanpy, anndata, PyDESeq2 and harmonypy cannot run in WebAssembly. Their
dependency chains (numcodecs, google-crc32c, numba, torch) need compiled extensions that have no
wasm build. What Pyodide does give us is numpy, scipy, scikit-learn, igraph and h5py — so the steps
below are not reimplementations of those packages, they are the same algorithms called directly from
the same libraries the desktop app uses underneath scanpy:

    PCA             sklearn.decomposition.PCA          (scanpy calls this)
    neighbours      sklearn.neighbors.NearestNeighbors (exact, where scanpy approximates)
    clustering      igraph.community_leiden            (the Leiden algorithm itself)
    markers         rank-sum, vectorised                (what sc.tl.rank_genes_groups computes)

UMAP is the exception and is done in JavaScript, from the neighbour graph computed here.

Anything needing a model, a reference or a network — DESeq2, Harmony, scVI, CellTypist, GSEA — is not
here and is not approximated. A worse substitute presented as the real thing is the one outcome worth
avoiding.
"""
from __future__ import annotations

import io
import json

import numpy as np
import scipy.sparse as sp

STATE: dict = {}


def _log(msg: str) -> None:
    print(msg, flush=True)


# ---------------------------------------------------------------- reading data
def load_mtx(mtx_bytes: bytes, features_text: str, barcodes_text: str) -> dict:
    """10x Matrix Market triplet. Stored genes × cells, transposed here to cells × genes."""
    from scipy.io import mmread

    m = mmread(io.BytesIO(mtx_bytes))
    genes = [ln.split("\t")[1] if "\t" in ln else ln
             for ln in features_text.strip().split("\n") if ln.strip()]
    cells = [ln.strip() for ln in barcodes_text.strip().split("\n") if ln.strip()]
    X = sp.csr_matrix(m).T.tocsr()                      # -> cells × genes
    if X.shape[1] != len(genes) and X.shape[0] == len(genes):
        X = X.T.tocsr()
    return _finish_load(X, np.array(genes[:X.shape[1]], dtype=object),
                        np.array(cells[:X.shape[0]], dtype=object))


def load_dense(text: str, sep: str = "") -> dict:
    """A plain counts table: genes as rows, samples/cells as columns (the usual GEO layout)."""
    import pandas as pd

    if not sep:
        head = text[:4000].split("\n")[0]
        sep = "\t" if head.count("\t") >= head.count(",") else ","
    df = pd.read_csv(io.StringIO(text), sep=sep, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").dropna(how="all", axis=1).fillna(0)
    X = sp.csr_matrix(df.values.T.astype(np.float32))   # -> cells × genes
    return _finish_load(X, df.index.to_numpy(dtype=object), df.columns.to_numpy(dtype=object))


def load_h5ad(buf: bytes) -> dict:
    """Enough of the .h5ad layout to read X, var_names and obs_names without anndata."""
    import h5py

    f = h5py.File(io.BytesIO(buf), "r")
    x = f["X"]
    if isinstance(x, h5py.Group):                        # sparse: CSR or CSC on disk
        enc = x.attrs.get("encoding-type", "csr_matrix")
        enc = enc.decode() if isinstance(enc, bytes) else enc
        shape = tuple(x.attrs["shape"])
        M = (sp.csr_matrix if "csr" in enc else sp.csc_matrix)(
            (x["data"][:], x["indices"][:], x["indptr"][:]), shape=shape)
        X = sp.csr_matrix(M)
    else:
        X = sp.csr_matrix(np.asarray(x))

    def names(group, fallback_n):
        g = f.get(group)
        if g is None:
            return np.array([f"{group[:3]}{i}" for i in range(fallback_n)], dtype=object)
        key = g.attrs.get("_index", b"_index")
        key = key.decode() if isinstance(key, bytes) else key
        col = g[key][:] if key in g else list(g.values())[0][:]
        return np.array([c.decode() if isinstance(c, bytes) else str(c) for c in col], dtype=object)

    genes = names("var", X.shape[1])
    cells = names("obs", X.shape[0])
    f.close()
    return _finish_load(X, genes, cells)


def _finish_load(X, genes, cells) -> dict:
    X = X.astype(np.float32)
    X.eliminate_zeros()
    # duplicate gene symbols are common in GEO tables and break every downstream lookup
    genes = np.asarray(genes, dtype=object)
    _, first = np.unique(genes, return_index=True)
    if len(first) < len(genes):
        keep = np.zeros(len(genes), bool)
        keep[np.sort(first)] = True
        X, genes = X[:, keep], genes[keep]
    STATE.clear()
    STATE.update(X=X, genes=genes, cells=np.asarray(cells, dtype=object), raw_shape=X.shape)
    return qc_metrics()


# ---------------------------------------------------------------- QC
def qc_metrics() -> dict:
    X = STATE["X"]
    genes = STATE["genes"]
    up = np.array([str(g).upper() for g in genes])
    n_genes = np.diff(X.indptr).astype(np.int32)
    total = np.asarray(X.sum(1)).ravel()
    mt = np.char.startswith(up.astype(str), "MT-")
    pct_mt = (np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(total, 1) * 100) if mt.any() else np.zeros(X.shape[0], np.float32)
    STATE.update(n_genes=n_genes, total=total, pct_mt=pct_mt, has_mt=bool(mt.any()))

    lg = np.log1p(n_genes.astype(np.float64))
    mad = lambda v: float(np.median(np.abs(v - np.median(v))) * 1.4826)   # noqa: E731
    auto_min = max(100, int(np.expm1(np.median(lg) - 5 * mad(lg))))
    auto_max = min(int(np.expm1(np.median(lg) + 5 * mad(lg))), int(n_genes.max()))
    auto_mt = float(np.clip(np.median(pct_mt) + 3 * mad(pct_mt), 5, 25)) if mt.any() else 100.0
    return {
        "n_cells": int(X.shape[0]), "n_genes_total": int(X.shape[1]),
        "median_genes": int(np.median(n_genes)), "median_umi": int(np.median(total)),
        "has_mt": bool(mt.any()),
        "auto": {"min_genes": auto_min, "max_genes": auto_max, "max_mt": round(auto_mt, 1)},
        "hist": {
            "genes": np.histogram(n_genes, bins=40)[0].tolist(),
            "genes_edges": np.histogram(n_genes, bins=40)[1].round(1).tolist(),
            "mt": np.histogram(pct_mt, bins=40)[0].tolist(),
            "mt_edges": np.histogram(pct_mt, bins=40)[1].round(2).tolist(),
        },
        "sample": {
            "n_genes": n_genes[:6000].tolist(),
            "total": total[:6000].round(1).tolist(),
            "pct_mt": pct_mt[:6000].round(2).tolist(),
        },
    }


def apply_qc(min_genes: int, max_genes: int, max_mt: float, min_cells: int = 3) -> dict:
    keep = (STATE["n_genes"] >= min_genes) & (STATE["n_genes"] <= max_genes) & (STATE["pct_mt"] <= max_mt)
    if keep.sum() < 50:
        raise ValueError(f"Only {int(keep.sum())} cells pass those cut-offs — loosen them.")
    X = STATE["X"][keep]
    gene_ok = np.asarray((X > 0).sum(0)).ravel() >= min_cells
    X = X[:, gene_ok]
    STATE.update(Xf=X.tocsr(), genes_f=STATE["genes"][gene_ok], cells_f=STATE["cells"][keep], keep=keep)
    return {"cells_kept": int(keep.sum()), "cells_removed": int((~keep).sum()),
            "genes_kept": int(gene_ok.sum()), "genes_removed": int((~gene_ok).sum())}


# ---------------------------------------------------------------- normalize, HVGs, PCA
def normalize_hvg_pca(n_hvg: int = 2000, n_pcs: int = 50) -> dict:
    from sklearn.decomposition import PCA

    X = STATE["Xf"].astype(np.float32)
    total = np.asarray(X.sum(1)).ravel()
    inv = (1e4 / np.maximum(total, 1)).astype(np.float32)
    Xn = sp.diags(inv) @ X                                # counts per 10,000
    Xn.data = np.log1p(Xn.data)                           # then log1p, as LogNormalize does
    STATE["Xn"] = Xn.tocsr()

    # Seurat-style dispersion HVGs: variance high for the expression level, binned so that highly
    # expressed genes are not automatically selected
    n = Xn.shape[0]
    mean = np.asarray(Xn.mean(0)).ravel()
    sq = np.asarray(Xn.multiply(Xn).mean(0)).ravel()
    var = np.maximum(sq - mean ** 2, 0) * n / max(n - 1, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        disp = np.where(mean > 0, var / mean, 0)
    ok = (mean > 0) & (disp > 0)
    bins = np.zeros(len(mean), int)
    if ok.sum() > 40:
        q = np.quantile(mean[ok], np.linspace(0, 1, 21))
        bins = np.clip(np.digitize(mean, np.unique(q)) - 1, 0, 19)
    z = np.zeros(len(mean))
    for b in np.unique(bins[ok]):
        m = ok & (bins == b)
        d = np.log1p(disp[m])
        z[m] = (d - d.mean()) / (d.std() + 1e-9)
    n_hvg = int(min(n_hvg, max(50, ok.sum() - 1)))
    hvg = np.zeros(len(mean), bool)
    hvg[np.argsort(-np.where(ok, z, -np.inf))[:n_hvg]] = True
    STATE["hvg"] = hvg

    # scale to unit variance (clipped at 10, as scanpy does) before PCA
    sub = Xn[:, hvg].toarray()
    sub -= sub.mean(0)
    sd = sub.std(0)
    sub /= np.where(sd > 0, sd, 1)
    np.clip(sub, -10, 10, out=sub)
    n_pcs = int(min(n_pcs, min(sub.shape) - 1))
    p = PCA(n_components=n_pcs, svd_solver="randomized", random_state=0)
    pcs = p.fit_transform(sub).astype(np.float32)
    STATE.update(pca=pcs, var_ratio=p.explained_variance_ratio_)

    vr = p.explained_variance_ratio_
    d = np.abs(np.arange(len(vr)) * (vr[-1] - vr[0]) - (len(vr) - 1) * (vr - vr[0]))
    knee = int(np.argmax(d)) + 1
    return {"n_hvg": int(hvg.sum()), "n_pcs": int(n_pcs),
            "knee": int(np.clip(knee + 5, 10, 40)),
            "variance_ratio": (vr * 100).round(3).tolist(),
            "top_hvg": [str(g) for g in STATE["genes_f"][np.argsort(-np.where(ok, z, -np.inf))[:20]]]}


# ---------------------------------------------------------------- neighbours + Leiden
def neighbors_and_cluster(n_pcs: int = 0, k: int = 15, resolution: float = 0.5) -> dict:
    from sklearn.neighbors import NearestNeighbors
    import igraph as ig

    P = STATE["pca"][:, :int(n_pcs)] if n_pcs else STATE["pca"]
    n = P.shape[0]
    k = int(min(k, n - 1))
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="brute", metric="euclidean").fit(P)
    dist, idx = nn.kneighbors(P)                         # exact; scanpy approximates with pynndescent
    STATE.update(knn_idx=idx[:, 1:], knn_dist=dist[:, 1:])

    rows = np.repeat(np.arange(n), k)
    cols = idx[:, 1:].ravel()
    # edge weight from distance, symmetrised: closer neighbours count for more
    d = dist[:, 1:].ravel()
    scale = np.maximum(dist[:, 1:].max(1), 1e-9).repeat(k)
    w = np.clip(1.0 - d / scale, 1e-3, 1.0)
    g = ig.Graph(n=n, edges=list(zip(rows.tolist(), cols.tolist())), directed=False)
    g.es["weight"] = w.tolist()
    g.simplify(combine_edges="max")
    part = g.community_leiden(objective_function="modularity", weights="weight",
                              resolution=float(resolution), n_iterations=3)
    lab = np.asarray(part.membership, dtype=np.int32)
    order = np.argsort(-np.bincount(lab))                # biggest cluster becomes 0, as scanpy does
    remap = np.zeros(order.max() + 1, np.int32)
    remap[order] = np.arange(len(order))
    lab = remap[lab]
    STATE["leiden"] = lab
    return {"n_clusters": int(lab.max() + 1),
            "sizes": np.bincount(lab).tolist(),
            "knn_idx": STATE["knn_idx"].astype(np.int32).ravel().tolist(),
            "knn_dist": STATE["knn_dist"].astype(np.float32).ravel().tolist(),
            "k": k}


def set_umap(coords_flat) -> None:
    xy = np.asarray(coords_flat, dtype=np.float32).reshape(-1, 2)
    STATE["umap"] = xy


# ---------------------------------------------------------------- marker genes
def markers(top: int = 8) -> dict:
    """Wilcoxon rank-sum of each cluster against the rest — what rank_genes_groups computes."""
    Xn, lab = STATE["Xn"], STATE["leiden"]
    n, n_g = Xn.shape
    genes = STATE["genes_f"]
    dense_ok = n_g <= 6000
    ranks = None
    if dense_ok:                                          # rank once, reuse for every cluster
        A = Xn.toarray()
        ranks = np.apply_along_axis(_rankdata, 0, A)
    out, table = {}, []
    for c in range(int(lab.max()) + 1):
        m = lab == c
        n1, n2 = int(m.sum()), int((~m).sum())
        if n1 < 3 or n2 < 3:
            continue
        if ranks is not None:
            R1 = ranks[m].sum(0)
        else:
            R1 = np.zeros(n_g)
            for j0 in range(0, n_g, 2000):               # chunked so memory stays bounded
                blk = Xn[:, j0:j0 + 2000].toarray()
                R1[j0:j0 + 2000] = np.apply_along_axis(_rankdata, 0, blk)[m].sum(0)
        U = R1 - n1 * (n1 + 1) / 2.0
        mu = n1 * n2 / 2.0
        sd = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
        z = (U - mu) / max(sd, 1e-9)
        in_mean = np.asarray(Xn[m].mean(0)).ravel()
        out_mean = np.asarray(Xn[~m].mean(0)).ravel()
        lfc = np.log2((np.expm1(in_mean) + 1e-9) / (np.expm1(out_mean) + 1e-9))
        pct_in = np.asarray((Xn[m] > 0).mean(0)).ravel() * 100
        pct_out = np.asarray((Xn[~m] > 0).mean(0)).ravel() * 100
        best = np.argsort(-z)[:top]
        out[str(c)] = [str(genes[i]) for i in best]
        for i in best[:5]:
            table.append([str(c), str(genes[i]), round(float(lfc[i]), 2),
                          round(float(pct_in[i]), 1), round(float(pct_out[i]), 1), round(float(z[i]), 1)])
    STATE["top_markers"] = out
    return {"per_cluster": out, "table": table}


def _rankdata(v):
    o = np.argsort(v, kind="stable")
    r = np.empty(len(v), np.float64)
    r[o] = np.arange(1, len(v) + 1)
    s = v[o]
    i = 0                                                 # average ranks within ties, as the test requires
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            r[o[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    return r


def gene_values(name: str) -> dict:
    g = np.where(np.array([str(x).upper() for x in STATE["genes_f"]]) == str(name).upper())[0]
    if not len(g):
        return {"found": False}
    v = np.asarray(STATE["Xn"][:, g[0]].todense()).ravel()
    return {"found": True, "gene": str(STATE["genes_f"][g[0]]),
            "values": np.round(v, 3).tolist(), "pct": round(float((v > 0).mean() * 100), 1)}


def embedding_payload() -> str:
    xy = STATE.get("umap")
    lab = STATE["leiden"]
    return json.dumps({
        "x": np.round(xy[:, 0], 3).tolist(), "y": np.round(xy[:, 1], 3).tolist(),
        "cluster": lab.tolist(),
        "n_genes": STATE["n_genes"][STATE["keep"]].tolist(),
        "pct_mt": np.round(STATE["pct_mt"][STATE["keep"]], 2).tolist(),
        "genes": [str(g) for g in STATE["genes_f"]],
        "top_markers": STATE.get("top_markers", {}),
    })
