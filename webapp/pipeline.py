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
    STATE["n_pcs_used"] = int(P.shape[1])
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


# ---------------------------------------------------------------- figures
# matplotlib is in Pyodide, so the web build can draw the same figures as the desktop app rather than
# hand-rolled SVG. Each returns a base64 PNG the page drops straight into an <img>.
def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 130, "savefig.dpi": 130, "font.size": 8,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#9aa8b0", "axes.labelcolor": "#33454f",
                         "xtick.color": "#6a7b84", "ytick.color": "#6a7b84",
                         "axes.titlesize": 9, "figure.facecolor": "white"})
    return plt


def _png(fig) -> str:
    import base64
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


ACCENT, BAD = "#0f766e", "#c9304f"


def _violin(ax, values, keep, title, lines=(), log=False):
    v = np.asarray(values, float)
    t = np.log10(v + 1) if log else v
    parts = ax.violinplot(t, positions=[0], widths=0.85, showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(ACCENT); b.set_edgecolor(ACCENT); b.set_alpha(0.25)
    rng = np.random.default_rng(0)
    pick = rng.choice(len(t), size=min(len(t), 3000), replace=False)
    jit = rng.uniform(-0.3, 0.3, len(pick))
    k = np.asarray(keep)[pick]
    ax.scatter(jit[k], t[pick][k], s=1.6, c=ACCENT, alpha=0.35, linewidths=0, rasterized=True)
    ax.scatter(jit[~k], t[pick][~k], s=3.2, c=BAD, alpha=0.85, linewidths=0, rasterized=True)
    q1, med, q3 = np.percentile(t, [25, 50, 75])
    ax.plot([0, 0], [q1, q3], color="#222", lw=2.6, solid_capstyle="butt")
    ax.scatter([0], [med], color="white", s=12, zorder=3, edgecolor="#222", linewidth=0.7)
    for ln in lines:
        ax.axhline(np.log10(ln + 1) if log else ln, color=BAD, ls="--", lw=0.9)
    ax.set_xticks([]); ax.set_title(title)
    ax.set_ylabel("log10(1+x)" if log else "")


def qc_figures(min_genes: int, max_genes: int, max_mt: float) -> dict:
    """The QC panel: what each cell looks like, and which ones the cut-offs remove."""
    plt = _style()
    ng, tot, mt = STATE["n_genes"], STATE["total"], STATE["pct_mt"]
    keep = (ng >= min_genes) & (ng <= max_genes) & (mt <= max_mt)
    has_mt = STATE["has_mt"]

    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.7))
    _violin(axs[0], ng, keep, "genes per cell", [min_genes, max_genes], log=True)
    _violin(axs[1], tot, keep, "counts per cell", [], log=True)
    _violin(axs[2], mt, keep, "mitochondrial %", [max_mt] if has_mt else [])
    fig.tight_layout()
    violins = _png(fig)

    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    rng = np.random.default_rng(0)
    pick = rng.choice(len(ng), size=min(len(ng), 8000), replace=False)
    s = ax.scatter(tot[pick], ng[pick], c=mt[pick], s=3, cmap="viridis", linewidths=0,
                   vmax=max(10, float(np.percentile(mt, 99))), rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("counts per cell"); ax.set_ylabel("genes per cell")
    ax.axhline(min_genes, color=BAD, ls="--", lw=0.9); ax.axhline(max_genes, color=BAD, ls="--", lw=0.9)
    ax.set_title("counts vs genes")
    fig.colorbar(s, ax=ax, fraction=0.04, pad=0.02, label="% mito")
    fig.tight_layout()
    scatter = _png(fig)

    # Highest-expressed genes: the ambient-RNA view. A gene taking a large share of all counts and
    # detected in nearly every cell is soup from lysed cells, not biology.
    X = STATE["X"]
    share = np.asarray(X.sum(0)).ravel()
    share = share / max(share.sum(), 1)
    det = np.asarray((X > 0).mean(0)).ravel()
    top = np.argsort(-share)[:20]
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    ax.barh(range(len(top)), share[top][::-1] * 100, color=ACCENT, alpha=0.85)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([str(STATE["genes"][i]) for i in top][::-1], fontsize=7)
    ax.set_xlabel("% of all counts"); ax.set_title("highest-expressed genes")
    fig.tight_layout()
    highest = _png(fig)

    soup = [{"gene": str(STATE["genes"][i]), "share": round(float(share[i] * 100), 2),
             "detected": round(float(det[i] * 100), 1)}
            # a true soup gene is in nearly every cell, not half of them: at 50% this flagged a
            # cluster marker, while real ambient (TTR in choroid plexus) sits at 9.5% of counts in 100%
            for i in top[:5] if share[i] > 0.01 and det[i] > 0.75]

    fails = {
        "too_few_genes": int((ng < min_genes).sum()),
        "too_many_genes": int((ng > max_genes).sum()),
        "high_mito": int((mt > max_mt).sum()) if has_mt else 0,
    }
    return {"violins": violins, "scatter": scatter, "highest": highest,
            "kept": int(keep.sum()), "removed": int((~keep).sum()), "fails": fails, "soup": soup}


# ---------------------------------------------------------------- per-cell scores
S_GENES = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7",
           "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
           "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6",
           "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"]
G2M_GENES = ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2", "NUF2",
             "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
             "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "CDC20",
             "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23",
             "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3"]
# immediate-early and heat-shock genes induced by warm dissociation (van den Brink et al. 2017)
DISSOCIATION = ["FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "IER2", "IER3", "DUSP1", "ZFP36",
                "SOCS3", "KLF6", "BTG2", "NR4A1", "PPP1R15A", "CEBPB", "RHOB", "SGK1", "NFKBIA",
                "HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "HSPH1", "HSP90AA1", "DNAJB1", "SQSTM1", "UBC"]


def _score(names: list[str]) -> tuple[np.ndarray, int]:
    """Mean scaled expression of a gene set minus a matched-expression control, as score_genes does.

    Also returns how many of the set were actually found. A score computed from two genes, or none, is
    not a weak result — it is no result, and the caller has to be able to tell the difference.
    """
    Xn = STATE["Xn"]
    up = np.array([str(g).upper() for g in STATE["genes_f"]])
    idx = np.where(np.isin(up, [g.upper() for g in names]))[0]
    if len(idx) < 5:
        return np.zeros(Xn.shape[0], np.float32), int(len(idx))
    mean = np.asarray(Xn.mean(0)).ravel()
    order = np.argsort(mean)
    rank = np.empty_like(order); rank[order] = np.arange(len(order))
    ctrl = []
    for i in idx:                                        # 50 control genes at a similar level
        lo = max(0, rank[i] - 60); hi = min(len(order), rank[i] + 60)
        ctrl.extend(order[lo:hi])
    ctrl = np.setdiff1d(np.unique(ctrl), idx)
    a = np.asarray(Xn[:, idx].mean(1)).ravel()
    b = np.asarray(Xn[:, ctrl].mean(1)).ravel() if len(ctrl) else 0.0
    return (a - b).astype(np.float32), int(len(idx))


def cell_scores() -> dict:
    """Cell cycle and dissociation stress, per cell and summarised per cluster."""
    lab = STATE["leiden"]
    s, n_s = _score(S_GENES)
    g2m, n_g = _score(G2M_GENES)
    diss, n_d = _score(DISSOCIATION)
    cc_ok = n_s >= 5 and n_g >= 5
    diss_ok = n_d >= 5
    # Without the genes there is no phase to assign. Calling every cell G2/M because both scores are
    # zero would be a confident answer built on nothing.
    phase = (np.where((s < 0) & (g2m < 0), "G1", np.where(s > g2m, "S", "G2M"))
             if cc_ok else np.full(len(s), "unknown", dtype=object))
    STATE.update(s_score=s, g2m_score=g2m, diss=diss, phase=phase)

    per = {}
    for c in range(int(lab.max()) + 1):
        m = lab == c
        per[str(c)] = {"cycling": round(float((phase[m] != "G1").mean() * 100), 1) if cc_ok else None,
                       "stress": round(float(diss[m].mean()), 3) if diss_ok else None}
    out = {"per_cluster": per, "cc_ok": cc_ok, "diss_ok": diss_ok,
           "genes_found": {"S": n_s, "G2M": n_g, "dissociation": n_d},
           "phase_counts": ({p: int((phase == p).sum()) for p in ("G1", "S", "G2M")} if cc_ok else None),
           "diss": np.round(diss, 3).tolist() if diss_ok else None}
    if diss_ok:
        vals = [v["stress"] for v in per.values()]
        spread = float(np.max(vals) - np.median(vals))
        out.update(stress_spread=round(spread, 3),
                   stress_cluster=max(per, key=lambda c: per[c]["stress"]),
                   stress_flag=bool(spread > 0.25))
    return out


def hvg_elbow_figures() -> dict:
    """Feature selection and how many components the elbow justifies."""
    plt = _style()
    Xn, hvg = STATE["Xn"], STATE["hvg"]
    n = Xn.shape[0]
    mean = np.asarray(Xn.mean(0)).ravel()
    sq = np.asarray(Xn.multiply(Xn).mean(0)).ravel()
    var = np.maximum(sq - mean ** 2, 0) * n / max(n - 1, 1)
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    ax.scatter(mean[~hvg], var[~hvg], s=2, c="#c7ccd1", linewidths=0, rasterized=True, label="other")
    ax.scatter(mean[hvg], var[hvg], s=2.5, c=BAD, linewidths=0, rasterized=True, label="selected")
    ax.set_xscale("symlog", linthresh=1e-3); ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_xlabel("mean expression"); ax.set_ylabel("variance")
    ax.set_title(f"{int(hvg.sum()):,} variable genes kept"); ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    hvg_png = _png(fig)

    vr = STATE["var_ratio"] * 100
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    ax.plot(np.arange(1, len(vr) + 1), vr, "o-", ms=3, color=ACCENT, lw=1)
    ax.axvline(STATE.get("n_pcs_used", 0) + 0.5, color="#888", ls="--", lw=0.9)
    ax.set_xlabel("principal component"); ax.set_ylabel("% variance")
    ax.set_title("scree plot")
    fig.tight_layout()
    return {"hvg": hvg_png, "elbow": _png(fig)}
