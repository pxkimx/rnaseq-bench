"""Single-cell RNA-seq pipeline built on Scanpy (follows the Scanpy/Seurat standard workflow
and the single-cell best-practices book: MAD-based QC, Scrublet, seurat_v3 HVGs, PCA,
kNN graph, UMAP, Leiden, Wilcoxon markers)."""
from __future__ import annotations

import json
import os
import re
from importlib.metadata import version
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from .common import PALETTE, Job, UserFacingError, commas, pct, style, log_exc
from .io_utils import load_single_cell
from .markers import CELL_TYPE_MARKERS, G2M_GENES, S_GENES

ROOT = Path(__file__).resolve().parent.parent

DEFAULTS = dict(min_genes=None, max_genes=None, max_mt=None, min_cells=3, doublets="auto",
                n_hvg=2000, resolution=0.5, n_pcs=None, batch_key=None, integrate=True,
                sample_key=None, condition_key=None, pb_covariates=None, reference=None, pseudobulk=True,
                # expert add-ons (each skipped gracefully when its resource is unavailable)
                celltypist=True, celltypist_model=None, pathways=True, trajectory=True, root_cluster=None, ccc=True,
                gsea=True, genes=None, subset=None, subset_key="leiden", subset_from_job=None)


def mad(x):
    return np.median(np.abs(x - np.median(x))) * 1.4826


def _dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)


def _violin(ax, values, keep, title, lines=(), log=False, color="#2b7a78"):
    v = np.asarray(values, float)
    t = np.log10(v + 1) if log else v
    parts = ax.violinplot(t, positions=[0], widths=0.8, showextrema=False)
    for b in parts["bodies"]:
        b.set_facecolor(color); b.set_edgecolor(color); b.set_alpha(0.25)
    rng = np.random.default_rng(0)
    idx = rng.choice(len(t), size=min(len(t), 4000), replace=False)
    jit = rng.uniform(-0.28, 0.28, len(idx))
    k = np.asarray(keep)[idx]
    ax.scatter(jit[k], t[idx][k], s=1.2, c=color, alpha=0.35, linewidths=0, rasterized=True)
    ax.scatter(jit[~k], t[idx][~k], s=2.5, c="#c9304f", alpha=0.8, linewidths=0, rasterized=True)
    q1, med, q3 = np.percentile(t, [25, 50, 75])
    ax.plot([0, 0], [q1, q3], color="#222", lw=3, solid_capstyle="butt")
    ax.scatter([0], [med], color="white", s=14, zorder=3, edgecolor="#222", linewidth=0.8)
    for ln in lines:
        if ln is None:
            continue
        y = np.log10(ln + 1) if log else ln
        ax.axhline(y, color="#c9304f", ls="--", lw=1)
        ax.text(0.47, y, f" {ln:,.0f}" if ln >= 10 else f" {ln:.1f}", color="#c9304f", va="bottom", ha="right", fontsize=8, transform=ax.get_yaxis_transform())
    ax.set_xticks([]); ax.set_title(title)
    if log:
        ticks = [t_ for t_ in [10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000] if np.log10(t_ + 1) <= t.max() * 1.02 and np.log10(t_ + 1) >= t.min() * 0.98]
        ax.set_yticks([np.log10(x + 1) for x in ticks]); ax.set_yticklabels([f"{x:,}" for x in ticks])
    ax.spines["bottom"].set_visible(False)


def annotate_clusters(adata, cluster_key="leiden", layer="log1p"):
    X = adata.layers[layer]
    upper = {g.upper(): i for i, g in enumerate(adata.var_names)}
    clusters = adata.obs[cluster_key].cat.categories
    means = np.vstack([np.asarray(X[(adata.obs[cluster_key] == c).values].mean(axis=0)).ravel() for c in clusters])
    frac = np.vstack([np.asarray((X[(adata.obs[cluster_key] == c).values] > 0).mean(axis=0)).ravel() for c in clusters])
    mu, sd = means.mean(0), means.std(0) + 1e-9
    z = (means - mu) / sd
    scores = {}
    for ct, genes in CELL_TYPE_MARKERS.items():
        ids = [upper[g] for g in genes if g in upper]
        if len(ids) < 2:
            continue
        # z-score across clusters, weighted by detection so absent genes don't dominate
        scores[ct] = (z[:, ids] * np.clip(frac[:, ids] * 2, 0, 1)).mean(1)
    if not scores:
        return {c: ("Unassigned", 0.0, []) for c in clusters}, None
    S = pd.DataFrame(scores, index=clusters)
    out, used = {}, {}
    for c in clusters:
        best = S.loc[c].idxmax(); val = float(S.loc[c].max())
        second = float(S.loc[c].drop(best).max()) if S.shape[1] > 1 else 0
        if val < 0.8:
            out[c] = ("Unassigned", val, [])
            continue
        ev = [g for g in CELL_TYPE_MARKERS[best] if g in upper and frac[list(clusters).index(c), upper[g]] > 0.3]
        label = best if val - second > 0.15 else f"{best} / {S.loc[c].drop(best).idxmax()}?"
        used[label] = used.get(label, 0) + 1
        out[c] = (label, val, ev)
    # number duplicates
    seen = {}
    for c in clusters:
        l = out[c][0]
        if l != "Unassigned" and used.get(l, 0) > 1:
            seen[l] = seen.get(l, 0) + 1
            out[c] = (f"{l} ({seen[l]})", out[c][1], out[c][2])
    return out, S


def _fig_of(ret):
    """Scanpy plotting functions return an Axes, a dict of Axes, or a plot object."""
    if ret is None:
        return plt.gcf()
    if isinstance(ret, dict):
        ret = next(v for v in ret.values() if v is not None)
    if hasattr(ret, "figure"):
        return ret.figure
    if hasattr(ret, "get_axes"):
        return ret.get_axes()["mainplot_ax"].figure
    return plt.gcf()


def _user_feature_plots(job, adata, genes):
    """UMAP feature plots for genes the user asked about."""
    from .sc_extra import _lognorm
    upper = {g.upper(): g for g in adata.var_names}
    want = [upper[g.strip().upper()] for g in (genes if isinstance(genes, list) else str(genes).split(",")) if g.strip().upper() in upper][:12]
    missing = [g for g in (genes if isinstance(genes, list) else str(genes).split(",")) if g.strip() and g.strip().upper() not in upper]
    if not want:
        job.flag("warn", f"None of the requested genes are in this dataset: {', '.join(missing[:6])}.")
        return
    b = _lognorm(adata); xy = b.obsm["X_umap"]
    cols = min(4, len(want)); rows = int(np.ceil(len(want) / cols))
    fig, axs = plt.subplots(rows, cols, figsize=(3.2 * cols, 2.9 * rows), squeeze=False)
    for ax, g in zip(axs.ravel(), want):
        v = b[:, g].X; v = np.asarray(v.toarray() if hasattr(v, "toarray") else v).ravel()
        o = np.argsort(v, kind="stable")
        sca = ax.scatter(xy[o, 0], xy[o, 1], c=v[o], s=2, cmap="viridis", linewidths=0, rasterized=True, vmax=np.percentile(v, 99) or 1)
        pct = 100 * (v > 0).mean()
        ax.set_title(f"{g} · {pct:.0f}% of cells", fontsize=8); ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(sca, ax=ax, fraction=0.04, pad=0.02).ax.tick_params(labelsize=6)
    for ax in axs.ravel()[len(want):]:
        ax.axis("off")
    fig.tight_layout()
    # per-label expression table
    lab = adata.obs["cell_type_simple"].astype(str) if "cell_type_simple" in adata.obs else adata.obs["leiden"].astype(str)
    rows_ = []
    for g in want:
        v = b[:, g].X; v = np.asarray(v.toarray() if hasattr(v, "toarray") else v).ravel()
        d = pd.DataFrame({"v": v, "l": lab.values}).groupby("l")["v"].agg(["mean", lambda x: 100 * (x > 0).mean()])
        best = d["mean"].idxmax()
        rows_.append([g, best, f"{d.loc[best, 'mean']:.2f}", f"{d.iloc[:, 1].loc[best]:.0f}%"])
    job.figure("genes_of_interest", "Your genes on the UMAP", fig, wide=len(want) > 4,
               how="Log-normalized expression per cell; brightest = 99th percentile.",
               yours=(f"Not found: {', '.join(missing[:6])}. " if missing else "") + "; ".join(f"<b>{r[0]}</b> highest in {r[1]}" for r in rows_) + ".")
    job.table("goi_table", "Where your genes are expressed", ["Gene", "Highest population", "Mean log expr.", "% cells expressing"], rows_)


def knee(y):
    y = np.asarray(y); x = np.arange(len(y))
    p1, p2 = np.array([0, y[0]]), np.array([len(y) - 1, y[-1]])
    v = p2 - p1
    pts = np.vstack([x, y]).T - p1
    d = np.abs(v[0] * pts[:, 1] - v[1] * pts[:, 0]) / np.linalg.norm(v)  # 2-D cross product (np.cross dropped 2-D support)
    return int(np.argmax(d)) + 1


def run(job: Job, files: list[Path], params: dict):
    P = {**DEFAULTS, **{k: v for k, v in (params or {}).items() if v not in (None, "")}}
    style()
    sc.settings.verbosity = 0
    sc.settings.figdir = str(job.fig_dir)
    job.step(3, "Reading data")
    adata, notes = load_single_cell(files, job.root)
    if adata.uns.get("empty_droplets_removed"):
        notes.append(f"{int(adata.uns['empty_droplets_removed']):,} barcodes with < 200 UMIs (empty droplets from a raw matrix) were dropped before QC.")
    # sub-clustering: keep only cells that a previous job put in the chosen clusters / labels
    if P["subset"] and P["subset_from_job"]:
        try:
            prev = Path(os.environ.get("TL_HOME", ROOT)) / "jobs" / str(P["subset_from_job"]) / "cell_metadata.csv"
            cm = pd.read_csv(prev, index_col=0, dtype=str)
            want = [x.strip() for x in (P["subset"] if isinstance(P["subset"], list) else str(P["subset"]).split(",")) if x.strip()]
            key = P["subset_key"] if P["subset_key"] in cm.columns else "leiden"
            keep_ids = set(cm.index[cm[key].astype(str).isin(want)])
            m = adata.obs_names.isin(keep_ids)
            if m.sum() < 50:
                raise UserFacingError(f"Only {int(m.sum())} cells match {key} ∈ {want} in job {P['subset_from_job']}.")
            adata = adata[m].copy()
            notes.append(f"Sub-clustering: kept {adata.n_obs:,} cells with {key} in {{{', '.join(want)}}} from job {P['subset_from_job']}; QC, HVGs, PCA and clustering are recomputed on this subset.")
        except UserFacingError:
            raise
        except Exception:  # noqa: BLE001
            raise UserFacingError(f"Could not apply the sub-clustering selection ({log_exc('subset')}).")
    N0, G0 = adata.n_obs, adata.n_vars

    # ---------------------------------------------------------------- QC metrics
    job.step(8, "Computing QC metrics")
    vn = adata.var_names.str.upper()
    adata.var["mt"] = vn.str.startswith("MT-")
    if "chrom" in adata.var:
        adata.var["mt"] = adata.var["mt"] | adata.var["chrom"].astype(str).isin(["MT", "chrM", "M"])
    adata.var["ribo"] = vn.str.match(r"^RP[SL]\d")
    adata.var["hb"] = vn.str.match(r"^HB[ABDEGMQZ]\d?$")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20], log1p=True, inplace=True)
    has_mt = bool(adata.var["mt"].any())
    species = "mouse" if (adata.var_names.str.match(r"^[A-Z][a-z]").mean() > 0.5) else "human"

    batch_key = P["batch_key"]
    if not batch_key:
        for k in ("sample", "plate", "plate_id", "batch", "orig.ident", "donor", "library", "sample_id", "patient", "gsm"):
            if k in adata.obs and 1 < adata.obs[k].nunique() <= 50:
                batch_key = k
                break
    if batch_key and batch_key not in adata.obs:
        batch_key = None
    nb = adata.obs[batch_key].nunique() if batch_key else 1

    lg = adata.obs["log1p_n_genes_by_counts"].values
    lt = adata.obs["log1p_total_counts"].values
    mtv = adata.obs["pct_counts_mt"].values
    auto_min = max(100, int(np.expm1(np.median(lg) - 5 * mad(lg))))
    auto_max = int(np.expm1(np.median(lg) + 5 * mad(lg)))
    auto_mt = float(np.clip(np.median(mtv) + 3 * mad(mtv), 5, 25)) if has_mt else 100.0
    min_g = int(P["min_genes"] or auto_min)
    max_g = int(P["max_genes"] or auto_max)
    max_mt = float(P["max_mt"] or auto_mt)
    top20 = adata.obs["pct_counts_in_top_20_genes"].values
    top_cut = np.median(top20) + 5 * mad(top20)

    low = adata.obs["n_genes_by_counts"] < min_g
    high = adata.obs["n_genes_by_counts"] > max_g
    mt_bad = adata.obs["pct_counts_mt"] > max_mt
    top_bad = pd.Series(top20 > top_cut, index=adata.obs_names)
    keep = ~(low | high | mt_bad | top_bad)
    adata.obs["qc_pass"] = keep.values
    job.result["params"] = {"min_genes": min_g, "max_genes": max_g, "max_mt": round(max_mt, 1),
                            "resolution": P["resolution"], "n_hvg": P["n_hvg"], "doublets": bool(P["doublets"]),
                            "batch_key": batch_key or "", "sample_key": P["sample_key"] or "", "condition_key": P["condition_key"] or "",
                            "obs_columns": [c for c in adata.obs.columns if adata.obs[c].dtype.kind in "OSUb" or adata.obs[c].nunique() < 30][:40], "auto": {"min_genes": auto_min, "max_genes": auto_max, "max_mt": round(auto_mt, 1)}}

    job.section("qc", "Step 1", "Quality control",
                "Each barcode is scored on library size, genes detected and mitochondrial fraction. Cut-offs default to "
                "median ± 5 MAD (genes) and median + 3 MAD (mito, 5–25%), the adaptive thresholds recommended in "
                "<i>Single-cell best practices</i> (Heumos et al., 2023).")
    fig, axs = plt.subplots(1, 3, figsize=(10, 3.6))
    kp = keep.values
    _violin(axs[0], adata.obs["n_genes_by_counts"], kp, "n_genes_by_counts", [min_g, max_g], log=True)
    _violin(axs[1], adata.obs["total_counts"], kp, "total_counts", [], log=True)
    _violin(axs[2], adata.obs["pct_counts_mt"], kp, "pct_counts_mt", [max_mt] if has_mt else [])
    fig.tight_layout()
    job.figure("qc_violin", "QC metrics per cell", fig, wide=True,
               sub="red points = removed · dashed = threshold · bar = IQR, dot = median",
               how="Low gene counts usually mean empty droplets or debris; an upper tail can be doublets. "
                   "A high mitochondrial fraction marks cells whose membrane ruptured, letting cytoplasmic mRNA escape.",
               yours=f"Median <b>{commas(np.median(adata.obs['n_genes_by_counts']))}</b> genes, "
                     f"<b>{commas(np.median(adata.obs['total_counts']))}</b> UMIs and "
                     f"<b>{np.median(mtv):.1f}%</b> mitochondrial reads per cell.")

    if batch_key and nb > 1:
        fig, axs = plt.subplots(1, 3, figsize=(min(16, max(9, 3 + nb * 1.2)), 3.4))
        for ax, key in zip(axs, ["n_genes_by_counts", "total_counts", "pct_counts_mt"]):
            groups = adata.obs[batch_key].astype(str)
            cats = sorted(groups.unique())
            data = [adata.obs.loc[groups == c, key].values for c in cats]
            parts = ax.violinplot(data, showextrema=False, showmedians=True)
            for i, b in enumerate(parts["bodies"]):
                b.set_facecolor(PALETTE[i % 20]); b.set_alpha(0.55)
            ax.set_xticks(range(1, len(cats) + 1)); ax.set_xticklabels(cats, rotation=60, ha="right", fontsize=8)
            ax.set_title(key)
            if key != "pct_counts_mt":
                ax.set_yscale("log")
        fig.tight_layout()
        job.figure("qc_by_sample", f"QC by {batch_key}", fig, wide=True,
                   how="Samples processed together should have similar distributions. A sample with a shifted distribution may need its own thresholds.",
                   yours=f"{nb} values of <b>{batch_key}</b> found.")

    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    o = np.argsort(mtv)
    s = ax.scatter(adata.obs["total_counts"].values[o], adata.obs["n_genes_by_counts"].values[o], c=mtv[o], s=3,
                   cmap="viridis", vmax=max(10, np.percentile(mtv, 99)), linewidths=0, rasterized=True)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("total_counts"); ax.set_ylabel("n_genes_by_counts")
    ax.axhline(min_g, color="#c9304f", ls="--", lw=1); ax.axhline(max_g, color="#c9304f", ls="--", lw=1)
    fig.colorbar(s, ax=ax, label="pct_counts_mt", shrink=0.8)
    job.figure("qc_scatter", "Counts vs genes", fig,
               how="Healthy cells follow one tight curve. Bright (high-mito) cells below the curve are damaged; a cloud in the bottom-left is ambient RNA or empty droplets.",
               yours=f"{pct((low | mt_bad).sum(), N0)} of barcodes fall in low-quality regions.")

    fig = _fig_of(sc.pl.highest_expr_genes(adata, n_top=20, show=False)); fig.set_size_inches(4.6, 4.4)
    top_share = float(np.sort(np.asarray(adata.X.sum(0)).ravel())[::-1][:20].sum() / adata.X.sum())
    job.figure("qc_top_genes", "Highest expressed genes", fig,
               how="Share of counts per cell taken by each gene. MALAT1, mitochondrial and ribosomal genes commonly top the list; a single gene dominating suggests contamination.",
               yours=f"The top 20 genes take <b>{100*top_share:.1f}%</b> of all counts; ribosomal genes average {np.median(adata.obs['pct_counts_ribo']):.1f}% per cell.")

    removed = int((~keep).sum())
    adata = adata[keep.values].copy()
    sc.pp.filter_genes(adata, min_cells=int(P["min_cells"]))

    # ---------------------------------------------------------------- doublets
    n_doublets = 0
    per_sample = adata.obs[batch_key].value_counts() if batch_key else pd.Series([adata.n_obs])
    plate_based = bool(per_sample.max() <= 400 and np.median(adata.obs["n_genes_by_counts"]) > 3000)
    if P["doublets"] == "auto":
        P["doublets"] = not plate_based
        if plate_based:
            notes.append("Cells per sample ≤ 400 with very deep libraries — this looks plate-based (Smart-seq style), so droplet doublet detection was skipped.")
    if P["doublets"] and P["doublets"] != "auto":
        job.step(18, "Detecting doublets (Scrublet)")
        try:
            sc.pp.scrublet(adata, batch_key=batch_key if nb > 1 else None, random_state=0)
            n_doublets = int(adata.obs["predicted_doublet"].sum())
            fig, ax = plt.subplots(figsize=(4.6, 3.4))
            ds = adata.obs["doublet_score"].values
            ax.hist(ds[~adata.obs["predicted_doublet"].values], bins=60, color="#2b7a78", alpha=0.8, label="singlet")
            ax.hist(ds[adata.obs["predicted_doublet"].values], bins=30, color="#c9304f", alpha=0.9, label="predicted doublet")
            ax.set_yscale("log"); ax.set_xlabel("Scrublet doublet score"); ax.set_ylabel("cells"); ax.legend()
            job.figure("qc_doublets", "Doublet scores", fig,
                       how="Scrublet simulates artificial doublets by summing random pairs of cells, then scores each real cell by how many simulated doublets surround it. A bimodal histogram with a small right-hand mode is typical.",
                       yours=f"<b>{commas(n_doublets)}</b> predicted doublets ({pct(n_doublets, adata.n_obs)}). "
                             f"10x Chromium expects ≈0.8% per 1,000 cells recovered, so ~{0.8*adata.n_obs/1000:.1f}% here.")
            adata = adata[~adata.obs["predicted_doublet"].values].copy()
        except Exception as e:  # noqa: BLE001
            notes.append(f"WARN:Doublet detection skipped ({log_exc('scrublet')}).")
    if adata.n_obs < 50:
        raise UserFacingError(f"Only {adata.n_obs} cells passed QC. Loosen the thresholds (lower min genes or raise max mito %).")

    # ---------------------------------------------------------------- normalize + HVG
    job.step(30, "Normalizing and selecting highly variable genes")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()
    n_hvg = int(min(P["n_hvg"], adata.n_vars - 1))
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat_v3", layer="counts",
                                    batch_key=batch_key if nb > 1 else None, span=0.5)
        hvg_method = "seurat_v3 (variance-stabilized, on raw counts)"
    except Exception:  # noqa: BLE001
        log_exc("hvg seurat_v3")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat", batch_key=batch_key if nb > 1 else None)
        hvg_method = "seurat (dispersion)"

    # cell cycle
    upper = {g.upper(): g for g in adata.var_names}
    s_g = [upper[g] for g in S_GENES if g in upper]
    g2m_g = [upper[g] for g in G2M_GENES if g in upper]
    has_cc = len(s_g) >= 10 and len(g2m_g) >= 10
    if has_cc:
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_g, g2m_genes=g2m_g, random_state=0)

    job.section("hvg", "Step 2", "Feature selection & dimensionality reduction",
                f"Counts were scaled to 10,000 per cell and log-transformed. {n_hvg:,} highly variable genes ({hvg_method}) were scaled and "
                "reduced with PCA; the kNN graph built on the top PCs drives UMAP and Leiden clustering.")
    hv = adata.var
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    xcol = "means"; ycol = "variances_norm" if "variances_norm" in hv else "dispersions_norm"
    m = hv["highly_variable"].values
    ax.scatter(hv[xcol][~m], hv[ycol][~m], s=2, c="#c7ccd1", linewidths=0, rasterized=True, label="other")
    ax.scatter(hv[xcol][m], hv[ycol][m], s=3, c="#2b7a78", linewidths=0, rasterized=True, label="highly variable")
    ax.set_xscale("log"); ax.set_xlabel("mean expression"); ax.set_ylabel("normalized variance")
    top = hv[m].sort_values(ycol, ascending=False).head(12)
    for g, r in top.iterrows():
        ax.annotate(g, (r[xcol], r[ycol]), fontsize=7, xytext=(3, 2), textcoords="offset points")
    ax.legend(markerscale=4, loc="upper left")
    job.figure("hvg", "Highly variable genes", fig,
               how="Genes are ranked by variance after accounting for the mean–variance trend. Variable genes carry cell-type signal; the rest is mostly sampling noise.",
               yours=f"Top-ranked: <b>{', '.join(top.index[:8])}</b>.")

    job.step(40, "Scaling and PCA")
    adata.raw = None
    ad_h = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(ad_h, max_value=10)
    n_comps = int(min(50, ad_h.n_obs - 1, ad_h.n_vars - 1))
    sc.tl.pca(ad_h, n_comps=n_comps, svd_solver="arpack", random_state=0)
    adata.obsm["X_pca"] = ad_h.obsm["X_pca"]
    adata.uns["pca"] = ad_h.uns["pca"]
    vr = ad_h.uns["pca"]["variance_ratio"]
    n_pcs = int(P["n_pcs"] or int(np.clip(knee(vr) + 5, 10, 40)))
    n_pcs = min(n_pcs, n_comps)
    loadings = pd.DataFrame(ad_h.varm["PCs"][:, :3], index=ad_h.var_names, columns=["PC1", "PC2", "PC3"])

    fig, ax = plt.subplots(figsize=(4.6, 3.6))
    ax.plot(np.arange(1, len(vr) + 1), vr * 100, "o-", ms=3, color="#2b7a78", lw=1)
    ax.axvline(n_pcs + 0.5, color="#888", ls="--", lw=1)
    ax.text(n_pcs + 1, vr.max() * 90, f"{n_pcs} PCs used", fontsize=8, color="#555")
    ax.set_yscale("log"); ax.set_xlabel("principal component"); ax.set_ylabel("% variance explained")
    job.figure("pca_elbow", "PCA elbow plot", fig,
               how="Variance captured by each component. Components before the elbow hold biological structure; the flat tail is noise. The cut-off was set 5 PCs past the elbow.",
               yours=f"PC1 explains <b>{100*vr[0]:.1f}%</b>; PCs 1–{n_pcs} together explain {100*vr[:n_pcs].sum():.1f}% of HVG variance.")

    fig, axs = plt.subplots(1, 3, figsize=(10, 3.4))
    for i, ax in enumerate(axs):
        l = loadings.iloc[:, i].sort_values()
        sel = pd.concat([l.head(6), l.tail(6)])
        ax.barh(range(len(sel)), sel.values, color=["#2a64b0"] * 6 + ["#c9304f"] * 6)
        ax.set_yticks(range(len(sel))); ax.set_yticklabels(sel.index, fontsize=7.5)
        ax.set_title(f"PC{i+1} loadings"); ax.axvline(0, color="#999", lw=0.6)
    fig.tight_layout()
    job.figure("pca_loadings", "Genes driving the top PCs", fig, wide=True,
               how="Genes with the largest positive (red) and negative (blue) loadings define each axis — often a PC contrasts two lineages. A PC driven by mitochondrial or ribosomal genes points to technical variation.",
               yours="PC1 separates <b>" + ", ".join(loadings["PC1"].sort_values().index[-3:]) + "</b> from <b>" + ", ".join(loadings["PC1"].sort_values().index[:3]) + "</b>.")

    # integration
    use_rep = "X_pca"
    integrated = False
    if batch_key and nb > 1 and P["integrate"]:
        try:
            import harmonypy  # noqa: F401
            job.step(50, f"Integrating {nb} batches with Harmony")
            import logging
            logging.getLogger("harmonypy").setLevel(logging.WARNING)
            ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, [batch_key], max_iter_harmony=20, random_state=0)
            Z = np.asarray(ho.Z_corr)
            adata.obsm["X_pca_harmony"] = Z if Z.shape[0] == adata.n_obs else Z.T
            use_rep = "X_pca_harmony"; integrated = True
        except Exception as e:  # noqa: BLE001
            notes.append(f"WARN:Harmony integration unavailable ({type(e).__name__}); clusters may reflect batch.")

    # ---------------------------------------------------------------- neighbors / umap / leiden
    job.step(55, "Building neighbour graph and UMAP")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, use_rep=use_rep, random_state=0)
    sc.tl.umap(adata, random_state=0)
    job.step(68, "Leiden clustering")
    res = float(P["resolution"])
    res_grid = sorted({0.25, 0.5, 1.0, res})
    for r in res_grid:
        sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}", flavor="igraph", n_iterations=2, random_state=0)
    adata.obs["leiden"] = adata.obs[f"leiden_{res}"]
    K = adata.obs["leiden"].nunique()

    # ---------------------------------------------------------------- markers
    job.step(75, "Finding marker genes (Wilcoxon)")
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", layer="log1p", use_raw=False, pts=True)
    ann, S = annotate_clusters(adata)
    adata.obs["cell_type"] = adata.obs["leiden"].map({c: f"{c}: {v[0]}" if v[0] != "Unassigned" else f"{c}: Unassigned" for c, v in ann.items()}).astype("category")
    adata.obs["cell_type_simple"] = adata.obs["leiden"].map({c: re.sub(r" \(\d+\)$", "", v[0]) for c, v in ann.items()}).astype("category")
    cl_colors = [PALETTE[i % 20] for i in range(K)]
    adata.uns["leiden_colors"] = cl_colors
    adata.uns["cell_type_colors"] = cl_colors

    # ---------------------------------------------------------------- overview tiles + flags
    med_g = np.median(adata.obs["n_genes_by_counts"]); med_u = np.median(adata.obs["total_counts"])
    job.tile(commas(N0), "barcodes loaded")
    job.tile(commas(adata.n_obs), "cells after QC")
    job.tile(commas(med_g), "median genes / cell")
    job.tile(commas(med_u), "median UMIs / cell")
    job.tile(commas(adata.n_vars), "genes in ≥3 cells")
    job.tile(str(K), f"Leiden clusters (res {res})")
    for n in notes:
        job.flag("warn", n[5:]) if n.startswith("WARN:") else job.flag("info", n)
    frac_rm = (removed + n_doublets) / N0
    job.flag("warn" if frac_rm > 0.3 else "ok",
             f"QC removed <b>{commas(removed)}</b> barcodes ({pct(removed, N0)}): {int(low.sum())} with < {min_g} genes, {int(high.sum())} with > {max_g:,} genes, "
             + (f"{int(mt_bad.sum())} over {max_mt:.1f}% mito, " if has_mt else "no mito filter (no MT genes), ")
             + f"{int(top_bad.sum())} dominated by a few genes"
             + (f"; Scrublet removed <b>{n_doublets}</b> doublets." if P["doublets"] else ".")
             + (" Over 30% removed is high — check for a failed capture or ambient RNA." if frac_rm > 0.3 else ""))
    if not has_mt:
        job.flag("warn", "No mitochondrial genes (MT-/mt- prefix) were found, so damaged cells couldn't be flagged. Convert Ensembl IDs to gene symbols if needed.")
    job.flag("warn" if med_g < 500 else "ok",
             f"Median <b>{commas(med_g)}</b> genes and <b>{commas(med_u)}</b> UMIs per cell after QC"
             + (" — shallow; rare populations may not resolve." if med_g < 500 else " — typical depth for droplet-based data."))
    named = [(c, v) for c, v in ann.items() if v[0] != "Unassigned"]
    job.flag("ok", f"Leiden found <b>{K}</b> clusters. Marker-based labels: " +
             (", ".join(f"{c} = {v[0]}" for c, v in named) if named else "none matched the reference panel confidently") + ".")
    unl = [c for c, v in ann.items() if v[0] == "Unassigned"]
    if unl:
        job.flag("warn", f"Clusters {', '.join(unl)} had no confident label — possibly a cell type outside the panel, residual doublets or stressed cells. Check their markers below.")
    if integrated:
        job.flag("info", f"{nb} batches (<b>{batch_key}</b>) were integrated with Harmony before clustering.")
    elif batch_key and nb > 1:
        job.flag("warn", f"{nb} batches found in <b>{batch_key}</b> but not integrated; check the UMAP coloured by batch for batch-driven clusters.")

    # mito / cc driven clusters
    per_cl = adata.obs.groupby("leiden", observed=True).agg(mt=("pct_counts_mt", "median"), ng=("n_genes_by_counts", "median"))
    hi_mt = per_cl.index[per_cl["mt"] > np.median(adata.obs["pct_counts_mt"]) + 2.5 * mad(adata.obs["pct_counts_mt"].values)].tolist()
    if hi_mt:
        job.flag("warn", f"Cluster(s) {', '.join(hi_mt)} have much higher mitochondrial fraction than the rest — they may be stressed/dying cells rather than a distinct type.")

    # ---------------------------------------------------------------- embedding figures
    job.step(82, "Plotting clusters")
    job.section("clusters", "Step 3", "Clustering & embedding",
                f"Leiden community detection on a 15-nearest-neighbour graph ({n_pcs} PCs{', Harmony-corrected' if integrated else ''}). "
                "UMAP is for visualisation only — trust which cells group together, not distances between islands.")
    job.widget("explorer", id="explorer", title="Interactive UMAP explorer", data="embedding.json")

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
    sc.pl.umap(adata, color="leiden", ax=axs[0], show=False, legend_loc="on data", legend_fontoutline=2, frameon=False, title=f"Leiden (resolution {res})", size=max(4, 120000 / adata.n_obs / 3))
    sc.pl.umap(adata, color="cell_type", ax=axs[1], show=False, frameon=False, title="Suggested cell types", size=max(4, 120000 / adata.n_obs / 3))
    fig.tight_layout()
    sizes = adata.obs["leiden"].value_counts()
    small = [c for c, n in sizes.items() if n / adata.n_obs < 0.02]
    job.figure("umap_clusters", "UMAP coloured by cluster and label", fig, wide=True,
               how="Each dot is a cell; colours are Leiden clusters (left) and the reference-panel label for each cluster (right).",
               yours=f"Largest cluster is <b>{sizes.index[0]}</b> ({pct(sizes.iloc[0], adata.n_obs)} of cells)."
                     + (f" Rare clusters (<2%): {', '.join(small)}." if small else ""))

    qc_keys = ["n_genes_by_counts", "total_counts", "pct_counts_mt"] + (["doublet_score"] if "doublet_score" in adata.obs else []) + (["phase"] if has_cc else []) + ([batch_key] if batch_key and nb > 1 else [])
    sc.pl.umap(adata, color=qc_keys, ncols=3, frameon=False, show=False, cmap="viridis", size=max(3, 80000 / adata.n_obs / 3), wspace=0.35)
    fig = plt.gcf()
    cc_txt = ""
    if has_cc:
        ph = adata.obs["phase"].value_counts(normalize=True)
        cc_txt = f" Cell-cycle: {100*ph.get('G1',0):.0f}% G1, {100*ph.get('S',0):.0f}% S, {100*ph.get('G2M',0):.0f}% G2M."
    job.figure("umap_qc", "UMAP coloured by QC covariates", fig, wide=True,
               how="If a technical metric lights up one island, that cluster may be an artefact (low-quality cells, doublets, batch) rather than biology.",
               yours=(f"Clusters with elevated mito: {', '.join(hi_mt)}." if hi_mt else "No cluster is dominated by a QC metric.") + cc_txt)

    from . import sc_extra
    sil = {r: v for r, v in zip(res_grid, sc_extra.silhouette(adata, [f"leiden_{r}" for r in res_grid], use_rep).values())}
    fig, axs = plt.subplots(1, len(res_grid), figsize=(3.3 * len(res_grid), 3.3))
    for ax, r in zip(np.atleast_1d(axs), res_grid):
        k = adata.obs[f"leiden_{r}"].nunique()
        adata.uns[f"leiden_{r}_colors"] = [PALETTE[i % 20] for i in range(k)]
        sc.pl.umap(adata, color=f"leiden_{r}", ax=ax, show=False, frameon=False, legend_loc="on data", legend_fontsize=7, legend_fontoutline=1.5, title=f"res {r} · {k} clusters", size=max(2, 40000 / adata.n_obs / 3))
    fig.tight_layout()
    job.figure("umap_resolution", "Clustering resolution sweep", fig, wide=True,
               how="Higher resolution splits the graph into more clusters. Stable clusters that persist across resolutions are more trustworthy; clusters that only appear at high resolution may be over-splitting.",
               yours=" · ".join(f"res {r}: <b>{adata.obs[f'leiden_{r}'].nunique()}</b>" for r in res_grid) + f". Resolution {res} is used below; change it in settings."
                     + (" Silhouette (cluster separation, −1…1): " + ", ".join(f"res {r}: {v:.2f}" for r, v in sil.items()) + " — higher is cleaner, but biology may still justify finer clusters." if sil else ""))

    # composition
    if batch_key and nb > 1:
        comp = pd.crosstab(adata.obs[batch_key], adata.obs["cell_type"], normalize="index") * 100
        fig, ax = plt.subplots(figsize=(min(12, 3 + nb * 0.6), 4))
        bottom = np.zeros(len(comp))
        for i, c in enumerate(comp.columns):
            ax.bar(comp.index.astype(str), comp[c], bottom=bottom, color=cl_colors[i % len(cl_colors)], label=c, width=0.8)
            bottom += comp[c].values
        ax.set_ylabel("% of cells"); ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
        job.figure("composition", f"Cell-type composition by {batch_key}", fig, wide=True,
                   how="Proportion of each cluster per sample. A cluster made almost entirely of one sample suggests a batch effect or a sample-specific population.",
                   yours="")
    else:
        fig, ax = plt.subplots(figsize=(5, max(2.5, K * 0.28)))
        cts = adata.obs["cell_type"].value_counts().reindex(adata.obs["cell_type"].cat.categories)
        ax.barh(range(K), cts.values, color=cl_colors)
        ax.set_yticks(range(K)); ax.set_yticklabels(cts.index, fontsize=8); ax.invert_yaxis()
        for i, v in enumerate(cts.values):
            ax.text(v, i, f" {v:,}", va="center", fontsize=7.5)
        ax.set_xlabel("cells")
        job.figure("composition", "Cells per cluster", fig,
                   how="Cluster sizes. Very small clusters deserve scrutiny: rare cell types and technical artefacts both look like this.",
                   yours=f"{K} clusters ranging from {cts.min():,} to {cts.max():,} cells.")

    # ---------------------------------------------------------------- marker figures
    job.step(88, "Plotting marker genes")
    job.section("markers", "Step 4", "Marker genes & annotation",
                "Wilcoxon rank-sum test of each cluster against all other cells on log-normalized expression (Benjamini–Hochberg). "
                "Labels come from scoring clusters against a panel of canonical markers — treat them as hypotheses to confirm.")
    sc.tl.dendrogram(adata, "leiden", use_rep=use_rep, n_pcs=n_pcs)
    fig = _fig_of(sc.pl.rank_genes_groups_dotplot(adata, n_genes=4, layer="log1p", use_raw=False, standard_scale="var", show=False, dendrogram=True, colorbar_title="scaled mean\nexpression"))
    mk = sc.get.rank_genes_groups_df(adata, None)
    top_by = {c: mk[(mk.group == c) & (mk.logfoldchanges > 0.5)].head(3).names.tolist() for c in adata.obs["leiden"].cat.categories}
    job.figure("dotplot_markers", "Top markers per cluster (dot plot)", fig, wide=True,
               sub="dot size = fraction of cells expressing · colour = mean expression scaled per gene",
               how="The standard Seurat/Scanpy marker figure. A good marker is a big dark dot in its own cluster's row and small pale dots elsewhere. Clusters are ordered by a dendrogram so similar ones sit together.",
               yours=" · ".join(f"<b>{c}</b>: {', '.join(g)}" for c, g in list(top_by.items())[:8]))

    panel = {}
    for ct, genes in CELL_TYPE_MARKERS.items():
        present = [upper[g] for g in genes if g in upper]
        if any(v[0].startswith(ct) for v in ann.values()) and present:
            panel[ct] = present[:4]
    if panel:
        seen = set()
        panel = {k: [g for g in v if not (g in seen or seen.add(g))] for k, v in panel.items()}
        panel = {k: v for k, v in panel.items() if v}
        fig = _fig_of(sc.pl.stacked_violin(adata, panel, groupby="cell_type", layer="log1p", use_raw=False, show=False, dendrogram=False, standard_scale="var", swap_axes=False))
        job.figure("violin_canonical", "Canonical markers behind each label", fig, wide=True,
                   how="Expression of the reference markers used for annotation, grouped by labelled cluster. Each label should be supported by several of its markers, not just one.",
                   yours="Evidence per cluster is listed in the annotation table below.")

    try:
        fig = _fig_of(sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, layer="log1p", use_raw=False, standard_scale="var", show=False, dendrogram=True, show_gene_labels=True, cmap="magma", figsize=(11, 6)))
        job.figure("heatmap_markers", "Marker heatmap (cells)", fig, wide=True,
                   how="Every column is a cell grouped by cluster; rows are the top 5 markers per cluster. Sharp blocks along the diagonal mean clean, well-separated clusters; smeared blocks suggest clusters that should be merged.",
                   yours="")
    except Exception:  # noqa: BLE001
        log_exc("optional figure")

    rows = []
    for c, (label, score, ev) in ann.items():
        sub = mk[mk.group == c].head(40)
        sub = sub[(sub.logfoldchanges > 0.5) & (sub.pvals_adj < 0.05)]
        rows.append([c, label, f"{score:.2f}", ", ".join(ev[:5]) or "—", ", ".join(sub.names.head(6)), commas(sizes.get(c, 0))])
    job.table("annotation", "Cluster annotation", ["Cluster", "Suggested label", "Score", "Supporting markers", "Top DE genes", "Cells"], rows,
              note="Score = mean z-score of the reference markers in that cluster (≥0.8 to assign). Confirm labels with CellTypist, Azimuth, or literature markers (CellMarker 2.0, PanglaoDB).")
    mk_out = mk.merge(pd.DataFrame(adata.uns["rank_genes_groups"]["pts"]).rename_axis("names").reset_index().melt("names", var_name="group", value_name="pct_in"), on=["names", "group"], how="left")
    mk_out.to_csv(job.root / "markers.csv", index=False)
    t_rows = []
    for c in adata.obs["leiden"].cat.categories:
        s2 = mk_out[(mk_out.group == c) & (mk_out.logfoldchanges > 0)].head(8)
        for _, r in s2.iterrows():
            t_rows.append([c, r.names, f"{r.logfoldchanges:.2f}", f"{100*r.pct_in:.0f}%" if pd.notna(r.pct_in) else "", f"{r.scores:.1f}", f"{r.pvals_adj:.1e}"])
    job.table("markers", "Top 8 markers per cluster", ["Cluster", "Gene", "log2FC", "% in cluster", "z-score", "adj. p"], t_rows,
              note="Full table (all genes, all clusters) in markers.csv.", csv="markers.csv")

    # ---------------------------------------------------------------- publication panels + automated annotation
    job.step(89, "Feature plots, split UMAPs, automated annotation")
    sc_extra.feature_plots(job, adata, "leiden", "cell_type_simple")
    if P["genes"]:
        _user_feature_plots(job, adata, P["genes"])
    ck_guess = P["condition_key"] or next((c for c in ("condition", "genotype", "treatment", "group", "disease", "status") if c in adata.obs), None)
    sk_guess = P["sample_key"] or batch_key
    sc_extra.umap_split(job, adata, ck_guess, sk_guess, "cell_type_simple")
    sc_extra.alluvial(job, adata, "leiden", "cell_type_simple", fid="alluvial_labels", title="Clusters → marker-panel labels")
    if ck_guess and ck_guess in adata.obs and adata.obs[ck_guess].nunique() <= 12:
        sc_extra.alluvial(job, adata, "cell_type_simple", ck_guess, fid="alluvial_condition", title=f"Cell types → {ck_guess}")
    if P["celltypist"]:
        job.step(90, "Automated annotation (CellTypist)")
        sc_extra.celltypist_annotate(job, adata, species, P["celltypist_model"], "leiden", "cell_type_simple")
    if P["pathways"]:
        job.step(91, "Pathway activity per population (PROGENy)")
        sc_extra.pathway_activity(job, adata, species, "cell_type_simple")
    if P["trajectory"]:
        job.step(92, "Trajectory: PAGA + diffusion pseudotime")
        sc_extra.trajectory(job, adata, "leiden", "cell_type_simple", P["root_cluster"], use_rep)
    if P["ccc"]:
        job.step(93, "Cell–cell communication (LIANA)")
        sc_extra.cell_communication(job, adata, species, "cell_type_simple", ck_guess if ck_guess in adata.obs else None)

    # ---------------------------------------------------------------- sample-level tests
    pb_info = None
    if P["pseudobulk"]:
        from . import pseudobulk as _pb
        sk, ck = _pb.detect_keys(adata.obs, P["sample_key"] or batch_key, P["condition_key"])
        if sk and ck:
            try:
                pb_info = _pb.run(job, adata, sk, ck, P["pb_covariates"] or [c for c in ("protocol", "digestion_protocol", "batch", "sex") if c in adata.obs and c not in (sk, ck)], P["reference"])
            except Exception as e:  # noqa: BLE001
                job.flag("warn", f"Sample-level tests failed ({type(e).__name__}: {e}).")
        else:
            job.flag("info", "No sample-level condition column was found, so composition and pseudobulk tests were skipped. Provide cell metadata with a sample column and a condition column (e.g. genotype) to enable them.")
    job.result["pseudobulk"] = pb_info
    if pb_info and P["gsea"] and (job.root / "pseudobulk_all_cells.csv").exists():
        try:
            from . import bulk_extra
            job.step(93, "GSEA on the pseudobulk comparison")
            r_all = pd.read_csv(job.root / "pseudobulk_all_cells.csv", index_col=0)
            bulk_extra.gsea(job, r_all, species, pb_info.get("alt", "alt"), pb_info.get("ref", "ref"))
        except Exception:  # noqa: BLE001
            job.flag("info", f"GSEA skipped ({log_exc('sc gsea')}).")

    # ---------------------------------------------------------------- outputs
    job.step(94, "Saving results")
    # per-sample sheet (columns that are constant within a sample) so the generated script can reproduce the run
    sk_out = P["sample_key"] or batch_key or ("sample" if "sample" in adata.obs else None)
    if sk_out and sk_out in adata.obs:
        g = adata.obs.groupby(sk_out, observed=True)
        const = [c for c in adata.obs.columns if c != sk_out and adata.obs[c].dtype.kind in "OUb" or str(adata.obs[c].dtype) == "category"]
        skip = {"qc_pass", "predicted_doublet", "leiden", "cell_type", "cell_type_simple", "phase", "doublet_score"}
        const = [c for c in const if c != sk_out and c not in skip and not c.startswith("leiden") and (g[c].nunique() <= 1).all()]
        if const:
            g[const].first().to_csv(job.root / "sample_metadata.csv")
    emb = adata.obsm["X_umap"]
    rng = np.random.default_rng(1)
    idx = np.sort(rng.choice(adata.n_obs, size=min(adata.n_obs, 30000), replace=False))
    json.dump({
        "x": np.round(emb[idx, 0], 3).tolist(), "y": np.round(emb[idx, 1], 3).tolist(),
        "cluster": adata.obs["leiden"].cat.codes.values[idx].tolist(),
        "clusters": [f"{c}: {ann[c][0]}" for c in adata.obs["leiden"].cat.categories],
        "colors": cl_colors,
        "numeric": {k: np.round(adata.obs[k].values[idx].astype(float), 3).tolist() for k in ["n_genes_by_counts", "total_counts", "pct_counts_mt"] + (["doublet_score"] if "doublet_score" in adata.obs else [])},
        "categorical": {k: {"levels": list(map(str, adata.obs[k].astype("category").cat.categories)), "codes": adata.obs[k].astype("category").cat.codes.values[idx].tolist()} for k in ([batch_key] if batch_key and nb > 1 else []) + (["phase"] if has_cc else [])},
        "index": idx.tolist(),
        "genes": adata.var_names[np.argsort(-np.asarray(adata.layers["log1p"].mean(0)).ravel())[:4000]].tolist(),
        "top_markers": {c: g for c, g in top_by.items()},
    }, open(job.root / "embedding.json", "w"))

    out = adata.copy()
    out.X = out.layers.pop("log1p")
    for k in list(out.obs.columns):
        if out.obs[k].dtype == bool:
            out.obs[k] = out.obs[k].astype(str)
    out.uns.pop("rank_genes_groups", None)
    out.write_h5ad(job.root / "analyzed.h5ad", compression="gzip")
    job.download("Processed AnnData (.h5ad)", "analyzed.h5ad")
    job.download("Marker genes (.csv)", "markers.csv")
    if (job.root / "sample_metadata.csv").exists():
        job.download("Sample sheet (.csv)", "sample_metadata.csv")
    job.download("Cell metadata: QC, cluster, cell type (.csv)", "cell_metadata.csv")
    adata.obs.to_csv(job.root / "cell_metadata.csv")

    job.result["summary"] = {"cells": adata.n_obs, "clusters": K, "species": species}
    job.method("Input", f"{N0:,} barcodes × {G0:,} genes. Detected species: {species}.")
    job.method("Quality control", f"sc.pp.calculate_qc_metrics; kept cells with {min_g}–{max_g:,} genes, ≤{max_mt:.1f}% mitochondrial counts and top-20-gene share ≤ median + 5 MAD. Genes kept if detected in ≥{P['min_cells']} cells.")
    if P["doublets"] is True:
        job.method("Doublets", "Scrublet (sc.pp.scrublet), predicted doublets removed.")
    job.method("Normalization", "sc.pp.normalize_total(target_sum=1e4) + sc.pp.log1p — equivalent to Seurat LogNormalize.")
    job.method("Feature selection", f"sc.pp.highly_variable_genes, flavor={hvg_method.split()[0]}, n_top_genes={n_hvg}" + (f", batch_key={batch_key}" if nb > 1 else "") + ".")
    if has_cc:
        job.method("Cell cycle", "sc.tl.score_genes_cell_cycle with Tirosh et al. (2016) S and G2/M gene lists.")
    job.method("PCA", f"Scaled to unit variance (clipped at 10); ARPACK PCA with {n_comps} components; {n_pcs} used (elbow + 5).")
    if integrated:
        job.method("Integration", f"Harmony (harmonypy) on PCA, batch = {batch_key}.")
    job.method("Graph & UMAP", "sc.pp.neighbors(n_neighbors=15) + sc.tl.umap (min_dist 0.5).")
    job.method("Clustering", f"Leiden (igraph), resolution {res}; sweep at {', '.join(map(str, res_grid))}.")
    job.method("Markers", "sc.tl.rank_genes_groups(method='wilcoxon') vs rest on log-normalized counts, BH-adjusted.")
    job.method("Annotation", "Cluster mean z-scores of a curated panel of ~30 canonical cell-type marker sets, weighted by detection rate. Heuristic; validate with reference mapping.")
    job.result["versions"] = {p: version(p) for p in ["scanpy", "anndata", "numpy", "scipy", "leidenalg", "umap-learn"]}
