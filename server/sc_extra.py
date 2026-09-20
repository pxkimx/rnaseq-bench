"""Expert-level add-ons for the single-cell pipeline. Every function adds to the job and returns True on
success; anything that needs a resource that cannot be fetched is skipped with a note, never fatal.

- feature_plots          static UMAP feature plots of the top markers (publication Figure 2 panels)
- umap_split             UMAP split by condition / sample (same axes, so shifts are comparable)
- alluvial               cluster -> label (or condition) flow diagram
- celltypist_annotate    automated reference-based cell-type labels (CellTypist, majority vote per cluster)
- pathway_activity       PROGENy pathway activity per cluster (decoupler)
- trajectory             PAGA graph + diffusion pseudotime from a root cluster
- cell_communication     ligand–receptor analysis between cell types (LIANA rank_aggregate)
- silhouette             clustering quality per resolution
"""
from __future__ import annotations

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from . import resources
from .common import PALETTE, Job, commas, log_exc


def _lognorm(adata) -> "sc.AnnData":
    """A light AnnData with log-normalized counts in X (the pipeline keeps raw counts in layers['counts'])."""
    if "log1p" in adata.layers:
        b = sc.AnnData(adata.layers["log1p"].copy(), obs=adata.obs.copy(), var=adata.var[[]].copy())
    else:
        b = sc.AnnData(adata.layers["counts"].copy(), obs=adata.obs.copy(), var=adata.var[[]].copy())
        sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    for k in ("X_umap", "X_pca"):
        if k in adata.obsm:
            b.obsm[k] = adata.obsm[k]
    return b


def _top_markers(adata, cluster_key: str, n: int = 2) -> dict[str, list[str]]:
    out = {}
    try:
        mk = sc.get.rank_genes_groups_df(adata, group=None)
        mk = mk[(mk.logfoldchanges > 1) & (mk.pvals_adj < 0.05)]
        for g, d in mk.groupby("group", observed=True):
            out[str(g)] = d.sort_values("scores", ascending=False).names.head(n).tolist()
    except Exception:  # noqa: BLE001
        log_exc("top markers")
    return out


# ---------------------------------------------------------------- figure-2 style panels
def feature_plots(job: Job, adata, cluster_key: str = "leiden", label_key: str = "cell_type_simple") -> bool:
    try:
        top = _top_markers(adata, cluster_key, 2)
        genes, seen = [], set()
        for c in adata.obs[cluster_key].cat.categories:
            for g in top.get(str(c), []):
                if g not in seen and len(genes) < 16:
                    genes.append(g); seen.add(g)
        if not genes:
            return False
        b = _lognorm(adata)
        n = len(genes); cols = 4; rows = int(np.ceil(n / cols))
        fig, axs = plt.subplots(rows, cols, figsize=(3.2 * cols, 2.9 * rows), squeeze=False)
        xy = b.obsm["X_umap"]
        order = np.random.default_rng(0).permutation(b.n_obs)
        for ax, g in zip(axs.ravel(), genes):
            v = b[:, g].X
            v = np.asarray(v.toarray() if hasattr(v, "toarray") else v).ravel()
            o = np.argsort(v[order], kind="stable"); idx = order[o]          # high values drawn last
            sca = ax.scatter(xy[idx, 0], xy[idx, 1], c=v[idx], s=2, cmap="viridis", linewidths=0, rasterized=True, vmax=np.percentile(v, 99) or 1)
            lab = ""
            for c, gl in top.items():
                if g in gl:
                    cl = adata.obs.loc[adata.obs[cluster_key].astype(str) == c, label_key]
                    lab = f" · cluster {c}" + (f" ({cl.iloc[0]})" if len(cl) and label_key in adata.obs else "")
                    break
            ax.set_title(g + lab, fontsize=8); ax.set_xticks([]); ax.set_yticks([])
            plt.colorbar(sca, ax=ax, fraction=0.04, pad=0.02).ax.tick_params(labelsize=6)
        for ax in axs.ravel()[n:]:
            ax.axis("off")
        fig.tight_layout()
        job.figure("feature_plots", "Feature plots: top markers on the UMAP", fig, wide=True,
                   how="Each panel colours cells by log-normalized expression of one gene (the strongest markers per cluster, 99th percentile clipped). "
                       "A marker should light up one island and be dark elsewhere; a gene glowing everywhere is not cluster-specific.",
                   yours=f"{n} genes shown: {', '.join(genes)}.")
        return True
    except Exception:  # noqa: BLE001
        log_exc("feature plots"); return False


def umap_split(job: Job, adata, condition_key: str | None, sample_key: str | None, label_key: str = "cell_type_simple") -> bool:
    try:
        key = condition_key if condition_key and condition_key in adata.obs else (sample_key if sample_key and sample_key in adata.obs else None)
        if key is None or adata.obs[key].nunique() < 2 or adata.obs[key].nunique() > 12:
            return False
        groups = [str(g) for g in pd.unique(adata.obs[key].astype(str))]
        labels = adata.obs[label_key].astype(str) if label_key in adata.obs else adata.obs["leiden"].astype(str)
        cats = sorted(labels.unique(), key=str)
        col = {c: PALETTE[i % 20] for i, c in enumerate(cats)}
        xy = adata.obsm["X_umap"]
        n = len(groups); cols = min(3, n); rows = int(np.ceil(n / cols))
        fig, axs = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.9 * rows), squeeze=False)
        for ax, g in zip(axs.ravel(), groups):
            ax.scatter(xy[:, 0], xy[:, 1], s=1.5, c="#e6e6e6", linewidths=0, rasterized=True)
            m = (adata.obs[key].astype(str) == g).values
            ax.scatter(xy[m, 0], xy[m, 1], s=2.5, c=[col[l] for l in labels[m]], linewidths=0, rasterized=True)
            ax.set_title(f"{key} = {g}  ({commas(int(m.sum()))} cells)", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        for ax in axs.ravel()[n:]:
            ax.axis("off")
        handles = [plt.Line2D([], [], marker="o", ls="", color=col[c], label=c, markersize=5) for c in cats]
        fig.legend(handles=handles, loc="lower center", ncol=min(6, len(cats)), fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.02))
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        job.figure("umap_split", f"UMAP split by {key}", fig, wide=True,
                   how="Same embedding, one panel per group, all other cells in grey. Populations present in one panel and empty in another are condition-specific; "
                       "the composition test below puts numbers and a p-value on such shifts.",
                   yours=f"{n} groups of <b>{key}</b>.")
        return True
    except Exception:  # noqa: BLE001
        log_exc("umap split"); return False


def alluvial(job: Job, adata, left: str, right: str, fid="alluvial", title=None) -> bool:
    """Flow diagram between two categorical obs columns (cluster -> cell type, or cell type -> condition)."""
    try:
        if left not in adata.obs or right not in adata.obs:
            return False
        ct = pd.crosstab(adata.obs[left].astype(str), adata.obs[right].astype(str))
        if ct.shape[0] > 40 or ct.shape[1] > 40:
            return False
        ct = ct.loc[ct.sum(1).sort_values(ascending=False).index, ct.sum(0).sort_values(ascending=False).index]
        L, R = ct.index.tolist(), ct.columns.tolist()
        total = ct.values.sum(); gap = 0.012 * total
        def spans(sizes):
            out, y = {}, 0.0
            for k, v in sizes.items():
                out[k] = (y, y + v); y += v + gap
            return out, y
        ls, hl = spans(ct.sum(1)); rs, hr = spans(ct.sum(0))
        H = max(hl, hr)
        fig, ax = plt.subplots(figsize=(7.5, max(4, 0.28 * max(len(L), len(R)) + 2)))
        colL = {k: PALETTE[i % 20] for i, k in enumerate(L)}
        for k, (a, b) in ls.items():
            ax.add_patch(plt.Rectangle((0, a), 0.06, b - a, color=colL[k])); ax.text(-0.02, (a + b) / 2, k, ha="right", va="center", fontsize=7.5)
        for k, (a, b) in rs.items():
            ax.add_patch(plt.Rectangle((0.94, a), 0.06, b - a, color="#9aa7ad")); ax.text(1.02, (a + b) / 2, k, ha="left", va="center", fontsize=7.5)
        lpos = {k: a for k, (a, _) in ls.items()}; rpos = {k: a for k, (a, _) in rs.items()}
        xs = np.linspace(0.06, 0.94, 60)
        for l in L:
            for r in R:
                v = ct.loc[l, r]
                if v == 0:
                    continue
                y0, y1 = lpos[l], rpos[r]
                t = (xs - 0.06) / 0.88; sm = t * t * (3 - 2 * t)
                lo = y0 + (y1 - y0) * sm; hi = lo + v
                ax.fill_between(xs, lo, hi, color=colL[l], alpha=0.45, linewidth=0)
                lpos[l] += v; rpos[r] += v
        ax.set_xlim(-0.35, 1.35); ax.set_ylim(-gap, H); ax.axis("off")
        ax.set_title(title or f"{left} → {right}", fontsize=10)
        fig.tight_layout()
        job.figure(fid, title or f"Flow: {left} → {right}", fig,
                   how="Ribbon width = number of cells. Reading left to right shows how each cluster distributes over the categories on the right.",
                   yours=f"{len(L)} × {len(R)} categories, {commas(int(total))} cells.")
        return True
    except Exception:  # noqa: BLE001
        log_exc("alluvial"); return False


# ---------------------------------------------------------------- CellTypist
DEFAULT_MODELS = {"human": "Immune_All_Low.pkl", "mouse": "Immune_All_Low.pkl"}


def celltypist_annotate(job: Job, adata, species: str, model_name: str | None, cluster_key: str = "leiden", marker_key: str = "cell_type_simple") -> bool:
    try:
        import celltypist
        name = model_name or DEFAULT_MODELS.get(species, "Immune_All_Low.pkl")
        path = resources.celltypist_model(name)
        if path is None:
            job.flag("info", f"CellTypist annotation skipped — model {name} could not be downloaded (needs internet once).")
            return False
        b = _lognorm(adata)
        if species == "mouse":
            b.var_names = b.var_names.str.upper()          # human models on mouse data: match by upper-cased symbol
            b.var_names_make_unique()
        model = celltypist.models.Model.load(path)
        pred = celltypist.annotate(b, model=model, majority_voting=True, over_clustering=adata.obs[cluster_key].astype(str).values)
        lab = pred.predicted_labels
        adata.obs["celltypist"] = lab["majority_voting"].astype(str).values
        adata.obs["celltypist_cell"] = lab["predicted_labels"].astype(str).values
        conf = pred.probability_matrix.max(axis=1).values if hasattr(pred, "probability_matrix") else np.full(adata.n_obs, np.nan)
        adata.obs["celltypist_conf"] = conf
        # UMAP coloured by CellTypist label
        xy = adata.obsm["X_umap"]
        cats = pd.Series(adata.obs["celltypist"]).value_counts().index.tolist()
        col = {c: PALETTE[i % 20] for i, c in enumerate(cats)}
        fig, ax = plt.subplots(figsize=(6.4, 5.2))
        ax.scatter(xy[:, 0], xy[:, 1], s=2, c=[col[c] for c in adata.obs["celltypist"]], linewidths=0, rasterized=True)
        for c in cats:
            m = (adata.obs["celltypist"] == c).values
            if m.sum() > 20:
                ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), c, fontsize=6.5, ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(f"CellTypist ({name.replace('.pkl', '')}), majority vote per cluster", fontsize=9)
        fig.tight_layout()
        # agreement table
        ct = pd.crosstab(adata.obs[cluster_key].astype(str), adata.obs["celltypist"])
        rows = []
        for c in adata.obs[cluster_key].cat.categories:
            c = str(c)
            if c not in ct.index:
                continue
            r = ct.loc[c]; best = r.idxmax(); frac = r.max() / r.sum()
            mk_lab = adata.obs.loc[adata.obs[cluster_key].astype(str) == c, marker_key].iloc[0] if marker_key in adata.obs else ""
            mc = float(np.nanmean(conf[(adata.obs[cluster_key].astype(str) == c).values]))
            rows.append([c, str(mk_lab), best, f"{100*frac:.0f}%", f"{mc:.2f}", commas(int(r.sum()))])
        job.figure("umap_celltypist", "Automated annotation (CellTypist)", fig,
                   how="A logistic-regression classifier trained on a curated reference atlas assigns each cell a type; the majority label per cluster is shown. "
                       "It is independent of the marker panel, so agreement between the two columns in the table is strong evidence; disagreement flags clusters to check by hand.",
                   yours=f"{len(cats)} labels assigned. " + "; ".join(f"cluster {r[0]}: {r[2]}" for r in rows[:6]) + ("…" if len(rows) > 6 else "."))
        job.table("celltypist_table", "Marker-panel label vs CellTypist label per cluster",
                  ["Cluster", "Marker panel", "CellTypist (majority)", "Agreement within cluster", "Mean confidence", "Cells"], rows,
                  note=f"Model: {name}. Human immune model applied to {'mouse (symbols upper-cased) — non-immune and mouse-specific types will be mislabelled; choose a tissue model in the settings if one exists' if species == 'mouse' else 'human data'}.")
        job.method("CellTypist", f"celltypist.annotate with model {name}, majority voting over Leiden clusters.")
        return True
    except Exception:  # noqa: BLE001
        job.flag("info", f"CellTypist annotation skipped ({log_exc('celltypist')}).")
        return False


# ---------------------------------------------------------------- pathway activity per cluster
def pathway_activity(job: Job, adata, species: str, group_key: str = "cell_type_simple") -> bool:
    try:
        import decoupler as dc
        net = resources.progeny(species)
        if net is None:
            job.flag("info", "Pathway activity per cluster skipped — PROGENy could not be downloaded (needs internet once).")
            return False
        b = _lognorm(adata)
        key = group_key if group_key in adata.obs else "leiden"
        groups = adata.obs[key].astype(str)
        X = b.X
        means = pd.DataFrame({g: np.asarray(X[(groups == g).values].mean(axis=0)).ravel() for g in pd.unique(groups)}, index=b.var_names).T
        means = means - means.mean(axis=0)
        means.columns = means.columns.astype(str)
        if species == "mouse":
            means.columns = means.columns.str.upper()
        means = means.loc[:, ~means.columns.duplicated()]
        es, pv = dc.mt.ulm(data=means, net=net, tmin=5)
        es.to_csv(job.root / "pathway_activity_by_cluster.csv")
        job.download("Pathway activity by cluster (.csv)", "pathway_activity_by_cluster.csv")
        fig, ax = plt.subplots(figsize=(max(5, 0.5 * es.shape[0] + 3), 4.8))
        lim = np.nanpercentile(np.abs(es.values), 98) or 1
        im = ax.imshow(es.T.values, aspect="auto", cmap="RdBu_r", vmin=-lim, vmax=lim)
        ax.set_yticks(range(es.shape[1])); ax.set_yticklabels(es.columns, fontsize=8)
        ax.set_xticks(range(es.shape[0])); ax.set_xticklabels(es.index, rotation=60, ha="right", fontsize=8)
        for i in range(es.shape[0]):
            for j in range(es.shape[1]):
                if pv.iloc[i, j] < 0.01:
                    ax.text(i, j, "•", ha="center", va="center", fontsize=7, color="black")
        plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="activity (ULM)")
        ax.set_title("PROGENy pathway activity per population", fontsize=10)
        fig.tight_layout()
        hi = [(pw, es[pw].idxmax()) for pw in es.columns if pv.loc[es[pw].idxmax(), pw] < 0.01]
        job.figure("pathway_activity", "Signalling pathway activity per population", fig, wide=True,
                   how="Mean expression per population (centred across populations) scored against PROGENy's pathway-responsive genes. Red = the pathway's "
                       "downstream programme is up in that population relative to the others; • = p < 0.01. Footprint-based, so it reflects signalling output, not receptor expression.",
                   yours="; ".join(f"<b>{pw}</b> highest in {g}" for pw, g in hi[:6]) or "No population shows a significant pathway signature.")
        job.method("Pathway activity", "decoupler ULM with PROGENy (top 500 responsive genes per pathway) on population-mean log-normalized expression.")
        return True
    except Exception:  # noqa: BLE001
        job.flag("info", f"Pathway activity skipped ({log_exc('pathway activity')}).")
        return False


# ---------------------------------------------------------------- trajectory
def trajectory(job: Job, adata, cluster_key: str = "leiden", label_key: str = "cell_type_simple", root: str | None = None, use_rep: str = "X_pca") -> bool:
    try:
        if adata.n_obs < 200 or adata.obs[cluster_key].nunique() < 3:
            return False
        sc.tl.paga(adata, groups=cluster_key)
        # root: the user's choice, else the most proliferative cluster (progenitors usually cycle), else cluster 0
        cats = [str(c) for c in adata.obs[cluster_key].cat.categories]
        if root is not None and str(root) in cats:
            root_c = str(root); why = "chosen by you"
        elif "phase" in adata.obs:
            cyc = adata.obs.groupby(cluster_key, observed=True)["phase"].apply(lambda s: (s != "G1").mean())
            root_c = str(cyc.idxmax()); why = f"the most proliferative cluster ({100*cyc.max():.0f}% cycling)"
        else:
            root_c = cats[0]; why = "default (cluster 0)"
        m = (adata.obs[cluster_key].astype(str) == root_c).values
        rep = adata.obsm[use_rep] if use_rep in adata.obsm else adata.obsm["X_pca"]
        centre = rep[m].mean(axis=0)
        adata.uns["iroot"] = int(np.where(m)[0][np.argmin(((rep[m] - centre) ** 2).sum(1))])
        sc.tl.diffmap(adata, n_comps=10)
        sc.tl.dpt(adata, n_dcs=10)
        pt = adata.obs["dpt_pseudotime"].replace([np.inf], np.nan)
        adata.obs["dpt_pseudotime"] = pt
        # figure: PAGA graph + pseudotime UMAP
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.8))
        sc.pl.paga(adata, ax=axs[0], show=False, threshold=0.05, node_size_scale=1.2, fontsize=7, frameon=False, title="PAGA connectivity between clusters")
        xy = adata.obsm["X_umap"]
        ok = pt.notna().values
        sca = axs[1].scatter(xy[ok, 0], xy[ok, 1], c=pt.values[ok], s=2, cmap="viridis", linewidths=0, rasterized=True)
        axs[1].scatter(xy[[adata.uns["iroot"]], 0], xy[[adata.uns["iroot"]], 1], s=60, marker="*", c="red", edgecolor="black", zorder=5)
        axs[1].set_xticks([]); axs[1].set_yticks([]); axs[1].set_title(f"Diffusion pseudotime from cluster {root_c} (★ root)", fontsize=9)
        plt.colorbar(sca, ax=axs[1], fraction=0.04, pad=0.02, label="pseudotime")
        fig.tight_layout()
        # ordering of clusters along pseudotime
        order = adata.obs.groupby(cluster_key, observed=True)["dpt_pseudotime"].median().sort_values()
        lab = adata.obs.groupby(cluster_key, observed=True)[label_key].first() if label_key in adata.obs else None
        rows = [[str(c), str(lab[c]) if lab is not None else "", f"{v:.2f}"] for c, v in order.items() if pd.notna(v)]
        job.section("trajectory", "Trajectory", "Cluster connectivity & pseudotime",
                    "PAGA (Wolf 2019) draws an edge between clusters whose cells are more connected in the neighbour graph than expected — a coarse map of which "
                    "populations are continuous with which. Diffusion pseudotime (Haghverdi 2016) then orders cells by graph distance from a root. "
                    "Only meaningful if the data actually contain a continuous process (differentiation, activation); discrete cell types will still get a pseudotime, so read it with care.")
        job.figure("paga_dpt", "PAGA graph and pseudotime", fig, wide=True,
                   how="Left: thick edges = strongly connected clusters (threshold 0.05). Right: cells coloured by pseudotime from the root. Root = " + why + ".",
                   yours="Cluster order along pseudotime: " + " → ".join(f"{r[0]}{(' ' + r[1]) if r[1] else ''}" for r in rows[:8]) + ". Change the root in the analysis settings if your biology starts elsewhere.")
        job.table("pseudotime_order", "Median pseudotime per cluster", ["Cluster", "Label", "Median pseudotime"], rows)
        job.method("Trajectory", f"sc.tl.paga on {cluster_key}; sc.tl.diffmap (10 components) and sc.tl.dpt with the root cell nearest the centroid of cluster {root_c}.")
        return True
    except Exception:  # noqa: BLE001
        job.flag("info", f"Trajectory analysis skipped ({log_exc('trajectory')}).")
        return False


# ---------------------------------------------------------------- cell–cell communication
def cell_communication(job: Job, adata, species: str, group_key: str = "cell_type_simple", condition_key: str | None = None) -> bool:
    try:
        import liana as li
        key = group_key if group_key in adata.obs else "leiden"
        if adata.obs[key].nunique() < 2:
            return False
        b = _lognorm(adata)
        b.obs[key] = b.obs[key].astype(str)
        resource = "mouseconsensus" if species == "mouse" else "consensus"
        li.mt.rank_aggregate(b, groupby=key, resource_name=resource, expr_prop=0.1, use_raw=False, n_perms=500, seed=0, verbose=False)
        r = b.uns["liana_res"].copy()
        r = r.sort_values("magnitude_rank")
        r.to_csv(job.root / "cell_communication.csv", index=False)
        job.download("Ligand–receptor interactions (.csv)", "cell_communication.csv")
        spec = r[r.specificity_rank < 0.05]
        top = spec.head(40) if len(spec) >= 10 else r.head(40)
        top = top.copy(); top["pair"] = top.ligand_complex + " → " + top.receptor_complex; top["st"] = top.source + " → " + top.target
        keep_st = top.st.value_counts().head(12).index          # the 12 busiest sender → receiver pairs keep the plot readable
        top = top[top.st.isin(keep_st)]
        pairs = top.pair.unique().tolist()[:25]; top = top[top.pair.isin(pairs)]; st = top.st.unique().tolist()
        # dot plot in the LIANA convention: size = specificity, colour = magnitude
        fig, ax = plt.subplots(figsize=(max(5, 0.5 * len(st) + 3.5), max(4, 0.3 * len(pairs) + 1.8)))
        for _, x in top.iterrows():
            i = pairs.index(x.pair); j = st.index(x.st)
            ax.scatter(j, i, s=30 + 100 * min(4, -np.log10(max(x.specificity_rank, 1e-4))), c=[1 - x.magnitude_rank], cmap="viridis", vmin=0, vmax=1, edgecolor="#333", linewidth=0.4)
        ax.scatter([], [], s=30 + 100 * 1, c="#888", label="specificity p≈0.1"); ax.scatter([], [], s=30 + 100 * 3, c="#888", label="p≈0.001")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=7, frameon=False, title="dot size", title_fontsize=7)
        import matplotlib as mpl
        fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1), cmap="viridis"), ax=ax, fraction=0.03, pad=0.12, label="magnitude (1 − rank)")
        ax.set_xticks(range(len(st))); ax.set_xticklabels(st, rotation=60, ha="right", fontsize=7.5)
        ax.set_yticks(range(len(pairs))); ax.set_yticklabels(pairs, fontsize=7.5)
        ax.set_xlim(-0.6, len(st) - 0.4); ax.set_ylim(len(pairs) - 0.4, -0.6)
        ax.set_title("Top ligand → receptor interactions (LIANA consensus)", fontsize=9)
        ax.grid(alpha=0.2)
        fig.tight_layout()
        job.section("ccc", "Cell–cell communication", "Ligand–receptor interactions between populations",
                    "LIANA scores every ligand–receptor pair in the consensus resource for each sender → receiver pair of populations, aggregating CellPhoneDB, "
                    "CellChat, NATMI, Connectome, SingleCellSignalR and logFC methods into a robust rank. Magnitude = how strongly ligand and receptor are expressed; "
                    "specificity = how particular the pair is to that sender/receiver combination. These are hypotheses about signalling, to be tested experimentally.")
        job.figure("ccc_dotplot", "Top interactions", fig, wide=True,
                   how="Colour = magnitude (brighter = higher combined ligand and receptor expression), dot size = specificity (bigger = more particular to that sender/receiver pair). Rows are ligand → receptor, columns sender → receiver; only the 12 busiest sender → receiver pairs are drawn.",
                   yours="; ".join(f"<b>{x.ligand_complex}→{x.receptor_complex}</b> ({x.source}→{x.target})" for _, x in top.head(4).iterrows()) + ".")
        rows = [[x.source, x.target, x.ligand_complex, x.receptor_complex, f"{x.magnitude_rank:.2g}", f"{x.specificity_rank:.2g}"] for _, x in top.iterrows()]
        job.table("ccc_table", "Top interactions", ["Sender", "Receiver", "Ligand", "Receptor", "Magnitude rank", "Specificity rank"], rows,
                  note="Lower rank = stronger. Full table in cell_communication.csv.", csv="cell_communication.csv")
        # per-condition comparison
        if condition_key and condition_key in adata.obs and adata.obs[condition_key].nunique() == 2:
            conds = sorted(adata.obs[condition_key].astype(str).unique())
            per = {}
            for c in conds:
                bc = b[b.obs[condition_key].astype(str) == c].copy()
                if bc.obs[key].value_counts().min() < 5 or bc.obs[key].nunique() < 2:
                    continue
                li.mt.rank_aggregate(bc, groupby=key, resource_name=resource, expr_prop=0.1, use_raw=False, n_perms=200, seed=0, verbose=False)
                d = bc.uns["liana_res"]; d = d[d.specificity_rank < 0.05]
                per[c] = set(zip(d.source, d.target, d.ligand_complex, d.receptor_complex))
            if len(per) == 2:
                a_, b_ = conds
                only_b = per[b_] - per[a_]; only_a = per[a_] - per[b_]
                rows = [[b_, *t] for t in sorted(only_b)[:12]] + [[a_, *t] for t in sorted(only_a)[:12]]
                job.table("ccc_condition", f"Interactions specific to one {condition_key}", ["Only in", "Sender", "Receiver", "Ligand", "Receptor"], rows,
                          note=f"Specific (specificity rank < 0.05) in one condition and not the other. {b_}: {len(only_b)} unique, {a_}: {len(only_a)} unique, shared: {len(per[a_] & per[b_])}.")
        job.method("Cell–cell communication", f"LIANA rank_aggregate ({resource} resource, expr_prop 0.1, 500 permutations) on log-normalized counts grouped by {key}.")
        return True
    except Exception:  # noqa: BLE001
        job.flag("info", f"Cell–cell communication skipped ({log_exc('liana')}).")
        return False


# ---------------------------------------------------------------- clustering quality
def silhouette(adata, keys: list[str], use_rep: str = "X_pca", n: int = 5000) -> dict[str, float]:
    try:
        from sklearn.metrics import silhouette_score
        rep = adata.obsm[use_rep] if use_rep in adata.obsm else adata.obsm["X_pca"]
        rng = np.random.default_rng(0)
        idx = rng.choice(adata.n_obs, min(n, adata.n_obs), replace=False)
        out = {}
        for k in keys:
            lab = adata.obs[k].astype(str).values[idx]
            out[k] = float(silhouette_score(rep[idx], lab)) if len(set(lab)) > 1 else float("nan")
        return out
    except Exception:  # noqa: BLE001
        log_exc("silhouette"); return {}
