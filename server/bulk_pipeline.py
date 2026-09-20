"""Bulk RNA-seq pipeline: PyDESeq2 (Python port of DESeq2) for normalization, dispersion
estimation, Wald tests and apeGLM-style LFC shrinkage; VST for PCA / heatmaps; optional
Enrichr over-representation via GSEApy."""
from __future__ import annotations

import re
import warnings
from importlib.metadata import version
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import DOWN, NS, PALETTE, UP, Job, UserFacingError, commas, fmt_n, log_exc, pct, style
from .io_utils import collapse_technical_replicates, infer_groups, infer_sex, load_bulk, map_gene_ids

warnings.filterwarnings("ignore")
REF_RE = re.compile(r"ctrl|control|untreat|vehicle|^wt$|wild|mock|dmso|naive|baseline|normal|uninf|non.?infect|sham|pbs|^no|healthy", re.I)


def _mad(x):
    return np.median(np.abs(x - np.median(x))) * 1.4826


def run(job: Job, files: list[Path], params: dict):
    import seaborn as sns
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    style()
    params = params or {}
    job.step(3, "Reading counts")
    counts, meta, notes = load_bulk(files)
    counts = map_gene_ids(counts, notes)
    cr = params.get("collapse_replicates", "auto")
    counts, meta, notes, collapsed = collapse_technical_replicates(
        counts, meta, notes, force=None if cr in (None, "auto") else bool(cr))
    counts.to_csv(job.root / "counts_used.csv")
    S, G0 = counts.shape[1], counts.shape[0]
    if meta is None:
        meta = pd.DataFrame({"condition": infer_groups(list(counts.columns))}, index=counts.columns)
        notes.append("No metadata file — groups were inferred from sample names. Upload a metadata CSV (sample, condition, batch…) for full control.")
    meta = meta.copy()
    meta.index = meta.index.astype(str)
    numeric_cols = []
    meta = meta.replace({"": np.nan, "NA": np.nan, "nan": np.nan, "None": np.nan})
    meta = meta.dropna(axis=1, how="all")
    for c in meta.columns:
        col = pd.to_numeric(meta[c], errors="coerce")
        if col.notna().all() and col.nunique() > 5:
            meta[c] = col.astype(float)          # continuous covariate (age, passage, RIN…)
            numeric_cols.append(c)
        else:
            meta[c] = meta[c].where(meta[c].isna(), meta[c].astype(str))
    if params.get("groups"):
        g = [x.strip() for x in params["groups"].split(",")]
        if len(g) != S:
            raise UserFacingError(f"Enter exactly {S} group labels (you entered {len(g)}).")
        meta["condition"] = g

    has_sex_col = any(c.lower() in ("sex", "gender") for c in meta.columns)
    sex_inferred = False
    if params.get("infer_sex", True) and not has_sex_col:
        sx = infer_sex(counts)
        if sx is not None:
            sex_inferred = True
            meta["sex_inferred"] = sx.reindex(meta.index).values
            notes.append("Sex was inferred from XIST vs Y-linked gene expression (column sex_inferred): "
                         + ", ".join(f"{k} {v}" for k, v in sx.value_counts().items()) + ".")
    meta.to_csv(job.root / "metadata_used.csv")

    # choose design factor: prefer complete columns; samples lacking a value for the factor are dropped
    complete = [c for c in meta.columns if meta[c].notna().all()]
    cand = [c for c in meta.columns if c in numeric_cols or 1 < meta[c].nunique() < S]
    cand = sorted(cand, key=lambda c: (c not in complete, ))
    if not cand:
        raise UserFacingError("No metadata column splits the samples into groups with replicates. Add a 'condition' column with at least two levels.")
    if not [c for c in cand if c not in numeric_cols]:
        raise UserFacingError("The design factor must be a categorical column (e.g. condition), but only numeric columns were found.")
    cat_cand = [c for c in cand if c not in numeric_cols]
    pref = [c for c in cat_cand if not re.search(r"sex|gender|batch|donor|source|line|replicate|lane|protocol|provider|gsm|title", c, re.I) and meta[c].nunique() <= 8] or cat_cand
    def _score(c):
        return (bool(re.search(r"condition|treat|group|genotype|infection|disease|status", c, re.I)),
                int(meta[c].notna().sum()), -int(meta[c].nunique()))
    factor = params.get("factor") if params.get("factor") in cat_cand else max(pref, key=_score)
    covars = [c for c in (params.get("covariates") or []) if c in cand and c != factor]
    missing = meta[factor].isna()
    if missing.any():
        notes.append(f"{int(missing.sum())} sample(s) have no value for {factor} and were left out of this comparison: {', '.join(meta.index[missing][:8])}.")
        meta = meta[~missing]
        counts = counts[meta.index]
        S = counts.shape[1]
    for c in covars:
        if meta[c].isna().any():
            meta[c] = meta[c].fillna("unknown")
    for c in meta.columns:
        if c not in numeric_cols:
            meta[c] = meta[c].astype(str)
    # columns that are (nearly) unique per sample cannot be covariates and are not informative for PC association
    cand = [c for c in cand if c in numeric_cols or 1 < meta[c].nunique() <= max(2, S // 2) or c == factor]
    levels = sorted(meta[factor].unique(), key=lambda l: (not REF_RE.search(l), l))
    ref = params.get("reference") if params.get("reference") in levels else levels[0]
    alts = [l for l in levels if l != ref]
    alt = params.get("alternative") if params.get("alternative") in alts else alts[0]
    alpha = float(params.get("alpha") or 0.05)
    lfc_thr = float(params.get("lfc") or 1.0)
    job.result["params"] = {"factor": factor, "reference": ref, "alternative": alt, "covariates": covars,
                            "alpha": alpha, "lfc": lfc_thr, "columns": cand, "levels": {c: sorted(meta[c].unique()) for c in cand}}

    colors = {l: PALETTE[i % 20] for i, l in enumerate([ref] + alts)}
    scol = meta[factor].map(colors)

    # ---------------------------------------------------------------- filtering + DESeq2
    min_n = int(meta[factor].value_counts().min())
    keep = (counts >= 10).sum(axis=1) >= min_n
    cf = counts[keep]
    job.step(12, f"Running DESeq2 on {len(cf):,} genes × {S} samples")
    design = "~" + " + ".join(covars + [factor])
    # columns that describe the same split of the samples as the design factor: their effect can never be
    # separated from it, so they must not be offered (or advised) as covariates
    confounded = {c for c in cand if c != factor and c not in numeric_cols
                  and (meta.groupby(c, observed=True)[factor].nunique().max() <= 1
                       or meta.groupby(factor, observed=True)[c].nunique().max() <= 1)}
    inference = DefaultInference(n_cpus=4)
    # make the reference level the baseline of the design matrix (formulaic follows Categorical order), so the
    # contrast coefficient exists for LFC shrinkage regardless of alphabetical order
    meta[factor] = pd.Categorical(meta[factor], categories=[ref] + alts)
    # a covariate that never varies inside a level of the factor is indistinguishable from the factor itself;
    # DESeq2 would fail deep inside the GLM with "Singular matrix", which tells the user nothing
    dc = covars + [factor]
    X = pd.get_dummies(meta[dc].astype({c: str for c in dc if c not in numeric_cols}), drop_first=True, dtype=float)
    X.insert(0, "_intercept", 1.0)
    if np.linalg.matrix_rank(X.values) < X.shape[1]:
        bad = [c for c in covars if c in confounded]
        if bad:
            ex = bad[0]
            tab = pd.crosstab(meta[ex], meta[factor])
            shown = "; ".join(f"{i} → {', '.join(tab.columns[tab.loc[i] > 0])}" for i in tab.index[:4])
            raise UserFacingError(
                f"The design ~{' + '.join(dc)} cannot be fitted: <b>{ex}</b> and <b>{factor}</b> describe the same split of "
                f"your samples ({shown}), so no test can tell a {ex} effect apart from a {factor} effect. "
                f"Remove {ex} from “Adjust for”, or compare {factor} within a single {ex} group.")
        raise UserFacingError(
            f"The design ~{' + '.join(dc)} cannot be fitted: the covariates you chose overlap too much with {factor} (or "
            "with each other) for the model to separate them. Adjust for fewer variables — each one must vary within "
            f"every {factor} group.")
    dds = DeseqDataSet(counts=cf.T, metadata=meta, design=design, refit_cooks=True, inference=inference, quiet=True)
    dds.deseq2()
    job.step(45, "Variance-stabilizing transform")
    dds.vst(use_design=False)
    vst = pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names).T
    norm = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names).T
    size_factors = pd.Series(np.asarray(dds.obs["size_factors"] if "size_factors" in dds.obs else dds.obsm["size_factors"]), index=dds.obs_names)

    results = {}
    for a in alts:
        job.step(55, f"Wald test: {a} vs {ref}")
        ds = DeseqStats(dds, contrast=[factor, a, ref], alpha=alpha, inference=inference, quiet=True)
        ds.summary()
        res = ds.results_df.copy()
        res["log2FC_shrunk"] = np.nan
        try:
            coeff = next(c for c in dds.varm["LFC"].columns if factor in c and a in c)
            ds.lfc_shrink(coeff=coeff)
            res["log2FC_shrunk"] = ds.results_df["log2FoldChange"]
        except Exception:  # noqa: BLE001
            log_exc("lfc_shrink")
            res["log2FC_shrunk"] = res["log2FoldChange"]
        results[a] = res
    res = results[alt]
    res.index.name = "gene"
    sig = res[(res.padj < alpha) & (res.log2FoldChange.abs() >= lfc_thr)]
    n_fdr_only = int(((res.padj < alpha) & (res.log2FoldChange.abs() < lfc_thr)).sum())
    up, down = sig[sig.log2FoldChange > 0], sig[sig.log2FoldChange < 0]
    out = res.join(norm.add_prefix("norm_"))
    out.sort_values("pvalue").to_csv(job.root / "deseq2_results.csv")
    norm.to_csv(job.root / "normalized_counts.csv")
    vst.to_csv(job.root / "vst_counts.csv")
    job.download(f"DESeq2 results {alt} vs {ref} (.csv)", "deseq2_results.csv")
    job.download("Normalized counts (.csv)", "normalized_counts.csv")
    job.download("VST expression (.csv)", "vst_counts.csv")
    job.download("Sample sheet used (.csv)", "metadata_used.csv")
    job.download("Counts used (.csv)", "counts_used.csv")

    # ---------------------------------------------------------------- QC stats
    lib = counts.sum()
    det = (counts > 0).sum()
    corr = vst.corr(method="pearson")
    mc = (corr.sum() - 1) / (S - 1)
    outliers = mc[mc < np.median(mc) - max(0.01, 4 * _mad(mc.values))].index.tolist()
    topv = vst.loc[vst.var(axis=1).sort_values(ascending=False).index[:500]]
    X = (topv.T - topv.T.mean()).values
    U, sv, Vt = np.linalg.svd(X, full_matrices=False)
    pcs = U * sv
    pvar = sv ** 2 / (sv ** 2).sum()
    pc_df = pd.DataFrame(pcs[:, :min(6, S)], index=topv.columns, columns=[f"PC{i+1}" for i in range(min(6, S))])
    # R² of metadata vs PCs
    r2 = pd.DataFrame(index=cand, columns=pc_df.columns, dtype=float)
    for c in cand:
        for p in pc_df.columns:
            y = pc_df[p]
            if c in numeric_cols:
                r2.loc[c, p] = float(np.corrcoef(y.values, meta[c].values)[0, 1] ** 2)
            else:
                gm = y.groupby(meta[c]).transform("mean")
                r2.loc[c, p] = 1 - ((y - gm) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    fac_pc = r2.loc[factor, ["PC1", "PC2"]].astype(float)
    strong_other = [(c, p, r2.loc[c, p]) for c in cand if c != factor and c not in covars for p in ["PC1", "PC2"] if r2.loc[c, p] > 0.5 and pvar[int(p[2]) - 1] > 0.1]

    # ---------------------------------------------------------------- overview
    job.tile(str(S), "samples")
    job.tile(f"{alt} vs {ref}", f"contrast ({factor})")
    job.tile(fmt_n(np.median(lib)), "median library size")
    job.tile(commas(len(cf)), "genes tested")
    job.tile(commas(len(up)), f"up (padj<{alpha}, |LFC|≥{lfc_thr:g})")
    job.tile(commas(len(down)), "down")
    for n in notes:
        job.flag("warn", n[5:]) if n.startswith("WARN:") else job.flag("info", n)
    sexcol = next((c for c in meta.columns if c.lower() in ("sex", "gender", "sex_inferred")), None)
    if sexcol and sexcol not in covars and sexcol != factor:
        ct = pd.crosstab(meta[factor], meta[sexcol])
        if ct.shape[1] > 1 and (ct.div(ct.sum(1), axis=0).max(1) > 0.85).any():
            job.flag("warn", f"Sex is unevenly distributed across {factor} groups ({'; '.join(f'{g}: ' + ', '.join(f'{v} {sx}' for sx, v in r.items()) for g, r in ct.iterrows())}). Add <b>{sexcol}</b> as a covariate to avoid sex-linked genes appearing as hits.")
    ratio = lib.max() / max(lib.min(), 1)
    job.flag("warn" if ratio > 3 else "ok", f"Library sizes range {fmt_n(lib.min())}–{fmt_n(lib.max())} ({ratio:.1f}-fold). DESeq2 size factors ({size_factors.min():.2f}–{size_factors.max():.2f}) correct for this"
             + ("; very shallow samples still give noisier estimates." if ratio > 3 else "."))
    job.flag("warn" if outliers else "ok", (f"<b>{', '.join(outliers)}</b> correlate(s) poorly with the other samples — a possible outlier, failed library or sample swap." if outliers
                                            else f"No outlier samples: mean correlation {mc.min():.3f}–{mc.max():.3f} on VST expression."))
    best_pc = fac_pc.idxmax()
    job.flag("ok" if fac_pc.max() > 0.5 else "warn",
             f"<b>{factor}</b> explains {100*fac_pc.max():.0f}% of the variance along {best_pc} ({100*pvar[int(best_pc[2])-1]:.0f}% of total)"
             + (" — the groups separate clearly." if fac_pc.max() > 0.5 else " — groups overlap; the effect is small relative to other variation."))
    for c, p, v in strong_other[:2]:
        # only suggest it if the design would still be fittable — a variable that tracks the factor cannot be adjusted for
        if c in confounded:
            job.flag("warn", f"<b>{c}</b> explains {100*v:.0f}% of {p}, but it describes the same split of the samples as "
                             f"{factor}, so it cannot be adjusted for — any {factor} difference you find may equally be a "
                             f"{c} difference. This is a limit of the experiment's design, not of the analysis.")
        else:
            job.flag("warn", f"<b>{c}</b> explains {100*v:.0f}% of {p}. Consider adding it as a covariate in the design (~ {c} + {factor}).")
    nA, nB = (meta[factor] == ref).sum(), (meta[factor] == alt).sum()
    if min(nA, nB) < 3:
        job.flag("warn", f"Only {min(nA, nB)} replicate(s) in a group — dispersion estimates are unreliable below 3; treat the gene list as exploratory.")
    top_up = up.sort_values("padj").head(6).index.tolist(); top_dn = down.sort_values("padj").head(6).index.tolist()
    job.flag("ok", f"DESeq2 ({design}): <b>{len(up):,}</b> genes up and <b>{len(down):,}</b> down in {alt} vs {ref} (padj < {alpha}, |log2FC| ≥ {lfc_thr:g}). "
             + (f"Top up: {', '.join(top_up)}. " if top_up else "") + (f"Top down: {', '.join(top_dn)}. " if top_dn else "")
             + (f"A further {n_fdr_only:,} genes are significant but change less than {2**lfc_thr:g}-fold." if n_fdr_only else ""))
    if len(alts) > 1:
        job.flag("info", "Other contrasts vs " + ref + ": " + "; ".join(f"{a}: {((r.padj < alpha) & (r.log2FoldChange.abs() >= lfc_thr)).sum():,} DE genes" for a, r in results.items() if a != alt) + ". Switch the contrast in settings.")

    # ---------------------------------------------------------------- QC figures
    job.section("qc", "Step 1", "Sample quality control", "Confirm every library was sequenced adequately and resembles its replicates before trusting any comparison.")
    order = meta.sort_values(factor).index
    fig, axs = plt.subplots(1, 2, figsize=(10, max(3, 0.24 * S + 1)))
    axs[0].barh(range(S), lib[order] / 1e6, color=scol[order]); axs[0].set_yticks(range(S)); axs[0].set_yticklabels(order, fontsize=8); axs[0].invert_yaxis(); axs[0].set_xlabel("library size (million reads)")
    axs[1].barh(range(S), size_factors[order], color=scol[order]); axs[1].set_yticks([]); axs[1].invert_yaxis(); axs[1].set_xlabel("DESeq2 size factor"); axs[1].axvline(1, color="#888", lw=0.8, ls="--")
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[l]) for l in [ref] + alts]
    axs[1].legend(handles, [ref] + alts, title=factor, loc="lower right", fontsize=8)
    fig.tight_layout()
    job.figure("libsize", "Library size & size factors", fig, wide=True,
               how="DESeq2's median-of-ratios size factors should roughly track library size. A size factor far from its library size means the sample's composition differs (e.g. a few genes soaking up reads).",
               yours=f"Median library <b>{fmt_n(np.median(lib))}</b> reads; {commas(np.median(det))} genes detected per sample.")

    fig, ax = plt.subplots(figsize=(5, max(3, 0.24 * S + 1)))
    ln = np.log2(norm[order] + 1)
    bp = ax.boxplot([ln[c].values for c in order], vert=False, showfliers=False, patch_artist=True, widths=0.6)
    for patch, s_ in zip(bp["boxes"], order):
        patch.set_facecolor(colors[meta.loc[s_, factor]]); patch.set_alpha(0.6)
    ax.set_yticks(range(1, S + 1)); ax.set_yticklabels(order, fontsize=8); ax.invert_yaxis(); ax.set_xlabel("log2(normalized count + 1)")
    meds = ln.median()
    job.figure("dist", "Normalized count distributions", fig,
               how="After normalization the boxes should line up. A shifted or squashed box suggests degraded RNA, contamination, or a failed library.",
               yours=f"Sample medians span <b>{meds.max()-meds.min():.2f}</b> log2 units" + (" — well aligned." if meds.max() - meds.min() < 0.5 else " — larger than usual; inspect that sample."))

    cg = sns.clustermap(corr, cmap="viridis", row_colors=scol, col_colors=scol, figsize=(min(10, 3 + S * 0.32),) * 2,
                        xticklabels=True, yticklabels=True, dendrogram_ratio=0.12, cbar_pos=(0.02, 0.85, 0.02, 0.12))
    cg.ax_heatmap.tick_params(labelsize=7)
    job.figure("corr", "Sample-to-sample correlation (VST, clustered)", cg.fig, wide=S > 10,
               how="Hierarchical clustering on Pearson correlation of variance-stabilized expression. Replicates should cluster with each other; the colour bar shows each sample's group.",
               yours=f"Lowest pairwise r = <b>{corr.values.min():.3f}</b>." + (f" Outlier(s): {', '.join(outliers)}." if outliers else ""))

    job.section("pca", "Step 2", "Sample structure", "PCA on the 500 most variable genes after DESeq2's variance-stabilizing transform (the plotPCA recipe).")
    others = [c for c in cand if c != factor]
    shape_by = others[0] if others else None
    markers = ["o", "s", "^", "D", "v", "P", "X"]
    fig, ax = plt.subplots(figsize=(5.4, 4.4))
    for s_ in pc_df.index:
        mkr = markers[sorted(meta[shape_by].unique()).index(meta.loc[s_, shape_by]) % 7] if shape_by else "o"
        ax.scatter(pc_df.loc[s_, "PC1"], pc_df.loc[s_, "PC2"], color=colors[meta.loc[s_, factor]], marker=mkr, s=70, edgecolor="white", linewidth=0.8)
    if S <= 24:
        for s_ in pc_df.index:
            ax.annotate(s_, (pc_df.loc[s_, "PC1"], pc_df.loc[s_, "PC2"]), fontsize=6.5, xytext=(4, 3), textcoords="offset points", color="#444")
    ax.set_xlabel(f"PC1 ({100*pvar[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({100*pvar[1]:.1f}%)")
    h = [plt.Line2D([], [], marker="o", ls="", color=colors[l], label=l) for l in [ref] + alts]
    if shape_by:
        h += [plt.Line2D([], [], marker=markers[i % 7], ls="", color="#666", label=f"{shape_by}={v}") for i, v in enumerate(sorted(meta[shape_by].unique()))]
    ax.legend(handles=h, fontsize=7, bbox_to_anchor=(1.01, 1), loc="upper left")
    job.figure("pca", "PCA of samples", fig,
               how="Distance ≈ overall transcriptome difference. Replicates should sit together and groups should separate along PC1 or PC2." + (f" Marker shape shows {shape_by}." if shape_by else ""),
               yours=f"PC1 {100*pvar[0]:.0f}% and PC2 {100*pvar[1]:.0f}% of variance; <b>{factor}</b> accounts for {100*fac_pc['PC1']:.0f}% of PC1 and {100*fac_pc['PC2']:.0f}% of PC2.")
    if covars:
        from . import bulk_extra as _bx
        _bx.adjusted_pca(job, vst, meta, factor, covars, colors)

    fig, ax = plt.subplots(figsize=(4.8, max(2.2, 0.45 * len(cand) + 1.2)))
    sns.heatmap(r2.astype(float), annot=True, fmt=".2f", cmap="rocket_r", vmin=0, vmax=1, ax=ax, cbar_kws={"label": "R²", "shrink": 0.7}, annot_kws={"fontsize": 8})
    ax.set_xticklabels([f"{p}\n{100*pvar[i]:.0f}%" for i, p in enumerate(r2.columns)], fontsize=8, rotation=0)
    job.figure("pc_assoc", "Which metadata drives each PC?", fig,
               how="R² from a one-way ANOVA of each PC on each metadata column. It reveals hidden batch effects: any strong, unmodelled factor on an early PC belongs in the design formula.",
               yours=("; ".join(f"<b>{c}</b> drives {p} (R² {v:.2f})" for c, p, v in strong_other[:3])
                      + (" — consider it as a covariate." if any(c not in confounded for c, _, _ in strong_other[:3])
                         else " — already in the design or inseparable from " + factor + ".")) if strong_other else f"Only <b>{factor}</b> strongly associates with the leading PCs.")

    # ---------------------------------------------------------------- DE figures
    job.step(70, "Plotting differential expression")
    job.section("de", "Step 3", f"Differential expression: {alt} vs {ref}",
                f"PyDESeq2 negative-binomial GLM, design <code>{design}</code>, Wald test with Benjamini–Hochberg FDR and independent filtering; "
                f"log2 fold changes shrunk with an apeGLM-style prior for ranking and MA plots. Positive = higher in {alt}.")
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    v = dds.var
    ax.scatter(v["_normed_means"] if "_normed_means" in v else norm.mean(1), v["genewise_dispersions"], s=3, c="#111", alpha=0.35, label="gene-est", rasterized=True)
    ax.scatter(v["_normed_means"] if "_normed_means" in v else norm.mean(1), v["dispersions"], s=3, c="#3a8fd1", alpha=0.5, label="final", rasterized=True)
    xm = (v["_normed_means"] if "_normed_means" in v else norm.mean(1)).values
    o = np.argsort(xm)
    ax.plot(xm[o], v["fitted_dispersions"].values[o], color="#c9304f", lw=1.5, label="fitted trend")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("mean of normalized counts"); ax.set_ylabel("dispersion"); ax.legend(markerscale=3, fontsize=8)
    job.figure("dispersion", "Dispersion estimates", fig,
               how="DESeq2 shrinks each gene's dispersion (black) toward the fitted trend (red) to get final values (blue). Dispersion should fall as expression rises; a cloud far above the trend signals outliers or unmodelled variation.",
               yours=f"Typical dispersion at high expression ≈ <b>{np.median(v['dispersions'].values[xm > np.percentile(xm, 75)]):.3f}</b> (≈{100*np.sqrt(np.median(v['dispersions'].values[xm > np.percentile(xm, 75)])):.0f}% biological CV)" + (" — low, typical of inbred models or cell lines." if np.median(v['dispersions'].values[xm > np.percentile(xm, 75)]) < 0.05 else "."))

    r = res.dropna(subset=["pvalue"]).copy()
    r["nlp"] = -np.log10(r["pvalue"].clip(lower=1e-300))
    r["sig"] = np.where((r.padj < alpha) & (r.log2FoldChange >= lfc_thr), "up", np.where((r.padj < alpha) & (r.log2FoldChange <= -lfc_thr), "down", "ns"))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    r.loc[(r.sig == "ns") & (r.padj < alpha), "sig"] = "fdr"
    cmap = {"up": UP, "down": DOWN, "ns": NS, "fdr": "#7d8a93"}
    labels = {"up": f"up ({(r.sig == 'up').sum():,})", "down": f"down ({(r.sig == 'down').sum():,})", "ns": "not significant", "fdr": f"padj < {alpha}, small fold change ({(r.sig == 'fdr').sum():,})"}
    for k in ["ns", "fdr", "down", "up"]:
        d = r[r.sig == k]
        ax.scatter(d.log2FoldChange, d.nlp, s=7 if k in ("up", "down") else 4, c=cmap[k], alpha=0.85 if k != "ns" else 0.5, linewidths=0, label=labels[k], rasterized=True)
    ax.axvline(lfc_thr, ls="--", lw=0.8, color="#777"); ax.axvline(-lfc_thr, ls="--", lw=0.8, color="#777")
    thr = r.loc[r.padj < alpha, "pvalue"].max()
    if pd.notna(thr):
        ax.axhline(-np.log10(thr), ls="--", lw=0.8, color="#777")
    lab = pd.concat([r[r.sig == "up"].nsmallest(8, "pvalue"), r[r.sig == "down"].nsmallest(8, "pvalue"), r.nsmallest(6, "pvalue")])
    lab = lab[~lab.index.duplicated()]
    try:
        from adjustText import adjust_text
        texts = [ax.text(x, y, g, fontsize=7.5) for g, x, y in zip(lab.index, lab.log2FoldChange, lab.nlp)]
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#999", lw=0.5))
    except Exception:
        for g, x, y in zip(lab.index, lab.log2FoldChange, lab.nlp):
            ax.annotate(g, (x, y), fontsize=7.5, xytext=(3, 2), textcoords="offset points")
    ax.set_xlabel(f"log2 fold change ({alt} / {ref})"); ax.set_ylabel("−log10 p-value"); ax.legend(fontsize=8, loc="upper left")
    lead = sig.reindex(sig.log2FoldChange.abs().sort_values(ascending=False).index).head(1)
    job.figure("volcano", "Volcano plot", fig, wide=True,
               how=f"Right = higher in {alt}; higher = more significant. Dashed lines mark |log2FC| = {lfc_thr:g} (a {2**lfc_thr:g}-fold change) and the p-value corresponding to FDR {alpha}. The most interesting genes sit in the upper corners.",
               yours=(f"Most significant: <b>{', '.join(r.nsmallest(6, 'pvalue').index)}</b>. Largest significant change: {lead.index[0]} ({2**abs(lead.log2FoldChange.iloc[0]):.1f}× {'up' if lead.log2FoldChange.iloc[0] > 0 else 'down'})." if len(lead) else
                      f"No genes pass padj < {alpha} and |log2FC| ≥ {lfc_thr:g}. Try a smaller fold-change cut-off or check replicate consistency."))

    fig, ax = plt.subplots(figsize=(5, 3.9))
    rm = res.dropna(subset=["baseMean"]).copy()
    rm = rm[rm.baseMean > 0]
    sgn = rm.index.isin(sig.index)
    ax.scatter(rm.baseMean[~sgn], rm.log2FC_shrunk[~sgn], s=3, c=NS, linewidths=0, rasterized=True)
    ax.scatter(rm.baseMean[sgn], rm.log2FC_shrunk[sgn], s=5, c=np.where(rm.log2FC_shrunk[sgn] > 0, UP, DOWN), linewidths=0, rasterized=True)
    ax.axhline(0, color="#555", lw=0.8); ax.set_xscale("log"); ax.set_xlabel("mean of normalized counts"); ax.set_ylabel("shrunken log2 fold change")
    job.figure("ma", "MA plot (shrunken LFC)", fig,
               how="Fold change against mean expression. Shrinkage pulls noisy fold changes of low-count genes toward zero, so remaining large effects on the left are trustworthy.",
               yours=f"{commas(len(sig))} significant genes; median baseMean of significant genes {fmt_n(sig.baseMean.median()) if len(sig) else '—'}.")

    fig, ax = plt.subplots(figsize=(4.8, 3.9))
    pv = res.pvalue.dropna()
    ax.hist(pv[res.baseMean.reindex(pv.index) > 1], bins=40, color="#2b7a78", edgecolor="white")
    ax.axhline((res.baseMean > 1).sum() / 40 * (1 - (pv < 0.05).mean()), ls="--", color="#888", lw=0.8)
    ax.set_xlabel("p-value"); ax.set_ylabel("genes")
    p05 = (pv < 0.05).mean()
    job.figure("pvals", "P-value histogram", fig,
               how="A well-behaved test gives a flat floor of true nulls plus a peak near zero. A hill in the middle or a spike at 1 indicates a mis-specified model (missing covariate, outliers).",
               yours=f"<b>{100*p05:.1f}%</b> of genes have p < 0.05 (5% expected by chance)" + (" — strong differential signal." if p05 > 0.15 else " — modest signal." if p05 > 0.07 else " — little evidence of differential expression."))

    top = sig.nsmallest(50, "padj").index if len(sig) >= 2 else res.nsmallest(50, "pvalue").index
    hm = vst.loc[top, order]
    cg = sns.clustermap(hm, z_score=0, cmap="RdBu_r", center=0, vmin=-2.5, vmax=2.5, col_cluster=False, col_colors=scol[order],
                        figsize=(min(11, 4 + S * 0.35), max(5, len(top) * 0.17 + 1.5)), yticklabels=True, xticklabels=True,
                        dendrogram_ratio=(0.08, 0.04), cbar_pos=(0.01, 0.9, 0.015, 0.08))
    cg.ax_heatmap.tick_params(labelsize=6.5)
    job.figure("heatmap", f"Top {len(top)} DE genes (row z-score of VST)", cg.fig, wide=True,
               how="Rows are genes scaled to their own mean, clustered by pattern; columns are samples ordered by group. Consistent colour blocks switching at the group boundary show a reproducible effect.",
               yours="Genes cluster into up- and down-regulated modules; look for single samples that break the pattern.")

    show = r.nsmallest(8, "pvalue").index
    fig, axs = plt.subplots(2, 4, figsize=(10, 4.6))
    for ax, g in zip(axs.ravel(), show):
        dfp = pd.DataFrame({"count": norm.loc[g], "group": meta[factor]})
        sns.stripplot(data=dfp, x="group", y="count", hue="group", palette=colors, order=[ref] + alts, ax=ax, size=5, legend=False, jitter=0.15)
        sns.boxplot(data=dfp, x="group", y="count", order=[ref] + alts, ax=ax, showfliers=False, boxprops=dict(facecolor="none"), width=0.5, linewidth=0.8)
        ax.set_yscale("log"); ax.set_title(g, fontsize=9); ax.set_xlabel(""); ax.set_ylabel("")
        ax.tick_params(axis="x", labelsize=7, rotation=30)
    axs[0, 0].set_ylabel("normalized count"); axs[1, 0].set_ylabel("normalized count")
    fig.tight_layout()
    job.figure("genecounts", "Top genes: per-sample counts", fig, wide=True,
               how="The equivalent of DESeq2's plotCounts. Always look at the raw points: a tiny p-value driven by one extreme sample is not a robust finding.",
               yours="")

    rows = [[g, f"{x.baseMean:,.0f}", f"{x.log2FoldChange:.2f}", f"{x.log2FC_shrunk:.2f}", f"{x.lfcSE:.2f}", f"{x.pvalue:.1e}", f"{x.padj:.1e}"]
            for g, x in res.dropna(subset=["padj"]).nsmallest(40, "padj").iterrows()]
    job.table("de_table", "Top 40 genes by adjusted p-value", ["Gene", "baseMean", "log2FC", "shrunk LFC", "lfcSE", "p-value", "padj"], rows,
              note="Full results for every gene in deseq2_results.csv.", csv="deseq2_results.csv")

    # ---------------------------------------------------------------- expert add-ons
    from . import bulk_extra
    mouse = np.mean([bool(re.match(r"^[A-Z][a-z0-9]+$", g)) for g in res.index[:500]]) > 0.5
    organism = "mouse" if mouse else "human"
    goi = params.get("genes")
    if goi:
        goi = goi if isinstance(goi, list) else [g for g in re.split(r"[,\s;]+", str(goi)) if g]
        bulk_extra.genes_of_interest(job, goi, norm, meta, factor, res, ref, alt)
    if params.get("gsea", True):
        job.step(80, "GSEA (pre-ranked)")
        bulk_extra.gsea(job, res, organism, alt, ref)
    if params.get("activity", True):
        job.step(85, "Pathway and TF activity (decoupler)")
        bulk_extra.activity(job, res, vst, meta, factor, ref, alt, organism)

    # ---------------------------------------------------------------- enrichment
    job.step(88, "Pathway enrichment (Enrichr)")
    enriched = False
    if params.get("enrichment", True) and (len(up) >= 5 or len(down) >= 5):
        try:
            import gseapy as gp
            libs = ["MSigDB_Hallmark_2020", "GO_Biological_Process_2023", "KEGG_2021_Human" if not mouse else "KEGG_2019_Mouse"]
            frames = []
            for name, genes in (("up", up.index), ("down", down.index)):
                if len(genes) < 5:
                    continue
                gl = [g.upper() for g in genes] if mouse and "Hallmark" in libs[0] else list(genes)
                e = gp.enrichr(gene_list=gl, gene_sets=libs, organism="human", outdir=None, background=[g.upper() for g in cf.index] if mouse else list(cf.index), cutoff=1, verbose=False)
                d = e.results.copy(); d["direction"] = name
                frames.append(d)
            er = pd.concat(frames)
            er.to_csv(job.root / "enrichment.csv", index=False)
            job.download("Pathway enrichment (.csv)", "enrichment.csv")
            job.section("enrich", "Step 4", "Pathway enrichment",
                        "Over-representation analysis of up- and down-regulated genes against MSigDB Hallmark, GO Biological Process and KEGG (Enrichr via GSEApy), using tested genes as background.")
            fig, axs = plt.subplots(1, 2, figsize=(11, 5))
            for ax, dname, col in zip(axs, ["up", "down"], [UP, DOWN]):
                d = er[(er.direction == dname)].nsmallest(12, "Adjusted P-value")
                if d.empty:
                    ax.axis("off"); continue
                d = d.iloc[::-1]
                ax.barh(range(len(d)), -np.log10(d["Adjusted P-value"]), color=col, alpha=0.85)
                ax.set_yticks(range(len(d))); ax.set_yticklabels([re.sub(r"\s*\(GO:\d+\)", "", t)[:48] for t in d.Term], fontsize=7.5)
                ax.set_xlabel("−log10 adj. p"); ax.set_title(f"{dname}-regulated in {alt}")
                ax.axvline(-np.log10(0.05), ls="--", color="#777", lw=0.8)
            fig.tight_layout()
            tu = er[(er.direction == "up") & (er["Adjusted P-value"] < 0.05)].nsmallest(3, "Adjusted P-value").Term.tolist()
            td = er[(er.direction == "down") & (er["Adjusted P-value"] < 0.05)].nsmallest(3, "Adjusted P-value").Term.tolist()
            job.figure("enrichment", "Enriched pathways", fig, wide=True,
                       how="Pathways containing more of your DE genes than expected by chance. Enrichment suggests biology to investigate, not proof of pathway activity.",
                       yours=(f"Up: <b>{'; '.join(tu)}</b>. " if tu else "No significant pathway among up genes. ") + (f"Down: <b>{'; '.join(td)}</b>." if td else "No significant pathway among down genes."))
            if tu or td:
                job.flag("ok", "Enriched pathways — up: " + ("; ".join(tu[:2]) or "none") + " · down: " + ("; ".join(td[:2]) or "none") + ".")
            enriched = True
        except Exception as e:  # noqa: BLE001
            job.flag("info", f"Pathway enrichment skipped — Enrichr could not be reached ({type(e).__name__}). It needs an internet connection.")

    job.method("Input", f"{G0:,} genes × {S} samples of raw counts" + (" after summing technical replicates" if collapsed else "") + f"; metadata columns: {', '.join(meta.columns)}"
               + (f" (continuous: {', '.join(numeric_cols)})" if numeric_cols else "") + ".")
    job.method("Filtering", f"Genes kept with ≥10 counts in at least {min_n} samples (size of the smallest group): {len(cf):,} genes.")
    job.method("DESeq2", f"PyDESeq2 {version('pydeseq2')}: median-of-ratios size factors, dispersion trend + MAP shrinkage, negative-binomial GLM with design {design}, Cooks outlier refitting, Wald test, BH FDR with independent filtering (alpha={alpha}).")
    job.method("LFC shrinkage", "DeseqStats.lfc_shrink (apeGLM-style Cauchy prior) on the contrast coefficient.")
    job.method("Transform", "Variance-stabilizing transformation (vst, blind to design) for PCA, correlation and heatmaps; PCA on the 500 most variable genes.")
    if enriched:
        job.method("Enrichment", "GSEApy Enrichr over-representation (Hallmark 2020, GO BP 2023, KEGG), background = tested genes, BH-adjusted.")
    job.method("Caveats", "Significance depends on replicates and a correct design. Interaction terms, paired designs beyond additive covariates, and time-course models need a custom analysis.")
    # ---------------------------------------------------------------- parameters panel
    try:
        job.result["params"].update({
            "collapse_replicates": cr if cr in ("auto", None) else str(bool(cr)).lower(),
            "infer_sex": bool(params.get("infer_sex", True)), "genes": params.get("genes") or "",
            "gsea": bool(params.get("gsea", True)), "activity": bool(params.get("activity", True)),
            "enrichment": bool(params.get("enrichment", True))})
        job.result["params_panel"] = _panel(params, {
            "S": S, "G0": G0, "collapsed": collapsed, "collapse": cr if cr in ("auto", None) else str(bool(cr)).lower(),
            "infer_sex": params.get("infer_sex", True), "sex_inferred": sex_inferred, "has_sex_col": has_sex_col,
            "columns": cand, "levels": {c: sorted(meta[c].unique()) for c in cand},
            "factor": factor, "ref": ref, "alt": alt, "covars": covars, "design": design,
            "min_n": min_n, "n_genes": len(cf), "alpha": alpha, "lfc": lfc_thr, "confounded": sorted(confounded),
            "n_up": len(up), "n_down": len(down), "flags": job.result.get("flags") or []})
    except Exception:  # noqa: BLE001
        log_exc("params panel")

    job.result["versions"] = {p: version(p) for p in ["pydeseq2", "numpy", "pandas", "scipy", "seaborn", "gseapy"]}
    job.result["summary"] = {"up": len(up), "down": len(down)}


def _panel(params, ctx):
    """Every choice this run made, with the reason it was made (see server/params_spec.py)."""
    from .params_spec import F, G, skipped, src
    cols, lv = ctx["columns"], ctx["levels"]
    flags = ctx["flags"]
    factor, ref, alt = ctx["factor"], ctx["ref"], ctx["alt"]
    groups = [
        G("input", "Input and samples", [
            F("collapse_replicates", "Sum technical replicates", ctx["collapse"],
              ("Columns that were the same library sequenced more than once were summed — technical replicates are not "
               "independent observations and DESeq2 expects one column per biological sample."
               if ctx["collapsed"] else
               "No technical replicates were detected, so every column was treated as its own biological sample. "
               "Set this if several columns are the same library run on different lanes."),
              type="select", source=src(params, "collapse_replicates"), options=["auto", "true", "false"]),
            F("groups", "Sample group labels", params.get("groups") or "",
              (f"Taken from your metadata file ({ctx['S']} samples, columns: {', '.join(cols)})." if not params.get("groups")
               else "You typed these labels, overriding the metadata."),
              type="text", source=src(params, "groups"), placeholder="only needed without a metadata file"),
            F("infer_sex", "Infer sex from XIST / Y genes", bool(ctx["infer_sex"]),
              ("No sex column was in your metadata, so it was inferred from XIST versus Y-linked expression and added as "
               "sex_inferred — worth adjusting for, since it is a common hidden source of variance."
               if ctx["sex_inferred"] else
               "Your metadata already has a sex column, so nothing was inferred." if ctx["has_sex_col"] else
               "Could not be inferred from these genes."),
              type="bool", source=src(params, "infer_sex"))],
            note=f"{ctx['S']} samples × {ctx['G0']:,} genes."),

        G("design", "Experimental design", [
            F("factor", "Design factor", factor,
              (f"You chose '{factor}'." if src(params, "factor") == "user" else
               f"'{factor}' was picked automatically: it splits the samples into {len(lv.get(factor, []))} groups with "
               "replicates, has a value for every sample, and its name matches the usual condition-like columns. "
               "This is the variable being tested — if it is the wrong one, nothing below is meaningful."),
              type="select", source=src(params, "factor", auto=True), options=cols),
            F("reference", "Reference (baseline) level", ref,
              (f"You set '{ref}' as the baseline." if src(params, "reference") == "user" else
               f"'{ref}' was taken as the baseline because its name looks like a control. ") +
              f"It is the denominator: a positive log2 fold change means higher in {alt} than in {ref}. Getting this "
              "backwards flips the sign of every result.",
              type="select", source=src(params, "reference", auto=True), options=lv.get(factor, [])),
            F("alternative", "Compared against", alt,
              f"The level being tested against {ref}. Only two levels are compared at a time; re-run with a different one "
              "to test another contrast." + (f" Other levels present: {', '.join(x for x in lv.get(factor, []) if x not in (ref, alt))}."
                                             if len(lv.get(factor, [])) > 2 else ""),
              type="select", source=src(params, "alternative", auto=True), options=lv.get(factor, [])),
            F("covariates", "Adjust for", ctx["covars"],
              (f"The model is {ctx['design']} — {', '.join(ctx['covars'])} "
               f"{'is' if len(ctx['covars']) == 1 else 'are'} held constant, so a difference in "
               f"{'it' if len(ctx['covars']) == 1 else 'them'} is not read as an effect of {factor}." if ctx["covars"] else
               "Nothing is adjusted for: the model is ~" + factor + ". If your samples came in batches, or differ in sex, "
               "age or RIN, add that column here — an unadjusted nuisance variable both hides real effects and creates "
               "false ones.") +
              " Only add a variable that varies within each group." +
              (f" Not offered here: {', '.join(ctx['confounded'])} — each splits your samples exactly the way {factor} does, "
               "so its effect could never be told apart from the one you are testing."
               if ctx["confounded"] else
               " A variable that tracks the condition would remove the very effect you are testing."),
              type="multiselect", source=src(params, "covariates"),
              options=[c for c in cols if c != factor and c not in ctx["confounded"]])],
            note=f"Design: {ctx['design']}"),

        G("filter", "Filtering and significance", [
            F("_filter", "Gene filter", f"≥ 10 counts in ≥ {ctx['min_n']} samples → {ctx['n_genes']:,} genes tested",
              f"{ctx['min_n']} is the size of your smallest group, so a gene expressed only in one group still survives. "
              "Genes below this carry no information and would only cost statistical power through multiple testing. Fixed.",
              type="fixed", source="auto"),
            F("alpha", "FDR threshold (padj)", ctx["alpha"],
              "The false discovery rate you accept among the genes you call significant: at 0.05, about 5% of them are "
              "expected to be false. This is also handed to DESeq2 for independent filtering, so changing it slightly "
              "changes which genes are testable, not just which are labelled.",
              type="float", source=src(params, "alpha"), step=0.01, min=0.001, max=0.5),
            F("lfc", "Minimum |log2 fold change|", ctx["lfc"],
              f"A gene must also change by at least {2 ** float(ctx['lfc']):.3g}-fold (|log2FC| ≥ {ctx['lfc']}) to be called. "
              "With enough replicates, tiny but consistent changes become statistically significant without being "
              "biologically interesting; this filter is about effect size, not confidence. It is applied after testing.",
              type="float", source=src(params, "lfc"), step=0.25, min=0, max=5),
            F("_shrink", "Log fold-change shrinkage", "apeGLM-style, on the tested contrast",
              "Low-count genes produce wild fold changes from very little evidence. Shrinkage pulls those toward zero while "
              "leaving well-measured genes alone, which is what makes the volcano and the ranking trustworthy. Fixed.",
              type="fixed")],
            note=f"{ctx['n_up']} up, {ctx['n_down']} down at padj < {ctx['alpha']} and |log2FC| ≥ {ctx['lfc']}."),

        G("addons", "Optional analyses", [
            F("genes", "Genes of interest", params.get("genes") or "",
              "Your own gene list gets a per-sample expression panel and its own statistics table, whether or not the genes "
              "passed the thresholds above. Comma-separated symbols.",
              type="text", source=src(params, "genes"), placeholder="e.g. TNF, IL6, CXCL10"),
            F("gsea", "GSEA (pre-ranked)", bool(params.get("gsea", True)),
              "Ranks all tested genes by the Wald statistic and asks which gene sets concentrate at either end (Hallmark, "
              "GO BP, KEGG, Reactome). Unlike the enrichment below it uses no significance cut-off, so coordinated small "
              "changes across a pathway are still detected." + skipped(flags, "GSEA"),
              type="bool", source=src(params, "gsea")),
            F("activity", "Pathway and TF activity (decoupler)", bool(params.get("activity", True)),
              "PROGENy scores 14 signalling pathways from their downstream response genes, and CollecTRI infers "
              "transcription-factor activity from the behaviour of each factor's targets — both read the consequences of "
              "activity rather than the expression of the pathway members, which is usually a poor proxy."
              + skipped(flags, "activity"), type="bool", source=src(params, "activity")),
            F("enrichment", "Over-representation (Enrichr)", bool(params.get("enrichment", True)),
              ("Asks whether your significant up- and down-gene lists contain more members of a pathway than chance, with "
               "the tested genes as background. Needs internet." + skipped(flags, "enrichment")
               if ctx["n_up"] >= 5 or ctx["n_down"] >= 5 else
               "Skipped: fewer than 5 significant genes in either direction, which is too few to over-represent anything."),
              type="bool", source=src(params, "enrichment"))],
            note="Each is skipped with a note rather than failing when its reference data cannot be reached."),
    ]
    return groups
