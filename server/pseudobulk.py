"""Sample-level statistics for single-cell data: composition tests and pseudobulk DESeq2.
Cells are never treated as replicates — the sample (donor, plate, library) is."""
from __future__ import annotations

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from .common import DOWN, NS, PALETTE, UP, Job, commas, log_exc, pct

COND_RE = re.compile(r"condition|genotype|treat|group|disease|status|infection|stim|dose|timepoint|age", re.I)


def detect_keys(obs: pd.DataFrame, sample_key: str | None, condition_key: str | None):
    """Pick a sample column and a condition column that is constant within samples."""
    cols = list(obs.columns)
    if not sample_key:
        for k in ("sample", "plate", "plate_id", "batch", "orig.ident", "donor", "library", "sample_id", "patient", "gsm"):
            if k in obs and 1 < obs[k].nunique() <= 200:
                sample_key = k
                break
    if not sample_key or sample_key not in obs:
        return None, None
    per_sample = obs.groupby(sample_key, observed=True)
    if condition_key and condition_key in obs:
        return sample_key, condition_key
    best = None
    for c in cols:
        if c == sample_key or obs[c].dtype.kind in "fc" and obs[c].nunique() > 10:
            continue
        nun = obs[c].nunique()
        if not (2 <= nun <= 6):
            continue
        if (per_sample[c].nunique() <= 1).all() and obs[c].value_counts().min() > 0:
            score = 2 if COND_RE.search(c) else 1
            if best is None or score > best[0]:
                best = (score, c)
    return sample_key, (best[1] if best else None)


def run(job: Job, adata, sample_key: str, condition_key: str, covariates: list[str] | None, reference: str | None,
        cluster_key: str = "cell_type_simple", min_cells: int = 20, min_samples: int = 3):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    obs = adata.obs
    smeta = obs.groupby(sample_key, observed=True).agg(
        **{condition_key: (condition_key, "first")}, n_cells=(sample_key, "size"),
        **{c: (c, "first") for c in (covariates or []) if c in obs and c != condition_key})
    smeta[condition_key] = smeta[condition_key].astype(str)
    levels = sorted(smeta[condition_key].unique(), key=lambda l: (not re.search(r"ctrl|control|wt|wild|untreat|vehicle|mock|healthy|normal|naive", l, re.I), l))
    ref = reference if reference in levels else levels[0]
    alt = next(l for l in levels if l != ref)
    covars = [c for c in (covariates or []) if c in smeta and smeta[c].nunique() > 1 and c != condition_key]
    n_per = smeta[condition_key].value_counts()
    job.section("pseudobulk", "Step 5", f"Sample-level comparison: {alt} vs {ref}",
                f"Cells are grouped by <b>{sample_key}</b> ({len(smeta)} samples: " + ", ".join(f"{v} {k}" for k, v in n_per.items()) +
                "). Composition is compared across samples, and counts are summed per sample for DESeq2 — the same "
                "negative-binomial test used for bulk RNA-seq, with samples as the replicates.")
    out = {"sample_key": sample_key, "condition_key": condition_key, "reference": ref, "alternative": alt, "covariates": covars}
    if n_per.min() < 2:
        job.flag("warn", f"Only {n_per.min()} sample(s) in one {condition_key} group — composition and pseudobulk tests were skipped (need ≥ 2, ideally ≥ 3 per group).")
        return out

    # ---------------------------------------------------------------- composition
    comp = pd.crosstab(obs[sample_key], obs[cluster_key], normalize="index") * 100
    comp = comp.loc[smeta.sort_values([condition_key]).index]
    cats = list(comp.columns)
    fig, ax = plt.subplots(figsize=(min(14, 4 + len(comp) * 0.3), 4.4))
    bottom = np.zeros(len(comp))
    for i, c in enumerate(cats):
        ax.bar(range(len(comp)), comp[c], bottom=bottom, color=PALETTE[i % 20], label=c, width=0.85)
        bottom += comp[c].values
    ax.set_xticks(range(len(comp)))
    ax.set_xticklabels([f"{s} · {smeta.loc[s, condition_key]}" for s in comp.index], rotation=90, fontsize=7)
    ax.set_ylabel("% of cells"); ax.set_ylim(0, 100)
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)
    fig.tight_layout()
    rows = []
    for c in cats:
        a = comp.loc[smeta[condition_key] == alt, c]; b = comp.loc[smeta[condition_key] == ref, c]
        try:
            p = float(mannwhitneyu(a, b).pvalue) if len(a) > 1 and len(b) > 1 else np.nan
        except Exception:  # noqa: BLE001
            p = np.nan
        rows.append([c, b.mean(), a.mean(), a.mean() - b.mean(), p])
    prop = pd.DataFrame(rows, columns=["cell_type", f"{ref}_pct", f"{alt}_pct", "diff", "p"])
    prop["padj"] = np.minimum(1, prop.p * prop.p.notna().sum())
    prop = prop.reindex(prop["diff"].abs().sort_values(ascending=False).index)
    prop.to_csv(job.root / "composition_test.csv", index=False)
    sig = prop[prop.padj < 0.05]
    top = prop.iloc[0]
    job.figure("composition_by_sample", f"Cell-type composition per {sample_key}", fig, wide=True,
               how="Each bar is one sample. Because samples are the replicates, a population that differs consistently across samples is a real shift; one that differs in a single sample is not.",
               yours=(f"<b>{top.cell_type}</b> is {top[f'{alt}_pct']:.1f}% of cells in {alt} vs {top[f'{ref}_pct']:.1f}% in {ref}"
                      + (f" (adjusted p = {top.padj:.1e})." if pd.notna(top.padj) else ".")
                      + (f" {len(sig)} population(s) differ at adjusted p &lt; 0.05." if len(sig) else " No population differs significantly after correction.")))
    job.table("composition", f"Composition test, {alt} vs {ref} (Mann–Whitney across samples, Bonferroni)",
              ["Cell type", f"% in {ref}", f"% in {alt}", "Difference", "p", "adj. p"],
              [[r.cell_type, f"{r[f'{ref}_pct']:.1f}", f"{r[f'{alt}_pct']:.1f}", f"{r['diff']:+.1f}", f"{r.p:.2e}" if pd.notna(r.p) else "—", f"{r.padj:.2e}" if pd.notna(r.padj) else "—"] for _, r in prop.iterrows()],
              note="Percentages are means across samples.", csv="composition_test.csv")
    if len(sig):
        job.flag("ok", f"Composition shift: <b>{sig.iloc[0].cell_type}</b> {sig.iloc[0][f'{ref}_pct']:.1f}% → {sig.iloc[0][f'{alt}_pct']:.1f}% ({ref} → {alt}), adjusted p = {sig.iloc[0].padj:.1e}.")
    job.download("Composition test (.csv)", "composition_test.csv")

    # ---------------------------------------------------------------- pseudobulk DE
    C = adata.layers["counts"]
    genes = adata.var_names
    inference = DefaultInference(n_cpus=4)

    def pseudobulk(mask):
        samples = [s for s in smeta.index if ((obs[sample_key] == s).values & mask).sum() >= min_cells]
        M = np.vstack([np.asarray(C[(obs[sample_key] == s).values & mask].sum(0)).ravel() for s in samples])
        return pd.DataFrame(np.round(M).astype(int), index=samples, columns=genes).T, smeta.loc[samples]

    def deseq(cnt, meta, tag):
        meta = meta.copy()
        for c in [condition_key] + covars:
            meta[c] = meta[c].astype(str)
        keep = (cnt >= 10).sum(1) >= max(2, int(meta[condition_key].value_counts().min()))
        cf = cnt[keep]
        design = "~" + " + ".join(covars + [condition_key])
        meta[condition_key] = pd.Categorical(meta[condition_key], categories=[ref] + [l for l in meta[condition_key].unique() if l != ref])
        dds = DeseqDataSet(counts=cf.T, metadata=meta, design=design, refit_cooks=True, inference=inference, quiet=True)
        dds.deseq2()
        ds = DeseqStats(dds, contrast=[condition_key, alt, ref], inference=inference, quiet=True)
        ds.summary()
        r = ds.results_df.copy()
        r["log2FC_shrunk"] = r["log2FoldChange"]
        try:
            coeff = next(c for c in dds.varm["LFC"].columns if condition_key in c and alt in c)
            ds.lfc_shrink(coeff=coeff)
            r["log2FC_shrunk"] = ds.results_df["log2FoldChange"]
        except Exception:  # noqa: BLE001
            log_exc("pseudobulk lfc_shrink")
        r = r.sort_values("pvalue"); r.index.name = "gene"
        r.to_csv(job.root / f"pseudobulk_{tag}.csv")
        s = r[(r.padj < 0.05) & (r.log2FoldChange.abs() >= 1)]
        return r, s, design, len(cf)

    def volcano(r, s, title):
        r = r.dropna(subset=["pvalue"]).copy(); r["nlp"] = -np.log10(r.pvalue.clip(lower=1e-300))
        sig = r.index.isin(s.index)
        fig, ax = plt.subplots(figsize=(7, 4.6))
        ax.scatter(r.log2FoldChange[~sig], r.nlp[~sig], s=3, c=NS, alpha=.5, linewidths=0, rasterized=True)
        ax.scatter(r.log2FoldChange[sig], r.nlp[sig], s=8, c=np.where(r.log2FoldChange[sig] > 0, UP, DOWN), linewidths=0, rasterized=True)
        for v in (-1, 1):
            ax.axvline(v, ls="--", lw=.8, color="#777")
        for g, x, y in zip(r[sig].head(14).index, r[sig].head(14).log2FoldChange, r[sig].head(14).nlp):
            ax.annotate(g, (x, y), fontsize=7.5, xytext=(3, 2), textcoords="offset points")
        ax.set_xlabel(f"log2 fold change ({alt} / {ref})"); ax.set_ylabel("−log10 p"); ax.set_title(title, fontsize=10)
        return fig

    job.step(90, "Pseudobulk DESeq2")
    cnt, meta = pseudobulk(np.ones(adata.n_obs, bool))
    r_all, s_all, design, ng = deseq(cnt, meta, "all_cells")
    up, down = s_all[s_all.log2FoldChange > 0], s_all[s_all.log2FoldChange < 0]
    job.figure("pseudobulk_volcano_all", f"Pseudobulk volcano, all cells: {alt} vs {ref}", volcano(r_all, s_all, f"All cells · {len(meta)} samples · {design}"), wide=True,
               how="Counts summed per sample, then DESeq2 exactly as for bulk. This avoids the inflated p-values that come from treating each cell as an independent replicate.",
               yours=f"<b>{commas(len(up))}</b> genes up and <b>{commas(len(down))}</b> down (padj &lt; 0.05, ≥2-fold) across {ng:,} tested genes."
                     + (f" Top up: {', '.join(up.head(6).index)}." if len(up) else "") + (f" Top down: {', '.join(down.head(6).index)}." if len(down) else ""))
    job.table("pseudobulk_all", f"Top genes, all cells ({alt} vs {ref})", ["Gene", "log2FC", "baseMean", "p", "padj"],
              [[g, f"{x.log2FoldChange:.2f}", f"{x.baseMean:,.0f}", f"{x.pvalue:.1e}", f"{x.padj:.1e}"] for g, x in r_all.dropna(subset=["padj"]).head(25).iterrows()],
              note="Full results in pseudobulk_all_cells.csv.", csv="pseudobulk_all_cells.csv")
    job.download(f"Pseudobulk DESeq2, all cells (.csv)", "pseudobulk_all_cells.csv")
    job.flag("ok", f"Pseudobulk DESeq2 ({design}, {len(meta)} samples): <b>{commas(len(up))}</b> up, <b>{commas(len(down))}</b> down in {alt} vs {ref}."
             + (f" Top up: {', '.join(up.head(5).index)}." if len(up) else ""))
    out["all_cells"] = {"up": len(up), "down": len(down), "design": design, "n_samples": len(meta)}
    out["ref"], out["alt"], out["condition_key"], out["sample_key"] = ref, alt, condition_key, sample_key

    # per cell type
    summary = []
    ct_counts = pd.crosstab(obs[sample_key], obs[cluster_key])
    for ct in ct_counts.columns:
        ok = ct_counts.index[ct_counts[ct] >= min_cells]
        g = smeta.loc[ok, condition_key].value_counts() if len(ok) else pd.Series(dtype=int)
        if len(g) < 2 or g.min() < min_samples:
            continue
        mask = (obs[cluster_key] == ct).values
        try:
            cnt_c, meta_c = pseudobulk(mask)
            r_c, s_c, design_c, ng_c = deseq(cnt_c, meta_c, re.sub(r"[^A-Za-z0-9]+", "_", str(ct)).strip("_"))
        except Exception as e:  # noqa: BLE001
            summary.append([ct, len(ok), "—", "—", f"failed ({type(e).__name__})"]); continue
        u, d = s_c[s_c.log2FoldChange > 0], s_c[s_c.log2FoldChange < 0]
        summary.append([ct, len(meta_c), len(u), len(d), ", ".join(u.head(5).index) or "—"])
        if len(s_c) >= 5:
            job.figure(f"pseudobulk_volcano_{re.sub(r'[^A-Za-z0-9]+', '_', str(ct))}", f"{ct}: {alt} vs {ref}", volcano(r_c, s_c, f"{ct} · {len(meta_c)} samples"),
                       how="Same test, restricted to one cell type — the change within that population rather than a change in its abundance.",
                       yours=f"{len(u)} up, {len(d)} down." + (f" Top up: {', '.join(u.head(4).index)}." if len(u) else ""))
    if summary:
        job.table("pseudobulk_by_type", "Pseudobulk DESeq2 per cell type", ["Cell type", "Samples", "Up", "Down", "Top up-regulated"], summary,
                  note=f"Only cell types with ≥{min_cells} cells in ≥{min_samples} samples per group are tested. Per-type CSVs are in the job folder.")
    out["per_type"] = summary
    job.method("Composition", f"Per-sample cell-type percentages compared between {condition_key} groups with a Mann–Whitney test, Bonferroni-adjusted across cell types.")
    job.method("Pseudobulk DE", f"Raw counts summed per {sample_key} (and per {sample_key} × cell type), PyDESeq2 with design {design}, Wald test, BH FDR. Threshold padj &lt; 0.05 and |log2FC| ≥ 1.")
    return out
