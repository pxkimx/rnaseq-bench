"""Expert-level add-ons for the bulk pipeline: GSEA (pre-ranked), pathway / TF activity (decoupler),
genes-of-interest panel, covariate-adjusted PCA. Each function adds figures/tables to the job and
returns True on success; failures are logged and reported as notes, never fatal."""
from __future__ import annotations

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from . import resources
from .common import DOWN, NS, UP, Job, log_exc

ORG_COLORS = {"Hallmark": "#2b7a78", "GO BP": "#8c564b", "KEGG": "#9467bd", "Reactome": "#ff7f0e"}


def _clean(term: str) -> str:
    t = re.sub(r"\s*\(GO:\d+\)|\s*R-HSA-\d+|\s*WP\d+", "", term)
    t = re.sub(r"^HALLMARK_", "", t).replace("_", " ")
    return t[:60]


# ---------------------------------------------------------------- GSEA
def gsea(job: Job, res: pd.DataFrame, organism: str, alt: str, ref: str) -> bool:
    """Pre-ranked GSEA on the Wald statistic against Hallmark, GO BP, KEGG and Reactome."""
    try:
        import gseapy as gp
        libs = resources.all_gene_sets(organism)
        if not libs:
            job.flag("info", "GSEA skipped — gene-set libraries could not be downloaded (needs internet once; they are cached afterwards).")
            return False
        r = res.dropna(subset=["stat"]).copy()
        r.index = r.index.astype(str)
        if organism == "mouse":
            r.index = r.index.str.upper()
        r = r[~r.index.duplicated()]
        rnk = r["stat"].sort_values(ascending=False)
        rnk.to_csv(job.root / "ranked_genes.rnk", sep="\t", header=False)
        frames = []
        for short, sets in libs.items():
            pre = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=1000, min_size=15, max_size=500, threads=4,
                             seed=0, outdir=None, no_plot=True, verbose=False)
            d = pre.res2d.copy()
            d["library"] = short
            d["_pre"] = short
            frames.append(d)
            if short == "Hallmark":
                hall = pre
        gr = pd.concat(frames, ignore_index=True)
        for c in ("NES", "FDR q-val", "NOM p-val", "ES"):
            gr[c] = pd.to_numeric(gr[c], errors="coerce")
        gr = gr.rename(columns={"FDR q-val": "FDR", "NOM p-val": "pval", "Lead_genes": "leading_edge"})
        gr["term"] = gr["Term"].map(_clean)
        keep = ["library", "Term", "term", "NES", "ES", "pval", "FDR", "Gene %", "leading_edge"]
        gr = gr[[c for c in keep if c in gr.columns]].sort_values("FDR")
        gr.to_csv(job.root / "gsea_results.csv", index=False)
        job.download("GSEA results (.csv)", "gsea_results.csv")
        job.download("Ranked gene list for GSEA (.rnk)", "ranked_genes.rnk")

        sig = gr[gr.FDR < 0.25]
        up = sig[sig.NES > 0].sort_values("NES", ascending=False)
        dn = sig[sig.NES < 0].sort_values("NES")
        job.section("gsea", "Step 5", f"Gene-set enrichment (GSEA): {alt} vs {ref}",
                    "Pre-ranked GSEA (Subramanian 2005) on every tested gene ordered by the DESeq2 Wald statistic — no fold-change or "
                    "p-value cut-off, so coordinated small shifts across a pathway are detected. NES > 0 = higher in "
                    f"{alt}. FDR < 0.25 is the conventional GSEA threshold; FDR < 0.05 is strong.")
        fig, axs = plt.subplots(1, 2, figsize=(12, 5.6))
        for ax, d, ttl in zip(axs, (up.head(15), dn.head(15)), (f"higher in {alt}", f"lower in {alt}")):
            if d.empty:
                ax.text(0.5, 0.5, "no gene set at FDR < 0.25", ha="center", va="center", transform=ax.transAxes, color="#777"); ax.axis("off"); continue
            d = d.iloc[::-1]
            ax.barh(range(len(d)), d.NES, color=[ORG_COLORS.get(l, "#888") for l in d.library], alpha=0.9)
            ax.set_yticks(range(len(d))); ax.set_yticklabels([f"{t}  [{l}]" for t, l in zip(d.term, d.library)], fontsize=7.5)
            ax.set_xlabel("normalized enrichment score (NES)"); ax.set_title(ttl)
            for i, (q, n) in enumerate(zip(d.FDR, d.NES)):
                ax.text(n, i, f" q={q:.2g}" if q >= 1e-3 else " q<0.001", va="center", ha="left" if n > 0 else "right", fontsize=6.5, color="#333")
        fig.tight_layout()
        tu = up.head(3).term.tolist(); td = dn.head(3).term.tolist()
        job.figure("gsea_bars", "Top enriched gene sets", fig, wide=True,
                   how="Bar length = NES (enrichment strength normalized for set size), colour = library. A set whose genes sit "
                       "at the top of the ranking gets a positive NES, at the bottom a negative one.",
                   yours=(f"Up in {alt}: <b>{'; '.join(tu)}</b>. " if tu else f"Nothing enriched among genes up in {alt}. ")
                         + (f"Down: <b>{'; '.join(td)}</b>." if td else "Nothing enriched among down genes."))
        # classic enrichment plots for the two strongest Hallmark hits
        try:
            top2 = gr[(gr.library == "Hallmark")].sort_values("FDR").head(2)
            if len(top2):
                figs = []
                for _, row in top2.iterrows():
                    rr = hall.results[row.Term]
                    axs_ = gp.gseaplot(rank_metric=hall.ranking, term=_clean(row.Term), hits=rr["hits"], nes=rr["nes"], pval=rr["pval"], fdr=rr["fdr"], RES=rr["RES"], ofname=None, figsize=(4.8, 4.6))
                    f_ = axs_[0].figure if isinstance(axs_, (list, tuple)) else plt.gcf()
                    figs.append((row, f_))
                for i, (row, f_) in enumerate(figs):
                    job.figure(f"gsea_plot{i}", f"Enrichment plot: {_clean(row.Term)}", f_,
                               how="Top: running enrichment score as the ranked list is walked from most-up to most-down. Middle: where the set's genes fall. "
                                   "A peak near the left means the set is concentrated among up-regulated genes.",
                               yours=f"NES = <b>{row.NES:.2f}</b>, FDR = {row.FDR:.2g}; leading edge: {str(row.leading_edge)[:90]}…")
        except Exception:  # noqa: BLE001
            log_exc("gsea plots")
        rows = [[x.library, x.term, f"{x.NES:.2f}", f"{x.FDR:.2g}", str(x.leading_edge)[:70]] for _, x in gr.head(30).iterrows()]
        job.table("gsea_table", "Top 30 gene sets by FDR", ["Library", "Gene set", "NES", "FDR", "Leading-edge genes"], rows,
                  note="All sets in gsea_results.csv. Leading edge = the genes that drive the enrichment; those are the ones to look at.", csv="gsea_results.csv")
        n_sig = int((gr.FDR < 0.05).sum())
        job.flag("ok" if n_sig else "info", f"GSEA: <b>{n_sig}</b> gene sets at FDR < 0.05 ({int((gr.FDR < 0.25).sum())} at FDR < 0.25)."
                 + (f" Strongest: {gr.iloc[0].term} (NES {gr.iloc[0].NES:.2f})." if len(gr) else ""))
        job.method("GSEA", f"GSEApy prerank on the Wald statistic; libraries {', '.join(libs)}; 1000 permutations, set size 15–500; FDR by gene-set permutation.")
        return True
    except LookupError:
        log_exc("gsea")
        job.flag("info", "GSEA skipped — fewer than 15 genes of any gene set were found among your genes. Usually the gene IDs are not symbols "
                         "(Ensembl IDs are converted automatically when the mapping table is available) or the species was mis-detected.")
        return False
    except Exception:  # noqa: BLE001
        job.flag("info", f"GSEA skipped ({log_exc('gsea')}).")
        return False


# ---------------------------------------------------------------- pathway + TF activity
def activity(job: Job, res: pd.DataFrame, vst: pd.DataFrame, meta: pd.DataFrame, factor: str, ref: str, alt: str, organism: str) -> bool:
    """decoupler ULM: PROGENy pathway and CollecTRI TF activities from the contrast statistic and per sample."""
    try:
        import decoupler as dc
        prog = resources.progeny(organism); tri = resources.collectri(organism)
        if prog is None and tri is None:
            job.flag("info", "Pathway/TF activity skipped — PROGENy/CollecTRI networks could not be downloaded (needs internet once).")
            return False
        stat = res["stat"].dropna()
        stat.index = stat.index.astype(str)
        stat = stat[~stat.index.duplicated()]
        mat = pd.DataFrame([stat.values], index=[f"{alt} vs {ref}"], columns=stat.index)
        job.section("activity", "Step 6", "Pathway & transcription-factor activity",
                    "Instead of asking which genes changed, infer what is driving them: PROGENy scores 14 signalling pathways from their "
                    "downstream responsive genes, and CollecTRI scores transcription factors from their known targets (decoupler, "
                    "univariate linear model on the DESeq2 statistic). Positive = more active in "
                    f"{alt}. This is footprint-based, so a pathway can be active even when its own components are not differentially expressed.")
        done = False
        if prog is not None:
            es, pv = dc.mt.ulm(data=mat, net=prog, tmin=5)
            s = es.iloc[0].sort_values(); p = pv.iloc[0].reindex(s.index)
            fig, ax = plt.subplots(figsize=(5.8, 4.2))
            ax.barh(range(len(s)), s.values, color=[UP if v > 0 else DOWN for v in s.values], alpha=0.9)
            ax.set_yticks(range(len(s))); ax.set_yticklabels(s.index, fontsize=8.5)
            ax.axvline(0, color="#444", lw=0.8); ax.set_xlabel("activity score (ULM t-value)")
            for i, (v, q) in enumerate(zip(s.values, p.values)):
                if q < 0.05:
                    ax.text(v, i, " *", va="center", ha="left" if v > 0 else "right", fontsize=9)
            ax.set_title(f"PROGENy pathways, {alt} vs {ref}")
            fig.tight_layout()
            top_up = [f"{k} ({v:.1f})" for k, v in s[::-1].head(3).items() if v > 0 and p[k] < 0.05]
            top_dn = [f"{k} ({v:.1f})" for k, v in s.head(3).items() if v < 0 and p[k] < 0.05]
            job.figure("progeny_contrast", "Pathway activity (PROGENy)", fig,
                       how="Each pathway's score is the slope of a regression of the DE statistic on the pathway's gene weights; * = p < 0.05.",
                       yours=(f"More active in {alt}: <b>{', '.join(top_up)}</b>. " if top_up else "") + (f"Less active: <b>{', '.join(top_dn)}</b>." if top_dn else "") or "No pathway reaches p < 0.05.")
            # per-sample activities on VST
            v = vst.T.copy(); v.columns = v.columns.astype(str)
            if organism == "mouse":
                v.columns = v.columns.str.upper()
            v = v.loc[:, ~v.columns.duplicated()]
            v = v - v.mean(axis=0)
            es_s, _ = dc.mt.ulm(data=v, net=prog, tmin=5)
            es_s.to_csv(job.root / "pathway_activity_per_sample.csv")
            order = meta.sort_values(factor).index
            es_o = es_s.reindex(order)
            fig, ax = plt.subplots(figsize=(max(5, 0.32 * len(order) + 2.5), 4.6))
            im = ax.imshow(es_o.T.values, aspect="auto", cmap="RdBu_r", vmin=-np.nanpercentile(np.abs(es_o.values), 98), vmax=np.nanpercentile(np.abs(es_o.values), 98))
            ax.set_yticks(range(es_o.shape[1])); ax.set_yticklabels(es_o.columns, fontsize=8)
            ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=90, fontsize=6.5)
            grp = meta.loc[order, factor].astype(str).values
            for i in range(1, len(grp)):
                if grp[i] != grp[i - 1]:
                    ax.axvline(i - 0.5, color="black", lw=1)
            for g in pd.unique(grp):
                idx = np.where(grp == g)[0]
                ax.text(idx.mean(), -0.9, g, ha="center", va="bottom", fontsize=8, fontweight="bold")
            plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="activity")
            fig.tight_layout()
            # per-pathway group test
            tests = []
            for pw in es_s.columns:
                a_ = es_s.loc[meta.index[meta[factor].astype(str) == alt], pw].dropna(); b_ = es_s.loc[meta.index[meta[factor].astype(str) == ref], pw].dropna()
                if len(a_) >= 2 and len(b_) >= 2:
                    t, pp = stats.ttest_ind(a_, b_, equal_var=False)
                    tests.append([pw, f"{a_.mean():.2f}", f"{b_.mean():.2f}", f"{t:.2f}", f"{pp:.2g}"])
            job.figure("progeny_samples", "Pathway activity per sample", fig, wide=True,
                       how="Same 14 pathways scored on each sample's variance-stabilized expression (centred per gene). Consistent colour within a group and a switch between groups is what a real pathway shift looks like.",
                       yours="Per-sample scores are in pathway_activity_per_sample.csv; group comparison in the table below.")
            if tests:
                job.table("progeny_table", f"Pathway activity by group (Welch t-test, {alt} vs {ref})", ["Pathway", f"mean {alt}", f"mean {ref}", "t", "p"], tests,
                          note="Uncorrected p-values across 14 pathways; treat p < 0.004 (Bonferroni) as convincing.")
            job.download("Pathway activity per sample (.csv)", "pathway_activity_per_sample.csv")
            done = True
        if tri is not None:
            es, pv = dc.mt.ulm(data=mat, net=tri, tmin=5)
            s = es.iloc[0]; p = pv.iloc[0]
            df = pd.DataFrame({"TF": s.index, "score": s.values, "p": p.values}).sort_values("score")
            df.to_csv(job.root / "tf_activity.csv", index=False)
            job.download("TF activity (.csv)", "tf_activity.csv")
            top = pd.concat([df.head(12), df.tail(12)])
            fig, ax = plt.subplots(figsize=(6.2, 5.4))
            ax.barh(range(len(top)), top.score, color=[UP if v > 0 else DOWN for v in top.score], alpha=0.9)
            ax.set_yticks(range(len(top))); ax.set_yticklabels(top.TF, fontsize=8)
            ax.axvline(0, color="#444", lw=0.8); ax.set_xlabel("TF activity score"); ax.set_title(f"Transcription factors, {alt} vs {ref}")
            fig.tight_layout()
            sig_up = df[(df.p < 0.05) & (df.score > 0)].nlargest(5, "score").TF.tolist()
            sig_dn = df[(df.p < 0.05) & (df.score < 0)].nsmallest(5, "score").TF.tolist()
            job.figure("tf_activity", "Transcription-factor activity (CollecTRI)", fig,
                       how="A TF scores high when its known activated targets are up and repressed targets are down. It infers activity from targets, so it can flag TFs whose own mRNA does not change.",
                       yours=(f"Likely more active in {alt}: <b>{', '.join(sig_up)}</b>. " if sig_up else "") + (f"Less active: <b>{', '.join(sig_dn)}</b>." if sig_dn else ""))
            if sig_up or sig_dn:
                job.flag("ok", f"TF activity — up in {alt}: {', '.join(sig_up[:3]) or 'none'}; down: {', '.join(sig_dn[:3]) or 'none'}.")
            done = True
        job.method("Activity inference", "decoupler ULM (univariate linear model) with PROGENy (top 500 genes/pathway) and CollecTRI regulons on the Wald statistic (contrast) and on centred VST values (per sample).")
        return done
    except AssertionError:
        log_exc("activity")
        job.flag("info", "Pathway/TF activity skipped — too few network genes were found among your genes (gene IDs not symbols, or species mis-detected).")
        return False
    except Exception:  # noqa: BLE001
        job.flag("info", f"Pathway/TF activity skipped ({log_exc('activity')}).")
        return False


# ---------------------------------------------------------------- genes of interest
def genes_of_interest(job: Job, genes: list[str], norm: pd.DataFrame, meta: pd.DataFrame, factor: str, res: pd.DataFrame, ref: str, alt: str) -> bool:
    try:
        upper = {g.upper(): g for g in norm.index.astype(str)}
        found = [upper[g.strip().upper()] for g in genes if g.strip().upper() in upper]
        missing = [g for g in genes if g.strip() and g.strip().upper() not in upper]
        if not found:
            job.flag("warn", f"None of the requested genes were found in the count table: {', '.join(genes[:8])}.")
            return False
        found = found[:16]
        levels = [ref] + [l for l in meta[factor].astype(str).unique() if l != ref]
        n = len(found); cols = min(4, n); rows_ = int(np.ceil(n / cols))
        fig, axs = plt.subplots(rows_, cols, figsize=(3.1 * cols, 2.8 * rows_), squeeze=False)
        rng = np.random.default_rng(0)
        for ax, g in zip(axs.ravel(), found):
            vals = np.log2(norm.loc[g].astype(float) + 1)
            data = [vals[meta.index[meta[factor].astype(str) == l]].values for l in levels]
            ax.boxplot(data, widths=0.5, showfliers=False, medianprops=dict(color="#222"))
            for i, d in enumerate(data):
                ax.scatter(1 + i + rng.uniform(-0.15, 0.15, len(d)), d, s=14, color=[UP, DOWN, "#9467bd", "#ff7f0e", "#2b7a78"][i % 5], alpha=0.85, zorder=3)
            ax.set_xticks(range(1, len(levels) + 1)); ax.set_xticklabels(levels, fontsize=7.5)
            if g in res.index:
                r = res.loc[g]
                ax.set_title(f"{g}  log2FC {r.log2FoldChange:.2f}, padj {r.padj:.2g}" if pd.notna(r.padj) else f"{g}  (not tested)", fontsize=8)
            else:
                ax.set_title(g, fontsize=8)
            ax.set_ylabel("log2(norm counts + 1)", fontsize=7)
        for ax in axs.ravel()[n:]:
            ax.axis("off")
        fig.tight_layout()
        job.section("goi", "Genes of interest", "Your genes",
                    "Normalized counts (DESeq2 size-factor corrected) for the genes you asked about, every sample shown. The title carries the DESeq2 result for the main contrast.")
        job.figure("genes_of_interest", "Requested genes", fig, wide=n > 4,
                   how="Boxes show the group distribution; dots are individual samples, so you can see whether a fold change is driven by all replicates or by one.",
                   yours=(f"Not found: {', '.join(missing[:6])}. " if missing else "") + " ".join(
                       f"<b>{g}</b>: {'up' if res.loc[g].log2FoldChange > 0 else 'down'} {abs(res.loc[g].log2FoldChange):.1f} log2 (padj {res.loc[g].padj:.2g})." for g in found if g in res.index and pd.notna(res.loc[g].padj) and res.loc[g].padj < 0.05)[:400] or "None of them is significant in the main contrast.")
        return True
    except Exception:  # noqa: BLE001
        job.flag("info", f"Genes-of-interest panel skipped ({log_exc('goi')}).")
        return False


# ---------------------------------------------------------------- covariate-adjusted PCA
def adjusted_pca(job: Job, vst: pd.DataFrame, meta: pd.DataFrame, factor: str, covars: list[str], colors: dict) -> bool:
    """Regress the covariates out of the VST matrix (limma::removeBatchEffect style, keeping the factor) and redo PCA."""
    if not covars:
        return False
    try:
        import formulaic
        X = formulaic.model_matrix("~ " + " + ".join(covars + [factor]), meta.astype({c: str for c in covars + [factor] if meta[c].dtype == object or str(meta[c].dtype) == "category"}))
        X = np.asarray(X, dtype=float)
        cov_cols = [i for i, c in enumerate(X.columns if hasattr(X, "columns") else range(X.shape[1]))]
        Xm = formulaic.model_matrix("~ " + " + ".join(covars + [factor]), meta)
        names = list(Xm.columns); Xa = np.asarray(Xm, dtype=float)
        is_cov = np.array([any(n.startswith(c) for c in covars) for n in names])
        Y = vst.loc[:, meta.index].values.T                        # samples x genes
        beta, *_ = np.linalg.lstsq(Xa, Y, rcond=None)
        Yadj = Y - Xa[:, is_cov] @ beta[is_cov]
        top = np.argsort(Yadj.var(axis=0))[::-1][:500]
        Z = Yadj[:, top] - Yadj[:, top].mean(axis=0)
        U, S_, Vt = np.linalg.svd(Z, full_matrices=False)
        pcs = U * S_; pvar = S_**2 / (S_**2).sum()
        fig, ax = plt.subplots(figsize=(5.4, 4.4))
        for i, s_ in enumerate(meta.index):
            ax.scatter(pcs[i, 0], pcs[i, 1], color=colors.get(str(meta.loc[s_, factor]), "#888"), s=70, edgecolor="white", linewidth=0.8)
            if len(meta) <= 24:
                ax.annotate(s_, (pcs[i, 0], pcs[i, 1]), fontsize=6.5, xytext=(4, 3), textcoords="offset points", color="#444")
        ax.set_xlabel(f"PC1 ({100*pvar[0]:.1f}%)"); ax.set_ylabel(f"PC2 ({100*pvar[1]:.1f}%)")
        ax.set_title(f"PCA after removing {', '.join(covars)}")
        ax.grid(alpha=0.25)
        fig.tight_layout()
        # how much of PC1 does the factor explain now?
        grp = meta[factor].astype(str).values
        gm = pd.Series(pcs[:, 0]).groupby(grp).transform("mean").values
        r2 = 1 - ((pcs[:, 0] - gm) ** 2).sum() / ((pcs[:, 0] - pcs[:, 0].mean()) ** 2).sum()
        job.figure("pca_adjusted", "Covariate-adjusted PCA", fig,
                   how=f"The linear effect of {', '.join(covars)} was regressed out of every gene (as limma's removeBatchEffect does, while protecting {factor}) and PCA recomputed. Use this only for visualization — DESeq2 already accounts for the covariates in its model.",
                   yours=f"<b>{factor}</b> explains {100*r2:.0f}% of PC1 once {', '.join(covars)} is removed.")
        return True
    except Exception:  # noqa: BLE001
        log_exc("adjusted pca")
        return False
