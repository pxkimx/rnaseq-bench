"""Read user files into AnnData (single-cell) or count DataFrames (bulk)."""
from __future__ import annotations

import gzip
import re
import shutil
import tarfile
import zipfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc

from .common import UserFacingError

SC_EXT = (".h5ad", ".h5", ".loom", ".mtx", ".mtx.gz")
TABLE_EXT = (".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz", ".xlsx", ".xls", ".counts", ".counts.gz", ".tab", ".tab.gz")


def _lower(p: Path) -> str:
    return p.name.lower()


IGNORE_EXT = (".bam", ".bai", ".bw", ".bigwig", ".bed", ".bedgraph", ".bedgraph.gz", ".bed.gz", ".wig", ".fastq", ".fq",
              ".fastq.gz", ".fq.gz", ".pdf", ".png", ".jpg", ".jpeg", ".svg", ".html", ".md5", ".rds", ".rdata", ".rda",
              ".xml", ".json", ".pptx", ".docx", ".gtf", ".gtf.gz", ".gff", ".gff.gz", ".fa", ".fa.gz", ".vcf", ".vcf.gz")


def unpack_archives(files: list[Path], dest: Path, _depth=0) -> list[Path]:
    """Expand .zip/.tar/.tar.gz (recursively: GEO RAW.tar often holds per-sample .tar.gz), drop files the
    pipelines cannot use, and keep gzipped tables as they are (they are read compressed)."""
    out = []
    for f in files:
        n = _lower(f)
        if n.endswith(".zip"):
            d = dest / f"{f.stem}_unzipped"
            with zipfile.ZipFile(f) as z:
                z.extractall(d)
            inner = [p for p in d.rglob("*") if p.is_file()]
            out += unpack_archives(inner, d, _depth + 1) if _depth < 3 else inner
        elif n.endswith((".tar.gz", ".tgz", ".tar")):
            d = dest / (re.sub(r"\.(tar\.gz|tgz|tar)$", "", f.name, flags=re.I) + "_untar")
            with tarfile.open(f) as t:
                t.extractall(d, filter="data")
            inner = [p for p in d.rglob("*") if p.is_file()]
            out += unpack_archives(inner, d, _depth + 1) if _depth < 3 else inner
        else:
            out.append(f)
    keep = []
    for p in out:
        n = _lower(p)
        if p.name.startswith(("._", ".")) or "__macosx" in n or n.endswith(IGNORE_EXT):
            continue
        keep.append(p)
    return sorted(set(keep))


# Header names that are annotation, not samples. The numeric filter catches the text ones; these are the
# numeric annotation columns that would otherwise be analysed as a sample — "Copies" in HOMER/UCSC-style
# tables is the usual offender, and it skews the size factors badly because its scale is unrelated.
ANNOT_COLS = {"chr", "chrom", "chromosome", "start", "end", "strand", "length", "gene_length", "gene_name",
              "gene_symbol", "symbol", "gene_type", "gene_biotype", "biotype", "description", "ensembl",
              "ensembl_id", "entrez", "entrezid", "name", "type", "copies", "exons", "exon_count", "exoncount",
              "transcript", "transcript_id", "transcriptid", "refseq", "accession", "locus", "gene", "gene_id",
              "geneid", "tss", "cds_size", "utr_size", "annotation", "annotation/divergence", "peak score"}


def _peek_table(path: Path, nrows=300):
    """First rows of a text table, with the separator and which columns are numeric counts."""
    n = _lower(path)
    opener = gzip.open if n.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        head = fh.readline()
        while head.startswith("#") or head.startswith("!"):
            head = fh.readline()
    sep = "\t" if head.count("\t") >= head.count(",") else ","
    df = pd.read_csv(path, sep=sep, index_col=0, nrows=nrows, comment="#", low_memory=False, encoding_errors="replace")
    num = df.apply(pd.to_numeric, errors="coerce")
    ok = num.notna().mean() > 0.9
    numeric = [c for c in df.columns if ok[c] and str(c).strip().lower() not in ANNOT_COLS]
    return df, sep, numeric


def guess_kind(files: list[Path]) -> str:
    names = [_lower(f) for f in files]
    if any(n.endswith(SC_EXT) for n in names):
        return "sc"
    try:
        main = pick_counts_file(files)
        df, _, numeric = _peek_table(main)
    except Exception:  # noqa: BLE001
        return "bulk"
    cols = [str(c) for c in numeric]
    if len(cols) >= 500:
        return "sc"
    barcode_like = sum(bool(re.search(r"[ACGTN]{10,}", c)) for c in cols)
    if barcode_like > 0.5 * max(len(cols), 1):
        return "sc"
    if len(cols) < 60:
        return "bulk"
    # 60-500 numeric columns: bulk if the sampled rows look deeply sequenced (per-column sums in the thousands
    # from ~300 genes means millions of reads per library), single-cell if they look like a sparse UMI table.
    sub = df[numeric].apply(pd.to_numeric, errors="coerce").fillna(0)
    frac_zero = float((sub.values == 0).mean())
    per_col = sub.sum(axis=0).median()
    return "sc" if (frac_zero > 0.6 and per_col < 2000) else "bulk"


def is_series_matrix(f: Path) -> bool:
    return "series_matrix" in _lower(f)


def read_series_matrices(files: list[Path]):
    """Sample annotations from every series matrix in the upload, concatenated.

    A GEO series sequenced on more than one instrument ships one series matrix per platform, each
    covering only its own samples. Reading just the first one (they arrive in file-name order, so that
    can easily be the smaller platform) means most samples match nothing and the whole annotation is
    dropped. Returns (frame or None, the files it came from).
    """
    from .geo import parse_series_matrix
    frames, used = [], []
    for f in sorted(x for x in files if is_series_matrix(x)):
        try:
            g = parse_series_matrix(f)
            if g is not None and len(g):
                frames.append(g); used.append(f)
        except Exception:  # noqa: BLE001
            continue
    if not frames:
        return None, []
    out = pd.concat(frames) if len(frames) > 1 else frames[0]
    return out[~out.index.duplicated()], used


def pick_counts_file(files: list[Path]) -> Path:
    tabs = [f for f in files if _lower(f).endswith(TABLE_EXT) and not is_series_matrix(f)
            and not any(k in _lower(f) for k in ("readme", ".summary", "md5", "checksum"))]
    if not tabs:
        raise UserFacingError("No count table found in the files given. Expected a .csv/.tsv/.txt/.xlsx counts matrix (genes × samples), "
                              "a 10x matrix (matrix.mtx + features + barcodes), or an .h5ad / 10x .h5 file. "
                              f"Files seen: {', '.join(sorted(f.name for f in files)[:8]) or 'none'}.")
    meta = [f for f in tabs if any(k in _lower(f) for k in ("meta", "sample", "coldata", "design", "pheno"))]
    counts = [f for f in tabs if f not in meta] or tabs
    return max(counts, key=lambda f: f.stat().st_size)


def pick_metadata_file(files: list[Path], counts: Path) -> Path | None:
    tabs = [f for f in files if f != counts and not is_series_matrix(f) and _lower(f).endswith(TABLE_EXT)
            and not any(k in _lower(f) for k in ("readme", ".summary", "md5", "checksum"))]
    named = [f for f in tabs if any(k in _lower(f) for k in ("meta", "sample", "coldata", "design", "pheno", "annot"))]
    pool = named or tabs
    return min(pool, key=lambda f: f.stat().st_size) if pool else None


def read_table(path: Path, nrows=None) -> pd.DataFrame:
    n = _lower(path)
    if n.endswith((".xlsx", ".xls")):
        df = pd.read_excel(path, index_col=0, nrows=nrows)
        df.columns = [str(c).strip() for c in df.columns]
        return df
    opener = gzip.open if n.endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        head = fh.readline()
        while head.startswith("#") or head.startswith("!"):   # featureCounts / GEO comment lines
            head = fh.readline()
    sep = "\t" if head.count("\t") >= head.count(",") else ","
    df = pd.read_csv(path, sep=sep, index_col=0, nrows=nrows, comment="#", low_memory=False, encoding_errors="replace")
    df.columns = [str(c).strip().strip('"') for c in df.columns]
    df.index = df.index.astype(str).str.strip().str.strip('"')
    return df


def _sample_stem(f: Path) -> str:
    n = re.sub(r"\.(h5ad|h5|loom|csv|tsv|txt|mtx)(\.gz)?$", "", f.name, flags=re.I)
    return re.sub(r"[_\-.]?(filtered_feature_bc_matrix|raw_feature_bc_matrix|counts?|matrix|umi_?counts?|expression)$", "", n, flags=re.I) or n


def _stack(parts: list, notes: list[str], what: str):
    for a in parts:
        a.obs_names = a.obs_names.astype(str)
        a.var_names = a.var_names.astype(str)
        a.var_names_make_unique()
    adata = ad.concat(parts, join="outer", label=None, index_unique="-", fill_value=0)
    adata.obs["sample"] = adata.obs["sample"].astype(str).astype("category")
    if all(p.var_names.equals(parts[0].var_names) for p in parts):
        adata.var = parts[0].var.reindex(adata.var_names)
    notes.append(f"Stacked {len(parts)} {what} ({adata.n_obs:,} cells); the sample name comes from each file's prefix and is stored in obs['sample'].")
    return adata


def _open_text(p: Path):
    return gzip.open(p, "rt") if _lower(p).endswith(".gz") else open(p, "rt")


def _mtx_stem(p: Path) -> str:
    """'GSM123_ctrl_matrix.mtx.gz' -> 'GSM123_ctrl'; 'matrix.mtx' -> '' (or the enclosing folder name when
    the file sits in a per-sample folder such as GSM123_ctrl/filtered_feature_bc_matrix/matrix.mtx.gz)."""
    n = re.sub(r"\.(mtx|tsv|csv|txt)(\.gz)?$", "", p.name, flags=re.I)
    stem = re.sub(r"[_\-.]?(matrix|features|genes|barcodes)$", "", n, flags=re.I)
    if stem:
        return stem
    for parent in p.parents:
        nm = parent.name
        if nm and nm.lower() not in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix", "filtered_gene_bc_matrices",
                                      "raw_gene_bc_matrices", "hg19", "grch38", "mm10", "grcm38", "outs", "input", "untar") \
                and not nm.lower().endswith("_untar"):
            return nm
    return ""


def _read_one_mtx(mtx: Path, feat: Path, bc: Path):
    import scipy.io
    import scipy.sparse as sp
    with (gzip.open(mtx, "rb") if _lower(mtx).endswith(".gz") else open(mtx, "rb")) as fh:
        m = scipy.io.mmread(fh)
    with _open_text(feat) as fh:
        ft = pd.read_csv(fh, sep="\t", header=None, dtype=str)
    with _open_text(bc) as fh:
        bcs = pd.read_csv(fh, sep="\t", header=None, dtype=str)[0].values
    # features: [id, symbol, type] (v3) or [id, symbol] (v2) or [symbol]
    if ft.shape[1] >= 2:
        sym, gid = ft[1].values, ft[0].values
    else:
        sym, gid = ft[0].values, ft[0].values
    x = sp.csr_matrix(m.T if m.shape[0] == len(sym) else m)
    if x.shape != (len(bcs), len(sym)):
        raise UserFacingError(f"{mtx.name}: matrix is {m.shape} but {len(sym)} features and {len(bcs)} barcodes were given.")
    a = ad.AnnData(X=x, obs=pd.DataFrame(index=bcs), var=pd.DataFrame({"gene_ids": gid}, index=sym))
    if a.n_obs > 60000:
        # raw_feature_bc_matrix: millions of empty droplets. Drop barcodes with < 200 UMIs right away
        # (the QC step applies the real thresholds later) so memory stays sane.
        tot = np.asarray(a.X.sum(axis=1)).ravel()
        keep = tot >= 200
        if keep.sum() < a.n_obs:
            a = a[keep].copy()
            a.uns["empty_droplets_removed"] = int((~keep).sum())
    if ft.shape[1] >= 3:
        keep = ft[2].str.contains("Gene Expression", case=False).values
        if keep.any() and not keep.all():
            a = a[:, keep].copy()
    a.var_names_make_unique()
    return a


def _read_mtx_samples(files: list[Path], mtx: list[Path], notes: list[str]) -> ad.AnnData:
    """One or many 10x triplets (matrix/features|genes/barcodes), grouped by file-name stem (e.g. GSM prefix)."""
    def _find(kind, stem):
        cands = [f for f in files if kind in _lower(f) and _mtx_stem(f) == stem]
        return cands[0] if cands else None
    parts = []
    for m in sorted(mtx):
        stem = _mtx_stem(m)
        feat = _find("features", stem) or _find("genes", stem)
        bc = _find("barcodes", stem)
        if feat is None or bc is None:
            if len(mtx) == 1:  # single sample: accept any features/barcodes file
                feat = next((f for f in files if "features" in _lower(f) or "genes" in _lower(f)), None)
                bc = next((f for f in files if "barcodes" in _lower(f)), None)
            if feat is None or bc is None:
                raise UserFacingError(f"{m.name} needs matching features.tsv (or genes.tsv) and barcodes.tsv files "
                                      f"(looked for files starting with '{stem}').")
        a = _read_one_mtx(m, feat, bc)
        a.obs["sample"] = stem or "sample1"
        parts.append(a)
    if len(parts) == 1:
        a = parts[0]
        if not a.obs["sample"].iloc[0] or a.obs["sample"].iloc[0] == "sample1":
            a.obs.drop(columns="sample", inplace=True)
        return a
    adata = _stack(parts, notes, "10x samples")
    return adata


# ------------------------------------------------------------------ single-cell
def load_single_cell(files: list[Path], workdir: Path) -> tuple[ad.AnnData, list[str]]:
    """Returns AnnData with raw integer counts in .X and a list of notes."""
    notes: list[str] = []
    names = {_lower(f): f for f in files}
    h5ad = [f for f in files if _lower(f).endswith(".h5ad")]
    h5 = [f for f in files if _lower(f).endswith(".h5")]
    loom = [f for f in files if _lower(f).endswith(".loom")]
    mtx = [f for f in files if _lower(f).endswith((".mtx", ".mtx.gz"))]

    if h5ad:
        if len(h5ad) > 1:
            parts = []
            for f in sorted(h5ad):
                a = _recover_counts(sc.read_h5ad(f), notes)
                a.obs["sample"] = _sample_stem(f)
                parts.append(a)
            adata = _stack(parts, notes, "h5ad")
        else:
            adata = sc.read_h5ad(h5ad[0])
            adata = _recover_counts(adata, notes)
    elif h5:
        parts = []
        for f in sorted(h5):
            a = sc.read_10x_h5(f)
            a.var_names_make_unique()
            if len(h5) > 1:
                a.obs["sample"] = _sample_stem(f)
            parts.append(a)
        adata = parts[0] if len(parts) == 1 else _stack(parts, notes, "10x .h5")
    elif loom:
        adata = sc.read_loom(loom[0])
        adata = _recover_counts(adata, notes)
    elif mtx:
        adata = _read_mtx_samples(files, mtx, notes)
    else:
        tabs = [f for f in files if _lower(f).endswith((".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz"))
                and not is_series_matrix(f) and not any(k in _lower(f) for k in ("meta", "annot", "coldata", "pheno", "readme"))]
        gsm_tabs = [f for f in tabs if re.match(r"GSM\d+", f.name)]
        if len(gsm_tabs) >= 2 and len(gsm_tabs) == len(tabs):
            # one dense table per GEO sample: read each, stack, label by GSM prefix
            parts = []
            for f in sorted(gsm_tabs):
                a = _read_dense_matrix_streaming(f, notes)
                a.obs["sample"] = _sample_stem(f)
                parts.append(a)
            adata = _stack(parts, notes, "per-sample count tables")
            f = gsm_tabs[0]
        else:
            f = pick_counts_file(files)
            adata = _read_dense_matrix_streaming(f, notes)
        meta = pick_metadata_file(files, f)
        if meta is not None:
            m = read_table(meta)
            m.index = m.index.astype(str)
            shared = adata.obs_names.intersection(m.index)
            if len(shared) > 0.5 * adata.n_obs:
                adata.obs = adata.obs.join(m, how="left")
                notes.append(f"Joined cell metadata from {meta.name} ({m.shape[1]} columns).")

    # GEO series matrix: annotate cells through the sample/plate prefix of their barcode
    gm, sms = read_series_matrices(files)
    if gm is not None:
        sm_names = ", ".join(x.name for x in sms)
        try:
            from .geo import match_metadata_to_columns
            if "sample" in adata.obs and adata.obs["sample"].nunique() > 1:
                # stacked 10x files: the sample name (file prefix, usually containing the GSM id) is the key
                prefix = adata.obs["sample"].astype(str)
                via = "file prefix"
            else:
                prefix = pd.Series([str(c).split("_")[0].split("-")[0] for c in adata.obs_names], index=adata.obs_names)
                via = "barcode prefix"
            matched = match_metadata_to_columns(gm, sorted(prefix.unique()))
            if matched is not None:
                keep_cols = [c for c in matched.columns if c not in ("title", "description", "supplementary_file", "organism")
                             and matched[c].astype(str).str.len().max() < 40]
                for c in keep_cols:
                    adata.obs[c] = prefix.map(matched[c]).values
                adata.obs["sample"] = prefix.values
                notes.append(f"Cell annotations from the GEO series matrix via {via}: {', '.join(keep_cols[:6])}; 'sample' = {via}.")
            else:
                notes.append(f"WARN:{sm_names} was read but its samples could not be matched to the {via}es ({', '.join(sorted(prefix.unique())[:4])}…).")
        except Exception as e:  # noqa: BLE001
            notes.append(f"WARN:Could not use {sm_names}: {type(e).__name__}.")

    # Ensembl IDs -> symbols
    if adata.var_names.str.match(r"^ENS[A-Z]*G\d+").mean() > 0.5:
        from .geo import map_ensembl
        sym, ann, species = map_ensembl(adata.var_names)
        if ann is not None:
            adata.var["ensembl"] = adata.var_names
            adata.var["chrom"] = ann["chrom"].values if "chrom" in ann else ""
            adata.var["biotype"] = ann["biotype"].values if "biotype" in ann else ""
            n_map = int((sym.values != adata.var_names.values).sum())
            adata.var_names = sym.values
            notes.append(f"{n_map:,} Ensembl {species} gene IDs converted to symbols (Ensembl 104).")
        else:
            notes.append("WARN:Gene IDs are Ensembl accessions and no symbol table could be fetched; mitochondrial and marker detection will be limited.")

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    adata.X = adata.X.astype(np.float32)
    if adata.n_obs < 50:
        raise UserFacingError(f"Only {adata.n_obs} cells found. For fewer than 50 columns, run a bulk analysis instead.")
    return adata, notes


def _read_dense_matrix_streaming(f: Path, notes: list[str]) -> ad.AnnData:
    """Read a genes x cells (or cells x genes) dense text matrix row by row into a sparse
    matrix, so a multi-GB table never has to exist densely in memory. Non-integer values
    (RSEM expected counts) are rounded; genes seen in < 3 cells are dropped."""
    n = _lower(f)
    opener = gzip.open if n.endswith(".gz") else open
    with opener(f, "rt") as fh:
        head = fh.readline().rstrip("\n")
        sep = "\t" if head.count("\t") >= head.count(",") else ","
        cols = [c.strip('"') for c in head.split(sep)]
        first = fh.readline().rstrip("\n").split(sep)
        if len(cols) == len(first) - 1:
            cols = ["gene"] + cols
        fh.seek(0); fh.readline()
        names = cols[1:]
        # featureCounts / annotated tables: skip non-numeric annotation columns (Chr, Start, Length, gene_name…)
        keep = np.ones(len(names), dtype=bool)
        probe = first[1:len(names) + 1]
        for j, (cname, cval) in enumerate(zip(names, probe)):
            try:
                float(cval.strip('"'))
            except ValueError:
                keep[j] = False
            if cname.strip().lower() in ANNOT_COLS:
                keep[j] = False
        if not keep.all():
            notes.append(f"Ignored non-count columns: {', '.join(n for n, k in zip(names, keep) if not k)}.")
        names = [n for n, k in zip(names, keep) if k]
        rows, idx, val, indptr, nonint, seen = [], [], [], [0], 0, 0
        for line in fh:
            p = line.rstrip("\n").split(sep)
            if len(p) < 2:
                continue
            v = np.asarray([x.strip('"') or "0" for x, k in zip(p[1:len(keep) + 1], keep) if k], dtype=np.float32)
            if v.size != len(names):
                continue
            if seen < 200000:
                nonint += int(np.sum(np.abs(v - np.round(v)) > 1e-3)); seen += v.size
            v = np.round(v)
            nz = np.flatnonzero(v >= 1)
            if nz.size < 3:
                continue
            rows.append(p[0].strip('"')); idx.append(nz.astype(np.int32)); val.append(v[nz]); indptr.append(indptr[-1] + nz.size)
    X = sp.csr_matrix((np.concatenate(val), np.concatenate(idx), np.array(indptr)), shape=(len(rows), len(names)))
    if seen and nonint / seen > 0.05:
        notes.append("Values were not integers (e.g. RSEM expected counts) and were rounded before analysis.")
    # orientation: we assumed rows = genes. If far more rows than columns and rows look like barcodes, transpose.
    bc = sum(bool(__import__("re").match(r"^[ACGTN]{8,}", r)) for r in rows[:200])
    if bc > 100:
        X = X.T.tocsr(); rows, names = names, rows
        notes.append("Rows looked like cell barcodes, so the table was transposed.")
    else:
        X = X.T.tocsr()          # cells x genes
    return ad.AnnData(X.astype(np.float32), obs=pd.DataFrame(index=pd.Index(names).astype(str)),
                      var=pd.DataFrame(index=pd.Index(rows).astype(str)))


def _recover_counts(adata: ad.AnnData, notes: list[str]) -> ad.AnnData:
    def is_counts(X):
        s = X[: min(200, X.shape[0])]
        s = s.data if sp.issparse(s) else np.asarray(s).ravel()
        s = s[s != 0][:20000]
        return s.size > 0 and np.all(s >= 0) and np.allclose(s, np.round(s))

    if is_counts(adata.X):
        return adata
    for key in ("counts", "raw_counts", "count"):
        if key in adata.layers and is_counts(adata.layers[key]):
            notes.append(f"Used raw counts from layers['{key}'] (X was normalized).")
            new = adata.copy()
            new.X = adata.layers[key]
            return new
    if adata.raw is not None and is_counts(adata.raw.X):
        notes.append("Used raw counts from adata.raw (X was normalized).")
        new = ad.AnnData(adata.raw.X, obs=adata.obs.copy(), var=adata.raw.var.copy(),
                         obsm=dict(adata.obsm))
        return new
    raise UserFacingError(
        "This .h5ad has no raw counts: X is normalized and there is no layers['counts'] or .raw with integers. "
        "The pipeline needs raw UMI counts to do QC and normalization correctly. Re-export with counts saved, "
        "e.g. adata.layers['counts'] = adata.X.copy() before normalizing."
    )


# ------------------------------------------------------------------ bulk
def _merge_per_sample_counts(files: list[Path], notes: list[str]) -> pd.DataFrame | None:
    """GEO RAW.tar for bulk often holds one small file per sample (HTSeq / featureCounts / kallisto output).
    If ≥2 GSM-prefixed tables each carry a single count column, join them into one matrix."""
    tabs = [f for f in files if _lower(f).endswith(TABLE_EXT) and re.match(r"GSM\d+", f.name) and not is_series_matrix(f)]
    if len(tabs) < 2:
        return None
    cols = {}
    for f in sorted(tabs):
        try:
            df = read_table(f)
        except Exception:  # noqa: BLE001
            continue
        num = df.apply(pd.to_numeric, errors="coerce")
        num = num.loc[:, num.notna().mean() > 0.9]
        num = num.drop(columns=[c for c in num.columns if str(c).strip().lower() in ANNOT_COLS | {"tpm", "fpkm", "rpkm", "eff_length", "effective_length", "tpm_value"}], errors="ignore")
        if num.shape[1] == 0 or num.shape[1] > 3:
            return None                       # a real matrix, not per-sample files
        col = "expected_count" if "expected_count" in num.columns else ("est_counts" if "est_counts" in num.columns else num.columns[-1])
        ser = num[col]
        ser = ser[~ser.index.astype(str).str.startswith("__")]   # HTSeq __no_feature etc.
        cols[_sample_stem(f)] = ser
    if len(cols) < 2:
        return None
    mat = pd.concat(cols, axis=1, join="outer").fillna(0)
    notes.append(f"Merged {len(cols)} per-sample count files into one matrix (columns named after each file's GSM prefix).")
    return mat


def load_bulk(files: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame | None, list[str]]:
    notes = []
    merged = _merge_per_sample_counts(files, notes)
    if merged is not None:
        f, df = None, merged
    else:
        f = pick_counts_file(files)
        df = read_table(f)
    # drop annotation columns (gene names, lengths…)
    num = df.apply(pd.to_numeric, errors="coerce")
    num = num.loc[:, num.notna().mean() > 0.9]
    dropped = [c for c in df.columns if c not in num.columns]
    for c in list(num.columns):
        if str(c).strip().lower() in ANNOT_COLS:
            num = num.drop(columns=c)
            dropped.append(c)
    if dropped:
        notes.append(f"Ignored non-count columns: {', '.join(map(str, dropped[:6]))}.")
    if num.shape[0] < 200 and num.shape[1] > 2000:
        num = num.T
        notes.append("The table had samples as rows and genes as columns, so it was transposed.")
    if num.shape[1] < 2:
        raise UserFacingError("Bulk analysis needs a table with at least two sample columns of counts.")
    vals = num.values
    nonint = np.mean(np.abs(vals - np.round(vals)) > 1e-6)
    if nonint > 0.05:
        notes.append("WARN:Values are not integers — this looks like TPM/FPKM or normalized data. DESeq2 needs raw counts; results were computed on rounded values and may be unreliable.")
    num = num.fillna(0).round().astype(int).clip(lower=0)
    num.index = num.index.astype(str)
    num = num[~num.index.duplicated()]
    meta = None
    gm, sms = read_series_matrices(files)
    if gm is not None:
        from .geo import match_metadata_to_columns
        sm_names = ", ".join(x.name for x in sms)
        try:
            matched = match_metadata_to_columns(gm, list(num.columns))
            if matched is not None:
                keep_cols = [c for c in matched.columns if c not in ("title", "description", "supplementary_file", "organism")
                             and matched[c].astype(str).str.len().max() < 40]
                meta = matched[keep_cols]
                notes.append(f"Sample annotations taken from the GEO series matrix ({sm_names}): {', '.join(keep_cols[:6])}.")
            else:
                notes.append(f"WARN:{sm_names} was read but its samples could not be matched to the count columns.")
        except Exception as e:  # noqa: BLE001
            notes.append(f"WARN:Could not parse {sm_names}: {type(e).__name__}.")
    mf = pick_metadata_file(files, f) if f is not None else None
    if mf is not None and meta is None:
        try:
            m = read_table(mf)
        except Exception:  # noqa: BLE001
            m = pd.DataFrame()
        m.index = m.index.astype(str)
        if m.empty:
            notes.append(f"{mf.name} could not be read as a table and was ignored.")
        elif set(num.columns) <= set(m.index):
            meta = m.loc[num.columns]
            notes.append(f"Sample metadata from {mf.name}: columns {', '.join(map(str, m.columns[:6]))}.")
        elif set(num.columns) <= set(m.columns.astype(str)) and m.shape[0] < 30:
            meta = m.T.loc[num.columns]
        else:
            notes.append(f"WARN:{mf.name} was not used — its row names don't match the count column names.")
    return num, meta, notes


def collapse_technical_replicates(counts: pd.DataFrame, meta: pd.DataFrame | None, notes: list[str], force: bool | None = None):
    """Sum columns that are technical replicates of one library (names differing only by a
    trailing rep/pt/lane index). Auto mode collapses only when within-group correlation is
    clearly above between-group correlation."""
    import re as _re
    stem = [_re.sub(r"[_\-. ]?(rep|r|pt|lane|l|tech|t)\d+$", "", c, flags=_re.I) for c in counts.columns]
    groups = pd.Series(stem, index=counts.columns)
    sizes = groups.value_counts()
    if force is False or sizes.max() < 2 or sizes.min() < 2 and force is None and (sizes >= 2).sum() < len(sizes) * 0.8:
        return counts, meta, notes, False
    if force is None:
        lc = np.log2(counts / counts.sum() * 1e6 + 1)
        lc = lc[lc.var(axis=1) > 0]
        corr = lc.corr()
        within, between = [], []
        for c in counts.columns:
            same = [o for o in counts.columns if o != c and groups[o] == groups[c]]
            other = [o for o in counts.columns if groups[o] != groups[c]]
            if same:
                within.append(corr.loc[c, same].mean())
            if other:
                between.append(corr.loc[c, other].max())
        if not within or np.median(within) <= np.median(between) + 0.002:
            return counts, meta, notes, False
        notes.append(f"Columns ending in a replicate index correlate at r={np.median(within):.4f} within a library vs "
                     f"{np.median(between):.4f} to the closest other library, so they were treated as technical replicates and summed "
                     f"({counts.shape[1]} columns -> {sizes.size} samples). Turn this off in settings if they are biological replicates.")
    collapsed = counts.T.groupby(groups.values, sort=False).sum().T
    if meta is not None:
        meta = meta.groupby(groups.reindex(meta.index).values, sort=False).first()
        meta = meta.reindex(collapsed.columns)
    return collapsed, meta, notes, True


def infer_sex(counts: pd.DataFrame) -> pd.Series | None:
    """Expression-based sex call from XIST vs Y-linked genes (human or mouse symbols)."""
    idx = {g.upper(): g for g in counts.index}
    x = idx.get("XIST")
    ys = [idx[g] for g in ("RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF2S3Y") if g in idx]
    if x is None or len(ys) < 2:
        return None
    cpm = counts / counts.sum() * 1e6
    xi = np.log2(cpm.loc[x] + 1)
    yi = np.log2(cpm.loc[ys].sum() + 1)
    call = pd.Series(np.where(yi > xi, "M", "F"), index=counts.columns, name="sex_inferred")
    margin = (yi - xi).abs()
    if (margin < 1).any():
        call[margin < 1] = "ambiguous"
    return call


def map_gene_ids(counts: pd.DataFrame, notes: list[str]) -> pd.DataFrame:
    """Replace Ensembl or RefSeq accessions with gene symbols.

    Everything downstream that knows biology — marker panels, gene sets, PROGENy, CollecTRI — is keyed
    on symbols, so an unmapped index does not fail loudly, it just quietly finds nothing.
    """
    if counts.index.str.match(r"^[NX][MR]_\d+").mean() > 0.5:
        from .geo import map_refseq
        sym, species = map_refseq(counts.index)
        if species == "unknown":
            notes.append("WARN:Gene IDs are RefSeq accessions and no symbol table was available (needs internet "
                         "once); gene sets, pathway activity and marker-based interpretation will find nothing.")
            return counts
        n_mapped = int((sym.values != counts.index.values).sum())
        out = counts.copy()
        out.index = sym.values
        out = out.groupby(level=0).sum()
        notes.append(f"{n_mapped:,} RefSeq {species} transcript IDs converted to gene symbols (UCSC refGene); "
                     "transcripts of the same gene were summed.")
        return out
    if not counts.index.str.match(r"^ENS[A-Z]*G\d+").mean() > 0.5:
        return counts
    from .geo import map_ensembl
    sym, _, species = map_ensembl(counts.index)
    if species == "unknown" or (sym.values == counts.index.values).all():
        notes.append("WARN:Gene IDs are Ensembl accessions and no symbol table was available (needs internet once); marker-based interpretation will be limited.")
        return counts
    n_mapped = int((sym.values != counts.index.values).sum())
    out = counts.copy()
    out.index = sym.values
    out = out.groupby(level=0).sum()
    notes.append(f"{n_mapped:,} Ensembl {species} gene IDs converted to symbols (Ensembl 104); duplicates summed.")
    return out


def infer_groups(names: list[str]) -> tuple[list[str], str]:
    """Guess condition labels from sample names. Returns (labels, how it was guessed).

    `how` matters: a trailing-replicate or prefix pattern is a real signal, but the A/B fallback is just
    the sample order cut in half and means nothing biologically. The caller must not treat the two alike.
    """
    import re
    strip = [re.sub(r"[_\-. ]?(rep|r|s|sample|n)?\d+$", "", n, flags=re.I) for n in names]
    if 1 < len(set(strip)) < len(names):
        return strip, "suffix"
    first = [re.split(r"[_\-. ]", n)[0] for n in names]
    if 1 < len(set(first)) < len(names):
        return first, "prefix"
    h = (len(names) + 1) // 2
    return ["A"] * h + ["B"] * (len(names) - h), "arbitrary"


# ---------------------------------------------------------------- metadata hygiene
_YM = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*(?:y|yr|yrs|year|years)\s*(?:(\d+(?:\.\d+)?)\s*(?:m|mo|mos|month|months))?\s*$", re.I)
# A number optionally followed by a unit. Deliberately no letter *prefix*: "GM08399", "AG09599" and "S1"
# are sample and cell-line identifiers, and reading them as the numbers 8399, 9599 and 1 would offer an
# identifier as a continuous covariate.
_NUM_UNIT = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*([A-Za-z%]*)\s*$")
_UNITS = {"", "y", "yr", "yrs", "year", "years", "m", "mo", "mos", "month", "months", "w", "wk", "wks",
          "week", "weeks", "d", "day", "days", "h", "hr", "hrs", "hour", "hours", "min", "mins", "sec",
          "dpi", "dpf", "hpf", "dpc", "passage", "p", "ng", "ug", "mg", "g", "kg", "ml", "ul", "l",
          "nm", "um", "mm", "cm", "gy", "%", "percent", "x", "fold", "c", "k"}


def numeric_with_units(col: pd.Series) -> pd.Series | None:
    """Parse a column like age that carries units — '8yr', '2yr3mos', '45 years', '30' — into numbers.

    GEO writes ages however the submitter felt: one non-numeric entry is enough for pd.to_numeric to give
    up on the whole column, and a continuous variable then silently becomes a 77-level categorical that
    nothing can use. Years and months are combined into fractional years; any other unit is kept as its
    own number, provided every value that parses agrees on the unit. Returns None when this is not a
    numeric column after all.
    """
    s = col.astype(str).str.strip()
    plain = pd.to_numeric(s, errors="coerce")
    if plain.notna().mean() >= 0.99:
        return plain
    ym = s.str.extract(_YM)
    if ym[0].notna().mean() > 0:
        v = pd.to_numeric(ym[0], errors="coerce") + pd.to_numeric(ym[1], errors="coerce").fillna(0) / 12
        out = plain.where(plain.notna(), v)
        if out.notna().mean() >= 0.8:
            return out
    hit = s.str.extract(_NUM_UNIT)
    units = {u.lower() for u in hit[1].dropna() if u}
    if hit[0].notna().mean() >= 0.8 and len(units) <= 1 and units <= _UNITS:
        return pd.to_numeric(hit[0], errors="coerce")
    return None


def unify_level_case(col: pd.Series) -> pd.Series:
    """Merge levels that differ only in case or surrounding space ('Male' / 'male' / ' Male ').

    GEO characteristics are free text, so the same group routinely appears under several spellings. Left
    alone they become separate levels: a two-level covariate silently costs three degrees of freedom, and
    a design factor gets split into groups that are not really different.
    """
    s = col.astype(str).str.strip()
    key = s.str.casefold()
    if key.nunique() >= s.nunique():
        return col
    best = s.groupby(key).agg(lambda x: x.value_counts().idxmax())
    return key.map(best).where(col.notna(), col)
