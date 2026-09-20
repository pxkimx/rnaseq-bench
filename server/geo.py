"""GEO helpers: parse series-matrix metadata, list and download supplementary files,
and fetch Ensembl-to-symbol tables. Network calls run on the user's machine."""
from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

CACHE = Path.home() / ".rnaseq-bench"
_OLD_CACHE = Path.home() / ".transcript-lens"
if not CACHE.exists() and _OLD_CACHE.exists():
    try:
        _OLD_CACHE.rename(CACHE)          # keeps the API key and Ensembl tables from the old name
    except OSError:
        pass
CACHE.mkdir(exist_ok=True)
UA = {"User-Agent": "RNAseqBench/1.0 (research tool)"}


# ---------------------------------------------------------------- series matrix
def parse_series_matrix(path: Path) -> pd.DataFrame:
    """One row per GSM sample; characteristics become columns (key: value -> column key)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        lines = [l.rstrip("\n") for l in fh if l.startswith("!Sample_")]
    rows: dict[str, list[str]] = {}
    for l in lines:
        tag, _, rest = l.partition("\t")
        vals = [v.strip('"') for v in rest.split("\t")]
        rows.setdefault(tag[len("!Sample_"):], []).append(vals)
    if "geo_accession" not in rows:
        raise ValueError("Not a GEO series matrix (no !Sample_geo_accession line).")
    gsm = rows["geo_accession"][0]
    df = pd.DataFrame(index=gsm)
    df.index.name = "gsm"
    for key in ("title", "source_name_ch1", "organism_ch1", "description"):
        if key in rows:
            df[key.replace("_ch1", "")] = rows[key][0]
    for vals in rows.get("characteristics_ch1", []):
        keys = set()
        for v in vals:
            if ":" in v:
                keys.add(v.split(":", 1)[0].strip())
        for k in keys:
            col = re.sub(r"[^a-z0-9]+", "_", k.lower()).strip("_") or "characteristic"
            df[col] = [v.split(":", 1)[1].strip() if v.startswith(k + ":") else (df[col].loc[g] if col in df else "")
                       for g, v in zip(gsm, vals)]
    if "supplementary_file_1" in rows:
        df["supplementary_file"] = rows["supplementary_file_1"][0]
    return df


def match_metadata_to_columns(meta: pd.DataFrame, columns: list[str]) -> pd.DataFrame | None:
    """Try to line GEO samples up with count-matrix columns: exact title match, a token of the
    title (e.g. plate P014 -> S014), or GSM ids embedded in column names."""
    cols = list(map(str, columns))
    out = pd.DataFrame(index=cols)
    title = meta["title"].astype(str) if "title" in meta else pd.Series("", index=meta.index)
    text = meta.astype(str).agg(" | ".join, axis=1)
    hits = 0
    for c in cols:
        row = None
        exact = text.str.contains(rf"(?<![A-Za-z0-9_]){re.escape(c)}(?![A-Za-z0-9_])", regex=True)
        if c in set(title):
            row = meta[title == c].iloc[0]
        elif exact.sum() == 1:
            row = meta[exact].iloc[0]
        elif (pref := text.str.contains(rf"(?<![A-Za-z0-9_]){re.escape(c)}(?![0-9])", regex=True)).sum() >= 1:
            row = meta[pref].iloc[0]          # e.g. column WT1 vs library names WT1pt1..3 (technical replicates)
        else:
            m = re.search(r"GSM\d+", c)
            if m and m.group(0) in meta.index:
                row = meta.loc[m.group(0)]
            else:
                num = re.search(r"\d{2,}", c)
                if num:
                    cand = meta[title.str.contains(rf"(?<![0-9]){num.group(0)}(?![0-9])", regex=True)]
                    if len(cand) == 1:
                        row = cand.iloc[0]
        if row is not None:
            hits += 1
            for k, v in row.items():
                out.loc[c, k] = v
    return out if hits >= max(2, len(cols) // 2) else None


# ---------------------------------------------------------------- downloads
def _series_dir(acc: str) -> str:
    acc = acc.upper().strip()
    if not re.fullmatch(r"GSE\d+", acc):
        raise ValueError("Accession must look like GSE123456")
    stem = acc[:-3] + "nnn"
    base = os.environ.get("TL_GEO_BASE", "https://ftp.ncbi.nlm.nih.gov/geo/series/")
    return f"{base.rstrip('/')}/{stem}/{acc}/"


def _ssl_ctx():
    """macOS python.org builds ship without root certs; fall back to certifi's bundle."""
    import ssl
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except Exception:  # noqa: BLE001
        return ssl.create_default_context()


def _open(url: str, timeout: int):
    req = urllib.request.Request(url, headers=UA)
    try:
        return urllib.request.urlopen(req, timeout=timeout)
    except urllib.error.URLError as e:
        if "CERTIFICATE" in str(e.reason).upper():
            return urllib.request.urlopen(req, timeout=timeout, context=_ssl_ctx())
        raise


def list_geo_files(acc: str) -> dict:
    base = _series_dir(acc)
    files = []
    for sub in ("suppl/", "matrix/"):
        try:
            html = _open(base + sub, 30).read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            files.append({"name": f"(no {sub} folder on GEO — HTTP {e.code}; the series may not be public yet)", "url": "", "size": ""})
            continue
        except Exception as e:  # noqa: BLE001
            reason = getattr(e, "reason", e)
            files.append({"name": f"(could not reach GEO for {sub}: {type(e).__name__}: {reason})", "url": "", "size": ""})
            continue
        for name, size in re.findall(r'href="([^"?/][^"]*)"[^\n]*?(\d+(?:\.\d+)?[KMG]?)\s*$', html, flags=re.M):
            if name.startswith(("GSE", "GSM")) or name.endswith((".gz", ".tar", ".h5", ".h5ad", ".txt", ".csv", ".tsv")):
                files.append({"name": name, "url": base + sub + name, "size": size})
        if not any(f["url"].startswith(base + sub) for f in files):
            for name in re.findall(r'href="([^"?/][^"]*\.(?:gz|tar|h5|h5ad|txt|csv|tsv|xlsx))"', html):
                files.append({"name": name, "url": base + sub + name, "size": ""})
    return {"accession": acc.upper(), "files": files}


def download_geo_files(urls: list[str], dest: Path, progress=None) -> list[Path]:
    dest.mkdir(parents=True, exist_ok=True)
    out = []
    for i, u in enumerate(urls):
        name = u.rsplit("/", 1)[-1]
        target = dest / name
        if progress:
            progress(i, len(urls), name)
        with _open(u, 120) as r, open(target, "wb") as f:
            total = int(r.headers.get("Content-Length") or 0)
            done, last = 0, 0.0
            while chunk := r.read(1 << 20):
                f.write(chunk)
                done += len(chunk)
                if progress and (done - last) > (8 << 20):
                    last = done
                    progress(i, len(urls), f"{name} — {_mb(done)}{' of ' + _mb(total) if total else ''}")
        out.append(target)
    return out


def _mb(n: int) -> str:
    return f"{n / 1e9:.2f} GB" if n >= 1e9 else f"{n / 1e6:.0f} MB"


# ---------------------------------------------------------------- gene id mapping
MAPPING_URLS = {  # related-sciences/ensembl-genes output branches (Ensembl 104)
    "mouse": "https://raw.githubusercontent.com/related-sciences/ensembl-genes/output/mus_musculus_core_104_39/genes.tsv.gz",
    "human": "https://raw.githubusercontent.com/related-sciences/ensembl-genes/output/homo_sapiens_core_104_38/genes.tsv.gz",
}


# UCSC refGene: RefSeq transcript accession -> gene symbol. Many GEO tables (HOMER/featureCounts output,
# and anything quantified against a RefSeq annotation) are indexed by NM_/NR_ accessions, which every
# gene-set and network resource misses — without this, GSEA and pathway activity silently find nothing.
REFSEQ_URLS = {
    "human": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz",
    "mouse": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/refGene.txt.gz",
}


def refseq_table(species: str) -> pd.Series | None:
    """RefSeq accession (no version) -> symbol. Cached in ~/.rnaseq-bench."""
    p = CACHE / f"refseq_{species}.tsv.gz"
    if not p.exists():
        try:
            with _open(REFSEQ_URLS[species], 120) as r, open(p, "wb") as f:
                shutil.copyfileobj(r, f)
        except Exception:  # noqa: BLE001
            return None
    try:
        t = pd.read_csv(p, sep="\t", header=None, usecols=[1, 12], names=["acc", "symbol"], dtype=str)
        t["acc"] = t["acc"].str.replace(r"\.\d+$", "", regex=True)
        return t.dropna().drop_duplicates("acc").set_index("acc")["symbol"]
    except Exception:  # noqa: BLE001
        p.unlink(missing_ok=True)
        return None


def map_refseq(ids: pd.Index) -> tuple[pd.Series, str]:
    """Symbols aligned to ids, and the species whose table matched best.

    A RefSeq accession does not say which species it belongs to, so both tables are tried and the one
    that explains more of the index wins.
    """
    stripped = ids.str.replace(r"\.\d+$", "", regex=True)
    best, best_species, best_hits = None, "unknown", 0
    for sp in ("human", "mouse"):
        t = refseq_table(sp)
        if t is None:
            continue
        sym = t.reindex(stripped)
        hits = int(sym.notna().sum())
        if hits > best_hits:
            best, best_species, best_hits = sym, sp, hits
    if best is None or best_hits < 0.2 * len(ids):
        return pd.Series(ids, index=ids), "unknown"
    out = pd.Series([s if isinstance(s, str) and s else i for s, i in zip(best.values, ids)], index=ids)
    return out, best_species


def ensembl_table(species: str) -> pd.DataFrame | None:
    """Ensembl stable ID -> symbol / chromosome / biotype. Cached in ~/.rnaseq-bench."""
    p = CACHE / f"ensembl_{species}.tsv.gz"
    if not p.exists():
        try:
            with _open(MAPPING_URLS[species], 120) as r, open(p, "wb") as f:
                shutil.copyfileobj(r, f)
        except Exception:  # noqa: BLE001
            return None
    try:
        t = pd.read_csv(p, sep="\t", dtype=str, usecols=lambda c: c in ("ensembl_gene_id", "gene_symbol", "chromosome", "gene_biotype"))
        t = t.rename(columns={"ensembl_gene_id": "ensembl", "gene_symbol": "symbol", "chromosome": "chrom", "gene_biotype": "biotype"})
        if "ensembl" not in t or "symbol" not in t:
            return None
        return t.drop_duplicates("ensembl").set_index("ensembl")[[c for c in ("symbol", "chrom", "biotype") if c in t]]
    except Exception:  # noqa: BLE001
        return None


def map_ensembl(ids: pd.Index) -> tuple[pd.Series, pd.DataFrame | None, str]:
    """Return (symbols aligned to ids, annotation frame or None, species)."""
    stripped = ids.str.replace(r"\.\d+$", "", regex=True)
    species = "mouse" if stripped.str.startswith("ENSMUSG").mean() > .5 else "human" if stripped.str.startswith("ENSG").mean() > .5 else None
    if species is None:
        return pd.Series(ids, index=ids), None, "unknown"
    t = ensembl_table(species)
    if t is None:
        return pd.Series(ids, index=ids), None, species
    sym = t["symbol"].reindex(stripped).values
    out = pd.Series([s if isinstance(s, str) and s else i for s, i in zip(sym, ids)], index=ids)
    return out, t.reindex(stripped).set_index(ids), species


def save_settings(d: dict):
    p = CACHE / "config.json"
    cur = load_settings()
    cur.update(d)
    p.write_text(json.dumps(cur))


def load_settings() -> dict:
    p = CACHE / "config.json"
    try:
        return json.loads(p.read_text())
    except Exception:  # noqa: BLE001
        return {}
