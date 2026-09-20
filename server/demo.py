"""Example datasets. PBMC 3k is downloaded by Scanpy when online; otherwise a realistic
simulation with the same cell types is generated."""
from __future__ import annotations

import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


def pbmc3k(dest: Path) -> tuple[Path, str]:
    dest.mkdir(parents=True, exist_ok=True)
    out = dest / "pbmc3k.h5ad"
    cached = EXAMPLES / "pbmc3k_raw.h5ad"
    if cached.exists():
        shutil.copy(cached, out)
        return out, "10x PBMC 3k (Scanpy / Seurat tutorial data)"
    try:
        import scanpy as sc
        a = sc.datasets.pbmc3k()
        a.write_h5ad(cached)
        shutil.copy(cached, out)
        return out, "10x PBMC 3k (Scanpy / Seurat tutorial data)"
    except Exception:
        simulate_pbmc().write_h5ad(out)
        return out, "Simulated PBMCs, 6 donors, 3 IFN-stimulated (offline example)"


def airway(dest: Path) -> tuple[list[Path], str]:
    dest.mkdir(parents=True, exist_ok=True)
    files = []
    for f in ("GSE96870_counts.csv", "GSE96870_metadata.csv"):
        shutil.copy(EXAMPLES / f, dest / f)
        files.append(dest / f)
    return files, "GSE96870 mouse cerebellum, influenza vs uninfected"


def _background_genes(n: int, taken: set, rng) -> list[str]:
    """Real human protein-coding symbols when the cached Ensembl table is available (so pathway, ligand–receptor
    and CellTypist resources find genes in simulated data), otherwise placeholder names."""
    try:
        from .geo import ensembl_table
        t = ensembl_table("human")
        sym = t["symbol"].dropna().astype(str)
        if "biotype" in t:
            sym = sym[t["biotype"].astype(str).str.contains("protein_coding").values]
        pool = sorted(set(sym) - taken - {s for s in sym if s.startswith(("MT-", "RPL", "RPS"))})
        if len(pool) >= n:
            return list(rng.choice(pool, n, replace=False))
    except Exception:  # noqa: BLE001
        pass
    return [f"GENE{i:05d}" for i in range(n)]


def simulate_pbmc(seed=7) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    types = {
        "CD4 T": (900, ["IL7R", "CCR7", "LTB", "CD3E", "CD3D", "TCF7", "MAL", "LEF1", "CD4"]),
        "CD8 T": (380, ["CD8A", "CD8B", "GZMK", "CD3E", "CD3D", "CCL5", "GZMH"]),
        "B": (360, ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19", "HLA-DQA1", "PAX5"]),
        "CD14 Mono": (520, ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "VCAN", "CST3"]),
        "NK": (170, ["GNLY", "NKG7", "PRF1", "GZMB", "KLRD1", "KLRF1"]),
        "FCGR3A Mono": (160, ["FCGR3A", "MS4A7", "LST1", "CDKN1C", "IFITM3"]),
        "DC": (40, ["FCER1A", "CST3", "CLEC10A", "CD1C"]),
        "Platelet": (18, ["PPBP", "PF4", "GP9", "TUBB1"]),
    }
    mito = ["MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "MT-ND4", "MT-ND5", "MT-ATP6", "MT-CYB", "MT-ND3", "MT-ATP8", "MT-ND6", "MT-ND4L"]
    ribo = [f"RP{x}{i}" for x in "LS" for i in range(3, 30)]
    house = ["MALAT1", "ACTB", "GAPDH", "B2M", "TMSB4X", "EEF1A1", "FTL", "FTH1", "HLA-B", "HLA-A", "PTPRC", "CD74", "HLA-DRA", "TPT1"]
    isg = ["ISG15", "IFI6", "MX1", "IFIT1", "IFIT3", "OAS1", "IFI44L", "RSAD2", "STAT1", "IRF7"]
    cc = ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "TYMS", "CDK1", "UBE2C", "BIRC5", "TPX2", "NUSAP1", "HMGB2", "RRM2", "FEN1", "GINS2", "MCM4", "MCM6", "CCNB2", "CENPF", "CKS2", "SMC4", "AURKB", "CDC20"]
    mk = []
    for _, (_, g) in types.items():
        mk += [x for x in g if x not in mk]
    genes = house + mk + mito + ribo + [g for g in cc if g not in mk] + isg
    genes += _background_genes(12000 - len(genes), set(genes), rng)
    G = len(genes); gi = {g: i for i, g in enumerate(genes)}
    base = np.exp(rng.normal(-2.2, 1.6, G))
    base[[gi[g] for g in house]] = rng.uniform(15, 60, len(house))
    base[[gi[g] for g in ribo]] = rng.uniform(3, 10, len(ribo))
    base[[gi[g] for g in mito]] = rng.uniform(1.5, 5, len(mito))
    base[[gi[g] for g in mk]] = 0.03
    rows, obs = [], []
    for t, (n, marks) in types.items():
        prog = base.copy()
        prog[rng.choice(np.arange(len(house) + len(mk) + 80, G), 80, replace=False)] *= rng.uniform(3, 8, 80)
        for g in marks:
            prog[gi[g]] = rng.uniform(2, 12) if g not in ("LYZ", "S100A8", "S100A9", "CD74", "CST3") else rng.uniform(20, 60)
        for c in range(n):
            donor = rng.integers(6)
            stim = donor >= 3
            kind = rng.random()
            damaged, empty = kind < 0.05, kind > 0.97
            sf = np.exp(rng.normal(0, 0.35)) * (0.15 if empty else 0.6 if damaged else 1.0) * (1.25 if "Mono" in t else 0.5 if t == "Platelet" else 1)
            mu = prog * sf * 0.55
            if damaged:
                mu[[gi[g] for g in mito]] *= 10
            mu[:200] *= 1 + 0.04 * (donor % 3)
            if stim:
                mu[[gi[g] for g in isg]] *= 12          # interferon response in donors 4-6
                mu[[gi[g] for g in isg]] += 2
            if rng.random() < 0.03 and t in ("CD4 T", "CD8 T", "NK"):
                mu[[gi[g] for g in cc]] += rng.uniform(0.5, 2, len(cc))
            x = rng.poisson(rng.gamma(2.0, mu / 2.0))
            rows.append(sp.csr_matrix(x.astype(np.float32)))
            obs.append((t, f"donor{donor+1}", "IFN-stimulated" if stim else "control"))
    # doublets
    X = sp.vstack(rows).tocsr()
    nd = 70
    a, b = rng.integers(0, X.shape[0], nd), rng.integers(0, X.shape[0], nd)
    X = sp.vstack([X, X[a] + X[b]]).tocsr()
    obs += [("doublet", obs[i][1], obs[i][2]) for i in a]
    bc = ["".join(rng.choice(list("ACGT"), 16)) + "-1" for _ in range(X.shape[0])]
    adata = ad.AnnData(X, obs=pd.DataFrame({"sample": [o[1] for o in obs], "condition": [o[2] for o in obs]}, index=bc), var=pd.DataFrame(index=genes))
    return adata
