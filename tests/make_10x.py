"""Write simulated PBMC as 6 GSM-prefixed 10x v3 triplets + a GEO-style series matrix."""
import sys, gzip, pathlib, numpy as np, scipy.io, scipy.sparse as sp
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from server.demo import simulate_pbmc
out = pathlib.Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
a = simulate_pbmc()
print(a, a.obs.columns.tolist())
samples = sorted(a.obs["donor"].unique()) if "donor" in a.obs else sorted(a.obs["sample"].unique())
key = "donor" if "donor" in a.obs else "sample"
gsm_rows = []
for i, s in enumerate(samples):
    sub = a[a.obs[key] == s]
    gsm = f"GSM{9000001+i}"
    cond = sub.obs["condition"].iloc[0]
    X = sp.coo_matrix(sub.X.T)
    with gzip.open(out / f"{gsm}_{s}_matrix.mtx.gz", "wb") as f: scipy.io.mmwrite(f, X)
    with gzip.open(out / f"{gsm}_{s}_features.tsv.gz", "wt") as f:
        f.write("".join(f"ENSG{j:011d}\t{g}\tGene Expression\n" for j, g in enumerate(sub.var_names)))
    with gzip.open(out / f"{gsm}_{s}_barcodes.tsv.gz", "wt") as f:
        f.write("".join(f"{b.split('-')[0] if '-' in b else b}-1\n" for b in sub.obs_names))
    gsm_rows.append((gsm, s, cond))
# series matrix
lines = ['!Series_title\t"Simulated PBMC IFN"', '!Series_geo_accession\t"GSE999999"',
         "!Sample_title\t" + "\t".join(f'"{s} PBMC"' for _, s, _ in gsm_rows),
         "!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g, _, _ in gsm_rows),
         "!Sample_organism_ch1\t" + "\t".join('"Homo sapiens"' for _ in gsm_rows),
         "!Sample_characteristics_ch1\t" + "\t".join(f'"donor: {s}"' for _, s, _ in gsm_rows),
         "!Sample_characteristics_ch1\t" + "\t".join(f'"treatment: {c}"' for _, _, c in gsm_rows),
         "!Sample_supplementary_file_1\t" + "\t".join(f'"ftp://x/{g}_{s}_matrix.mtx.gz"' for g, s, _ in gsm_rows),
         "!series_matrix_table_begin", '"ID_REF"\t' + "\t".join(f'"{g}"' for g, _, _ in gsm_rows), "!series_matrix_table_end"]
with gzip.open(out / "GSE999999_series_matrix.txt.gz", "wt") as f: f.write("\n".join(lines) + "\n")
print("wrote", sorted(p.name for p in out.iterdir()))
