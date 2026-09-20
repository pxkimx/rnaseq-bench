import sys, gzip, pathlib, tarfile, shutil, numpy as np, scipy.sparse as sp, h5py, pandas as pd
sys.path.insert(0, "/home/claude/rnaseq-bench")
from server.demo import simulate_pbmc
a = simulate_pbmc()
# 1. RAW.tar of the 10x triplets + series matrix outside
d = pathlib.Path("data_rawtar"); shutil.rmtree(d, ignore_errors=True); d.mkdir()
with tarfile.open(d / "GSE999999_RAW.tar", "w") as tf:
    for p in sorted(pathlib.Path("data10x").iterdir()):
        if "series" not in p.name: tf.add(p, arcname=p.name)
shutil.copy("data10x/GSE999999_series_matrix.txt.gz", d)
# 2. dense csv.gz cells as columns, Ensembl-like? use symbols; plus metadata csv
d = pathlib.Path("data_dense"); shutil.rmtree(d, ignore_errors=True); d.mkdir()
X = pd.DataFrame(a.X.T.toarray().astype(int), index=a.var_names, columns=[f"{s}_{b}" for s, b in zip(a.obs['sample'], a.obs_names)])
X.to_csv(d / "GSE999999_counts.csv.gz")
a.obs.rename(columns={"condition": "treatment"}).set_axis(X.columns).to_csv(d / "GSE999999_cell_metadata.csv")
# 3. 10x .h5 (v3 format)
d = pathlib.Path("data_h5"); shutil.rmtree(d, ignore_errors=True); d.mkdir()
Xc = sp.csc_matrix(a.X.T)  # genes x cells, CSC as cellranger writes
with h5py.File(d / "GSM9000001_filtered_feature_bc_matrix.h5", "w") as f:
    g = f.create_group("matrix")
    g.create_dataset("data", data=Xc.data.astype(np.int32)); g.create_dataset("indices", data=Xc.indices.astype(np.int64))
    g.create_dataset("indptr", data=Xc.indptr.astype(np.int64)); g.create_dataset("shape", data=np.array(Xc.shape, dtype=np.int32))
    g.create_dataset("barcodes", data=np.array([b.encode() for b in a.obs_names]))
    ft = g.create_group("features")
    ft.create_dataset("id", data=np.array([f"ENSG{i:011d}".encode() for i in range(a.n_vars)]))
    ft.create_dataset("name", data=np.array([v.encode() for v in a.var_names]))
    ft.create_dataset("feature_type", data=np.array([b"Gene Expression"] * a.n_vars))
    ft.create_dataset("genome", data=np.array([b"GRCh38"] * a.n_vars))
    ft.create_dataset("_all_tag_keys", data=np.array([b"genome"]))
# 4. h5ad with normalized X but counts layer, ensembl var names
d = pathlib.Path("data_h5ad"); shutil.rmtree(d, ignore_errors=True); d.mkdir()
b = a.copy(); b.layers["counts"] = b.X.copy()
import scanpy as sc; sc.pp.normalize_total(b); sc.pp.log1p(b)
b.var["symbol"] = b.var_names; b.var_names = [f"ENSG{i:011d}" for i in range(b.n_vars)]
b.write_h5ad(d / "GSE999999_processed.h5ad")
print("ok")
