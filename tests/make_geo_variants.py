import sys, gzip, pathlib, tarfile, shutil, io, numpy as np, scipy.sparse as sp, scipy.io, h5py, pandas as pd
sys.path.insert(0, "/home/claude/rnaseq-bench")
from server.demo import simulate_pbmc
a = simulate_pbmc(); key="sample"
samples = sorted(a.obs[key].unique())
def write_h5(path, sub):
    Xc = sp.csc_matrix(sub.X.T)
    with h5py.File(path, "w") as f:
        g = f.create_group("matrix")
        g.create_dataset("data", data=Xc.data.astype(np.int32)); g.create_dataset("indices", data=Xc.indices.astype(np.int64))
        g.create_dataset("indptr", data=Xc.indptr.astype(np.int64)); g.create_dataset("shape", data=np.array(Xc.shape, dtype=np.int32))
        g.create_dataset("barcodes", data=np.array([b.encode() for b in sub.obs_names]))
        ft = g.create_group("features")
        ft.create_dataset("id", data=np.array([f"ENSG{i:011d}".encode() for i in range(sub.n_vars)]))
        ft.create_dataset("name", data=np.array([v.encode() for v in sub.var_names]))
        ft.create_dataset("feature_type", data=np.array([b"Gene Expression"] * sub.n_vars))
        ft.create_dataset("genome", data=np.array([b"GRCh38"] * sub.n_vars)); ft.create_dataset("_all_tag_keys", data=np.array([b"genome"]))
sm = pathlib.Path("data10x/GSE999999_series_matrix.txt.gz")
# A: RAW.tar of per-sample .tar.gz each holding filtered_feature_bc_matrix/{matrix,features,barcodes}
d = pathlib.Path("data_nested"); shutil.rmtree(d, ignore_errors=True); d.mkdir(); tmp = d / "tmp"; tmp.mkdir()
inner = []
for i, s in enumerate(samples):
    sub = a[a.obs[key] == s]; gsm = f"GSM{9000001+i}"
    fold = tmp / f"{gsm}_{s}" / "filtered_feature_bc_matrix"; fold.mkdir(parents=True)
    with gzip.open(fold / "matrix.mtx.gz", "wb") as f: scipy.io.mmwrite(f, sp.coo_matrix(sub.X.T))
    with gzip.open(fold / "features.tsv.gz", "wt") as f: f.write("".join(f"ENSG{j:011d}\t{g}\tGene Expression\n" for j, g in enumerate(sub.var_names)))
    with gzip.open(fold / "barcodes.tsv.gz", "wt") as f: f.write("".join(f"{b}\n" for b in sub.obs_names))
    tgz = tmp / f"{gsm}_{s}.tar.gz"
    with tarfile.open(tgz, "w:gz") as t: t.add(fold.parent, arcname=f"{gsm}_{s}")
    inner.append(tgz)
(tmp / "GSM9000001_donor1_peaks.bed.gz").write_bytes(gzip.compress(b"chr1\t1\t2\n"))
(tmp / "GSE999999_README.txt").write_text("hello\n")
with tarfile.open(d / "GSE999999_RAW.tar", "w") as t:
    for f in inner + [tmp / "GSM9000001_donor1_peaks.bed.gz", tmp / "GSE999999_README.txt"]: t.add(f, arcname=f.name)
shutil.rmtree(tmp); shutil.copy(sm, d)
# B: per-sample .h5 in RAW.tar
d = pathlib.Path("data_multih5"); shutil.rmtree(d, ignore_errors=True); d.mkdir(); tmp = d / "tmp"; tmp.mkdir()
for i, s in enumerate(samples):
    write_h5(tmp / f"GSM{9000001+i}_{s}_filtered_feature_bc_matrix.h5", a[a.obs[key] == s])
with tarfile.open(d / "GSE999999_RAW.tar", "w") as t:
    for f in sorted(tmp.iterdir()): t.add(f, arcname=f.name)
shutil.rmtree(tmp); shutil.copy(sm, d)
# C: per-sample dense csv.gz
d = pathlib.Path("data_multicsv"); shutil.rmtree(d, ignore_errors=True); d.mkdir()
for i, s in enumerate(samples):
    sub = a[a.obs[key] == s]
    pd.DataFrame(sub.X.T.toarray().astype(int), index=sub.var_names, columns=sub.obs_names).to_csv(d / f"GSM{9000001+i}_{s}_counts.csv.gz")
shutil.copy(sm, d)
print("ok")
