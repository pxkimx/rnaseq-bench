import pathlib, tarfile, gzip, shutil, numpy as np, pandas as pd
src = pathlib.Path("data_bulk"); c = pd.read_csv(src / "GSE96870_counts.csv", index_col=0); m = pd.read_csv(src / "GSE96870_metadata.csv", index_col=0)
# A: RAW.tar with per-sample HTSeq 2-col files (gene \t count, __no_feature rows), series matrix
d = pathlib.Path("data_htseq"); shutil.rmtree(d, ignore_errors=True); d.mkdir(); tmp = d / "tmp"; tmp.mkdir()
gsms = {s: f"GSM{2545000+i}" for i, s in enumerate(c.columns)}
for s in c.columns:
    ser = c[s].copy(); ser.loc["__no_feature"] = 12345; ser.loc["__ambiguous"] = 10
    with gzip.open(tmp / f"{gsms[s]}_{s}.htseq.counts.txt.gz", "wt") as f: f.write("".join(f"{g}\t{int(v)}\n" for g, v in ser.items()))
with tarfile.open(d / "GSE96870_RAW.tar", "w") as t:
    for f in sorted(tmp.iterdir()): t.add(f, arcname=f.name)
shutil.rmtree(tmp)
lines = ['!Series_title\t"x"', "!Sample_title\t" + "\t".join(f'"{s}"' for s in c.columns), "!Sample_geo_accession\t" + "\t".join(f'"{gsms[s]}"' for s in c.columns),
         "!Sample_characteristics_ch1\t" + "\t".join(f'"infection: {m.loc[s, "condition"]}"' for s in c.columns),
         "!Sample_characteristics_ch1\t" + "\t".join(f'"Sex: {m.loc[s, "sex"]}"' for s in c.columns),
         "!Sample_characteristics_ch1\t" + "\t".join(f'"time: {m.loc[s, "time"]}"' for s in c.columns),
         "!series_matrix_table_begin", '"ID_REF"\t' + "\t".join(f'"{gsms[s]}"' for s in c.columns), "!series_matrix_table_end"]
with gzip.open(d / "GSE96870_series_matrix.txt.gz", "wt") as f: f.write("\n".join(lines) + "\n")
# B: xlsx supplementary with an annotation column + a README
d = pathlib.Path("data_xlsx"); shutil.rmtree(d, ignore_errors=True); d.mkdir()
x = c.copy(); x.insert(0, "gene_name", x.index.str.upper()); x.insert(1, "Length", 1000)
x.to_excel(d / "GSE96870_raw_counts.xlsx"); m.to_csv(d / "GSE96870_sample_info.csv"); (d / "README.txt").write_text("Counts from featureCounts\n")
print("ok")
