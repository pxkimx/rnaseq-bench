# RNAseq Bench changelog

## 3.2.2
Found while reproducing a published figure from GEO: GSE264154.

- **A cluster that is essentially one sample is now flagged.** In that dataset one donor's library was
  sequenced about a third as deeply as the others, and Harmony cannot merge a difference in how much was
  detected. It formed its own clusters covering 54% of the cells — half the UMAP was one library. The
  run now says which sample, how much of the data it is, and how its depth compares with the rest.
- **The adaptive upper gene cut-off could land above the busiest cell**, which is not a threshold but a
  number that removes nothing (it read 260,047 on a dataset whose deepest cell had 19,057 genes). It is
  now capped at the observed maximum, and the Parameters section says when that happened — routine in
  single-nucleus data, where the spread of gene counts is wide.

## 3.2.1
Everything here was found by running v3.2.0: on real GEO data, and on the packaged Mac app rather than
from source.

Fixes — the Mac app
- **Downloads that need the network now work.** The app's Python ships no CA certificates and does not
  read the keychain, so every `urllib` download failed with a certificate error: PROGENy and CollecTRI
  reported "needs internet" on machines that had internet, while gene sets (which use a different HTTP
  library) worked. All TLS clients now use certifi.
- The launcher installs the expert add-on packages separately from the ones the app needs to run, so a
  failure there can no longer stop the app from starting, and it cannot downgrade numpy/pandas/scanpy to
  satisfy an add-on.
- `liana` (cell–cell communication) is unavailable on this stack: releases new enough for the current
  anndata require pandas 2, and the rest of the analysis runs on pandas 3. The step now says so plainly
  instead of failing with a type error; nothing else is affected.
- build_mac.sh wrote the versioned zip one directory above the app.

Fixes found by running real GEO data (GSE113957)
- **Only the first GEO series matrix was read.** A series sequenced on two instruments ships one per
  platform; the smaller one sorted first, matched almost nothing, and every sample annotation was
  discarded. All series matrices are now merged.
- **With no metadata, samples were split into arbitrary halves called A and B** and compared as though
  that meant something. That split is now refused with instructions, and a name-inferred grouping is
  flagged for checking rather than reported as fact.
- **RefSeq gene IDs (NM_/NR_) are now mapped to symbols** (UCSC refGene). Without this, gene sets,
  PROGENy and CollecTRI matched nothing and reported themselves as skipped.
- Numeric annotation columns such as `Copies` were analysed as if they were samples.
- GEO characteristics are free text: ages like `8yr` or `2yr3mos` made the whole column non-numeric, so
  a continuous variable became a 77-level categorical; `Male`/`male` were separate levels. Both are now
  cleaned, and the merge is reported.
- The generated script now sets its own certificate bundle, and each optional step runs under a guard so
  one unavailable resource prints "[skipped] …" instead of ending the script.
- Sub-clustering reported the number of cells left after quality control as though it were the number
  selected.
- Changing only the compared level so that it matched the reference was not caught before re-analysis.

## 3.2
New — expert analyses
- **Bulk**: pre-ranked GSEA (Hallmark, GO BP, KEGG, Reactome); PROGENy pathway and CollecTRI
  transcription-factor activity, for the contrast and per sample; a genes-of-interest panel;
  covariate-adjusted PCA.
- **Single-cell**: static feature plots, UMAP split by condition, cluster→label and label→condition
  alluvials, CellTypist annotation, PROGENy activity per population, PAGA + diffusion pseudotime,
  LIANA ligand–receptor analysis, silhouette per resolution, sub-clustering of a previous job's
  clusters, and GSEA on the pseudobulk comparison.

New — the Parameters section
- Every analysis now carries a **Parameters** section: each choice the pipeline made, the value it
  used, where that value came from (derived from your data / chosen by you / built-in default /
  fixed) and a plain-language reason. Change anything and press **Re-analyze**.
- Only what you change is sent, so the rest still adapts to your data; clearing a field returns it
  to automatic.
- The section is also an appendix in the PDF report, and the assistant reads it before discussing
  or changing a parameter.

Fixes
- Choosing a covariate that describes the same split of the samples as the design factor crashed
  DESeq2 with "Singular matrix". Such columns are now detected before the fit and explained, are
  not offered as covariates, and are no longer suggested by the PC-association advice.
- The generated script covers the expert add-ons that a job actually ran, and its requirements line
  lists the packages they need.
- The PDF report header said TRANSCRIPT LENS.
- Pipeline-derived cell columns (qc_pass, leiden, …) are no longer offered as batch, sample or
  condition columns.

## 3.1.3
- No more "running" dialog: the server runs in the background; stop it from the sidebar ("Quit RNAseq Bench") or with "Stop RNAseq Bench.command".

## 3.1.2
- The .app declares arm64 (LSArchitecturePriority) and the launcher re-executes natively if macOS started it under Rosetta — this was why package installs from the app tried to compile Intel wheels.

## 3.1.1
- Launcher: rebuilds a broken package environment instead of failing; only uses a Python matching the Mac's chip (no Rosetta builds); the rename no longer moves the venv.

## 3.1
- Renamed to **RNAseq Bench** (was Transcript Lens). Existing analyses, packages and the API key are
  moved automatically to ~/Library/Application Support/RNAseqBench and ~/.rnaseq-bench on first launch.

## 3.0.1
- Python 3.13 supported (tested); 3.14 explicitly rejected with an explanation (no harmonypy wheel, scanpy import error).

## 3.0 — rebuilt and re-tested from a clean install
Fixes
- Single-cell runs crashed at "Scaling and PCA" with *"Both input arrays must be (arrays of) 3-dimensional
  vectors"* on current NumPy (2.5+): the PCA knee finder used np.cross on 2-D points. Rewritten.
- Scrublet was silently skipped on fresh installs (needs scikit-image) — now installed; every skipped step
  is logged with its traceback.
- LFC shrinkage silently did nothing whenever the alternative level sorted before the reference
  (e.g. HGPS vs WT); the reference level is now the baseline of the design in bulk, pseudobulk and scripts.
- featureCounts-style tables (Chr/Start/End/Strand/Length columns) were mis-detected as single-cell and
  then failed to parse; type detection now looks at sequencing depth and column names, not column count.
- Generated scripts: multi-sample 10x loading, per-sample sheet (sample_metadata.csv), samples without a
  factor value dropped, reference baseline — the bulk script now reproduces the app's numbers exactly.
GEO import
- RAW.tar handled recursively (per-sample .tar.gz, per-sample folders, junk files ignored); per-sample 10x
  triplets / .h5 / .csv.gz / HTSeq count files are stacked and named by their GSM prefix.
- Series-matrix annotations join onto stacked samples through the GSM prefix, so genotype/treatment
  columns are available for pseudobulk without extra files.
- .xlsx count tables, transposed tables, comment lines and odd encodings are read; empty droplets are
  dropped early when a raw (unfiltered) matrix is given.
- Download progress in MB; better pre-ticked files; certificate fallback; clear error messages.
App
- Launcher rewritten: detects an older running copy and replaces it, re-installs packages when
  requirements change, self-tests the install, writes a log. "Run self-test.command" added.
- Settings: "Test key & list models" — the model list comes from your key; retired models fall back.
- Failed runs show a "Technical details" panel with traceback + package versions.

## 2.3
- Failed analyses now show a 'Technical details' panel with the full traceback and package versions (Copy button).

## 2.2
- 10x matrix loading fixed: reads matrix/features|genes/barcodes directly (gzipped or not, any file names);
  multiple GSM-prefixed samples are stacked with the prefix stored as `sample`.
- Version number shown in the sidebar and in the app/zip name.

## 2.1
- GEO import: falls back to certifi root certificates (macOS python.org builds); shows the real network error.

## 2.0
- Standalone app: modern UI, GEO import, Claude assistant, code viewer, PDF reports, Mac app bundle.
