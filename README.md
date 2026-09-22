# RNAseq Bench

A local app for bulk and single-cell RNA-seq: upload counts (or pull them straight from GEO), get the
standard Scanpy / DESeq2 analyses with publication-style figures, a plain-language reading of every
plot, a PDF report, the exact code that produced it — and a built-in Claude assistant that can explain,
re-run and fix analyses for you.

## Start

**Mac app** — unzip `RNAseq Bench vX.Y.zip`, run `xattr -cr` on the app once (see README-FIRST.txt),
double-click. First launch installs packages into `~/Library/Application Support/RNAseqBench`
(3–5 min); the launcher replaces an older running copy, re-installs when requirements change, and
writes `rnaseq-bench.log` there. `Run self-test.command` runs a 2-minute simulated analysis.

**From source (macOS / Linux)**
```bash
cd rnaseq-bench
./start.sh          # installs into .venv on first run, then opens http://localhost:8765
python -m server.selftest   # optional check of the installation
```
Needs Python 3.10–3.13 (not 3.14 yet — see below). To use the assistant, open **Settings & API key** in the sidebar, paste an
Anthropic API key (stored only in `~/.rnaseq-bench/config.json`) and click *Test key*.

## GEO import — what to tick
- `*_series_matrix.txt.gz` — always: sample annotations (genotype, treatment, age…) are joined onto
  your samples by name or GSM id and drive the comparisons.
- a counts/matrix supplementary file when there is one (`*counts*.csv.gz`, `*.h5ad`, `*.h5`, `*.xlsx`).
- otherwise `GSExxxxxx_RAW.tar` — per-sample 10x triplets, `.h5`, count tables or HTSeq files are
  unpacked (recursively), stacked, and named by their GSM prefix (`obs["sample"]`).
Bulk vs single-cell is detected from the data (depth, barcodes, column names); override it with the
data-type switch above the upload box.

## What it does

**Inputs** — `.h5ad`, 10x `.h5` / `matrix.mtx` folders, `.loom`, dense CSV/TSV (any size — streamed to
sparse), optional sample or cell metadata CSV, a **GEO series matrix** (sample annotations are matched to
your columns automatically), `.zip` / `.tar.gz`, or an accession typed into **Import from GEO**.
Ensembl IDs are mapped to symbols (human/mouse, cached after first use).

**Bulk (PyDESeq2)** — technical-replicate detection and summing, sex inference from XIST/chrY,
categorical *and continuous* covariates (age, RIN…), design formula you control, DESeq2 with LFC
shrinkage, VST → PCA, correlation, PC–metadata association (batch detection), dispersion, volcano, MA,
p-value histogram, heatmap, per-gene counts, Enrichr pathways, pre-ranked **GSEA** (Hallmark, GO BP,
KEGG, Reactome), **PROGENy** pathway and **CollecTRI** transcription-factor activity for the contrast and
per sample, and a panel for your own genes of interest.

**Single-cell (Scanpy)** — MAD-based QC, Scrublet (auto-skipped for plate-based data), HVGs, PCA,
Harmony, UMAP, Leiden with resolution sweep, Wilcoxon markers, marker-panel annotation (immune, tissue and
vascular panels), interactive UMAP with gene lookup — and when cells carry a sample and a condition,
**composition tests** and **pseudobulk DESeq2** (all cells and per cell type), so cells are never
mistaken for replicates. Also: feature plots, UMAP split by condition, cluster→label and label→condition
alluvials, **CellTypist** annotation, **PROGENy** activity per population, **PAGA + diffusion pseudotime**,
**LIANA** ligand–receptor analysis, silhouette per resolution, GSEA on the pseudobulk comparison, and
sub-clustering of any cluster from a previous run.

**Parameters** — every analysis ends with a Parameters section listing each choice it made: the value
used, where it came from (derived from your data, chosen by you, a built-in default, or not adjustable)
and a plain-language reason. Change anything and press **Re-analyze**; only what you change is overridden,
so everything else still adapts to your data, and clearing a field returns it to automatic. The same
account appears as an appendix in the PDF.

**Assistant** — floating *Ask Claude* button. It can read results, list and read files, re-run with new
parameters, write a corrected metadata file, run Python inside the job folder (custom plots, extra
statistics), and hand you the reproducible script. Tool calls are shown in the chat as expandable cards.

**View code** — every analysis has a `reproduce_analysis.py` written with plain Scanpy / PyDESeq2 calls
and the parameters used; the same viewer shows the pipeline source itself.

## Web version
`source/webapp/` is a browser-only build — QC, clustering and marker genes computed in the tab with nothing
uploaded, deployable to Cloudflare Pages as three static files. It deliberately omits differential
expression, integration, annotation and gene sets, which cannot run in WebAssembly. See `webapp/README.md`.

## Outputs
Per job in `jobs/<id>/`: `report.pdf`, `figures/*.png` (200 dpi), `deseq2_results.csv` /
`pseudobulk_*.csv`, `composition_test.csv`, `analyzed.h5ad`, `markers.csv`, `cell_metadata.csv`,
`sample_metadata.csv`, `metadata_used.csv`, `counts_used.csv`. See CHANGELOG.md for what changed per version.

## Layout
```
server/app.py            HTTP API            server/agent.py       Claude assistant + tools
server/bulk_pipeline.py  bulk analysis        server/codegen.py     reproducible scripts
server/sc_pipeline.py    single-cell          server/geo.py         GEO + Ensembl helpers
server/pseudobulk.py     sample-level tests   server/io_utils.py    readers, replicates, sex, ID mapping
server/report.py         PDF                  web/index.html        UI
```

## Caveats
Labels from the marker panel are hypotheses. Complex designs (interactions, time courses) need a custom
model — ask the assistant to write it. Everything runs on your machine; only the assistant talks to
Anthropic's API, and only with what you type plus the tool results it requests.
