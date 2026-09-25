# RNAseq Bench — notes for Claude Code

Local app for bulk + single-cell RNA-seq analysis (FastAPI + Scanpy/PyDESeq2, vanilla-JS UI, in-app Claude assistant).
Owner: Paul (biologist, learning the analyses — explain the *why* in plain language in every figure's "how"/"yours" text).

## Layout
- `server/app.py` HTTP API · `server/sc_pipeline.py` single-cell · `server/bulk_pipeline.py` bulk · `server/pseudobulk.py`
  sample-level tests · `server/sc_extra.py` / `server/bulk_extra.py` expert add-ons (v3.2) · `server/resources.py`
  cached downloads (gene sets, PROGENy, CollecTRI, CellTypist models) · `server/io_utils.py` readers · `server/geo.py`
  GEO + Ensembl · `server/agent.py` assistant tools · `server/codegen.py` reproducible script · `server/report.py` PDF ·
  `server/dossier.py` gene dossier (Biopython Entrez/UniProt/KEGG/alphafold_db; `/api/dossier/{gene}`, agent tool
  `gene_dossier`, UI `openDossier()`) · `server/selftest.py` · `web/index.html` the whole UI (single file) · `launch.sh` Mac launcher · `tests/` harness.
- Job output: `TL_HOME/jobs/<id>/` (result.json drives the UI; figures/*.png; CSVs). `TL_HOME` = ~/Library/Application Support/RNAseqBench.
- Every figure goes through `job.figure(id, title, fig, how=, yours=)`; sections via `job.section`. Keep that contract —
  the UI, PDF report and agent all read result.json.

## webapp/ — the browser-only version
`source/webapp/` is a separate, static app (index.html + app.js + pipeline.py) that runs QC → HVG → PCA →
neighbours → Leiden → markers inside a browser tab via Pyodide. It shares no code with `server/` on purpose:
scanpy, anndata, PyDESeq2 and harmonypy cannot run in WebAssembly (numcodecs, google-crc32c, numba, torch
have no wasm wheels — verified, not assumed), so it calls sklearn/igraph directly and does UMAP in JS.
Never add an approximation of DESeq2, Harmony or CellTypist there; those steps are absent by design and the
page says so. Deploys to Cloudflare Pages as static assets.

## Run / build / test
- `./dev.sh` — server from source with auto-reload (uses the app's venv). `./build_mac.sh` — refresh `../RNAseq Bench.app` in place + zip.
- `python -m server.selftest` — 2-minute simulated sc + bulk run. `tests/README.md` — synthetic GEO-style datasets and `tests/run_job.py`.
- Offline resources for tests: `RB_RESOURCE_DIR=tests/res_offline` (make with `tests/make_resources.py`).
- Bump `VERSION` + `CHANGELOG.md` for every change that ships; the launcher uses VERSION to replace a running older server.

## Conventions / gotchas
- No `np.cross` on 2-D vectors (NumPy 2.5 errors). pandas 3 is in use (string dtype). Python 3.10–3.13; 3.14 unsupported (harmonypy, scanpy).
- Anything needing internet (GEO, Enrichr, gene sets, PROGENy, CellTypist) must degrade to an `[info]` flag, never fail the run.
  That rule covers the generated script too — its optional steps run under `step()` in `codegen.py`.
- TLS: `server/__init__.py` points SSL_CERT_FILE/REQUESTS_CA_BUNDLE at certifi. The Mac app's python.org build
  has no CA bundle, so without it every `urllib` download (decoupler/omnipath, UCSC) fails with
  CERTIFICATE_VERIFY_FAILED while `requests`-based ones (GSEApy) work — a confusing half-offline state.
- GEO metadata is free text and one series matrix per platform: read them all (`io_utils.read_series_matrices`),
  clean level capitalisation and units (`unify_level_case`, `numeric_with_units`) before deciding what is
  continuous. Gene IDs may be Ensembl **or RefSeq** (`map_gene_ids`); unmapped IDs make every gene-set step
  silently find nothing.
- Never report an arbitrary grouping as a result: `infer_groups` returns how it guessed, and the A/B fallback
  is refused rather than analysed.
- Non-fatal errors: wrap in try/except and call `log_exc("where")` so the traceback lands in the server log.
- DESeq2: reference level must be first in the Categorical (else lfc_shrink silently does nothing).
- Cells are not replicates: single-cell DE between conditions is pseudobulk per sample.
- Gene dossier: every source runs under `guarded` + a 25 s `DEADLINE` (Entrez and alphafold_db open URLs with no
  timeout; Entrez retries 3× 15 s by default — set to 2× 2 s). Only complete dossiers are cached
  (`~/.rnaseq-bench/dossier/`). KEGG gene ids for hsa/mmu/rno are NCBI Gene ids (Biopython's `kegg_conv`
  rejects `ncbi-geneid` sources). `python tests/test_dossier.py` checks the offline/hung paths without network.
- Sister apps: MassSpec Bench :8766, Structure Bench :8767 (dossier links to `#uniprot=` / `#gene=&species=`),
  Clone Bench :8768.
- Mac app runs `launch.sh` (arm64 re-exec guard, venv rebuild guard, requirements hash → reinstall).

## v3.2 status (shipped 3.2.0; 3.2.1 is the real-data / packaged-app fix round)
Shipped and tested on the synthetic datasets in `tests/`: bulk GSEA prerank, PROGENy/CollecTRI activity
(contrast + per-sample), genes-of-interest, covariate-adjusted PCA; sc feature plots, UMAP split by condition,
alluvials, PROGENy per population, PAGA + DPT, LIANA, silhouette per resolution, sub-clustering
(`subset`, `subset_key`, `subset_from_job`), GSEA on pseudobulk. CellTypist is wired but its model needs
internet, so it is still untested here.

**Parameters section** (`server/params_spec.py`): each pipeline ends its run by writing
`result["params_panel"]` — groups of fields, each with `key`, `label`, `value`, `why`, `source`
(auto / user / default / data) and a control `type`. `web/index.html` renders it generically as the
"Parameters" report section, and `_panel()` in `sc_pipeline.py` / `bulk_pipeline.py` builds it.
Rules when you add a parameter:
- Every editable field's `key` must be a parameter the pipeline's `run(params=...)` accepts, or Re-analyze
  drops it silently. `type="fixed"` and keys starting with `_` are display-only.
- The `why` is generated at analysis time and must say what actually happened, not what was requested —
  use `params_spec.skipped()` for add-ons that were skipped, and name the derivation for `auto` values.
- The UI sends only the fields you changed; `null` clears one back to automatic (`/api/rerun` merges over
  the previous job's params). Never send the whole panel back, or auto-derived thresholds become pinned.
- Column choices come from `params_spec.real_obs()` (drops qc_pass, leiden, …); bulk covariate choices drop
  columns confounded with the design factor.
- `result["params"]` still carries the flat values for the agent, codegen and re-runs — keep both updated.

TODO: CellTypist on a real (online) run; RNA velocity (scVelo) is the one publication-checklist item not
covered, and needs spliced/unspliced counts most datasets do not ship.
