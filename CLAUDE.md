# RNAseq Bench — notes for Claude Code

Local app for bulk + single-cell RNA-seq analysis (FastAPI + Scanpy/PyDESeq2, vanilla-JS UI, in-app Claude assistant).
Owner: Paul (biologist, learning the analyses — explain the *why* in plain language in every figure's "how"/"yours" text).

## Layout
- `server/app.py` HTTP API · `server/sc_pipeline.py` single-cell · `server/bulk_pipeline.py` bulk · `server/pseudobulk.py`
  sample-level tests · `server/sc_extra.py` / `server/bulk_extra.py` expert add-ons (v3.2) · `server/resources.py`
  cached downloads (gene sets, PROGENy, CollecTRI, CellTypist models) · `server/io_utils.py` readers · `server/geo.py`
  GEO + Ensembl · `server/agent.py` assistant tools · `server/codegen.py` reproducible script · `server/report.py` PDF ·
  `server/selftest.py` · `web/index.html` the whole UI (single file) · `launch.sh` Mac launcher · `tests/` harness.
- Job output: `TL_HOME/jobs/<id>/` (result.json drives the UI; figures/*.png; CSVs). `TL_HOME` = ~/Library/Application Support/RNAseqBench.
- Every figure goes through `job.figure(id, title, fig, how=, yours=)`; sections via `job.section`. Keep that contract —
  the UI, PDF report and agent all read result.json.

## Run / build / test
- `./dev.sh` — server from source with auto-reload (uses the app's venv). `./build_mac.sh` — refresh `../RNAseq Bench.app` in place + zip.
- `python -m server.selftest` — 2-minute simulated sc + bulk run. `tests/README.md` — synthetic GEO-style datasets and `tests/run_job.py`.
- Offline resources for tests: `RB_RESOURCE_DIR=tests/res_offline` (make with `tests/make_resources.py`).
- Bump `VERSION` + `CHANGELOG.md` for every change that ships; the launcher uses VERSION to replace a running older server.

## Conventions / gotchas
- No `np.cross` on 2-D vectors (NumPy 2.5 errors). pandas 3 is in use (string dtype). Python 3.10–3.13; 3.14 unsupported (harmonypy, scanpy).
- Anything needing internet (GEO, Enrichr, gene sets, PROGENy, CellTypist) must degrade to an `[info]` flag, never fail the run.
- Non-fatal errors: wrap in try/except and call `log_exc("where")` so the traceback lands in the server log.
- DESeq2: reference level must be first in the Categorical (else lfc_shrink silently does nothing).
- Cells are not replicates: single-cell DE between conditions is pseudobulk per sample.
- Mac app runs `launch.sh` (arm64 re-exec guard, venv rebuild guard, requirements hash → reinstall).

## v3.2 status (in progress)
Done & tested on simulated data: bulk GSEA prerank, PROGENy/CollecTRI activity (contrast + per-sample), genes-of-interest,
covariate-adjusted PCA; sc feature plots, UMAP split by condition, alluvials, CellTypist (needs internet; untested here),
PROGENy per population, PAGA + DPT, LIANA cell–cell communication, silhouette per resolution, sub-clustering params
(`subset`, `subset_key`, `subset_from_job`), GSEA on pseudobulk.
TODO: (1) **Parameters panel** in the UI — list every parameter with the value used and *why* (auto-derived thresholds
explained), editable, one "Re-analyze" button (extend the tune drawer; sc params in `result.params`, bulk in `result.params`);
(2) expose new params in the tune drawer + agent `rerun_analysis` description (`genes`, `celltypist_model`, `root_cluster`,
`subset`, toggles `gsea/activity/celltypist/pathways/trajectory/ccc`); (3) codegen snippets for the add-ons; (4) PDF report
check with the new sections; (5) requirements: add celltypist, decoupler, liana; (6) run the full tests matrix, then ship.
