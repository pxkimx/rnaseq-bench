# RNAseq Bench Web

Single-cell QC, clustering and marker genes, computed **in the browser tab**. Three static files and no
server: data is read by the page and never leaves the machine.

    index.html    the page
    app.js        boot, file reading, UMAP, plotting
    pipeline.py   the analysis, run by Pyodide

## Why it does less than the desktop app

Not a scoping choice — a hard limit, measured rather than assumed. In a real browser:

| | |
|---|---|
| numpy · scipy · pandas · scikit-learn · h5py · **igraph** | available |
| **scanpy / anndata** | blocked by `numcodecs`, `google-crc32c` (no wasm wheels) |
| **PyDESeq2** | will not install |
| **harmonypy** | now depends on `torch` |
| umap-learn · scVI · CellTypist · decoupler · liana | need numba or torch |

So the steps here call the same libraries scanpy calls underneath — `sklearn` for PCA and neighbours,
`igraph.community_leiden` for the clustering, a vectorised rank-sum for markers. UMAP is the exception
and runs in JavaScript (`umap-js`), fed the neighbour graph scikit-learn already computed.

**Differential expression, batch integration, annotation and gene sets are absent, not approximated.**
Substituting something weaker and presenting it as the real thing is the one outcome worth avoiding.

Practical to roughly 20–30k cells; neighbours are exact rather than approximate, which is where the
time goes.

## Run locally

    python3 -m http.server 8792     # from this folder, then open localhost:8792

## Deploy to Cloudflare Pages

Nothing to build — it is already static.

**Dashboard:** Workers & Pages → Create → Pages → *Upload assets*, drag this folder in. Done.

**CLI** (needs Node):

    npm install -g wrangler
    wrangler pages deploy . --project-name rnaseq-bench-web

**From Git:** point Pages at the repo, set the build output directory to `source/webapp` and leave the
build command empty.

Pyodide is fetched from jsdelivr on first load (~30 MB, then cached by the browser). No cross-origin
isolation headers are required.
