/* RNAseq Bench Web — glue between the page and the Python running in this tab.
   Pyodide does the numerics; UMAP is the one step done in JavaScript, because umap-learn needs numba
   and numba has no WebAssembly build. It is fed the neighbour graph scikit-learn already computed, so
   the layout comes from the same neighbourhoods the clustering used. */
const $ = s => document.querySelector(s), $$ = s => [...document.querySelectorAll(s)];
const PALETTE = ["#00e0cf","#ffb020","#ff5c8a","#4ea3ff","#b27bff","#2ee6a8","#ff7ae0","#ffd23f","#38bdf8",
  "#a3e635","#fb923c","#c084fc","#22d3ee","#f43f5e","#84cc16","#e879f9","#facc15","#60a5fa"];
const S = {files: [], py: null, emb: null, hidden: new Set(), mode: "cluster", vals: null};
const toast = m => {const t = document.createElement("div"); t.className = "toast"; t.textContent = m;
  document.body.append(t); setTimeout(() => t.remove(), 3200)};

/* ---------------------------------------------------------------- boot */
(async () => {
  const bar = $("#bootBar"), msg = $("#bootMsg");
  const tick = (p, m) => {bar.style.width = p + "%"; if (m) msg.textContent = m};
  try {
    tick(8, "Starting Python…");
    const py = await loadPyodide();
    tick(30, "Loading numpy, scipy, scikit-learn, igraph, h5py…");
    await py.loadPackage(["numpy", "scipy", "pandas", "scikit-learn", "igraph", "h5py"]);
    tick(85, "Loading the pipeline…");
    py.FS.writeFile("/pipeline.py", await (await fetch("pipeline.py")).text());
    py.runPython("import sys; sys.path.insert(0,'/'); import pipeline");
    S.py = py;
    tick(100, "Ready.");
    $("#bootCard").classList.add("hidden");
    $("#loadCard").classList.remove("hidden");
  } catch (e) {
    msg.innerHTML = `Could not start: ${String(e).slice(0, 300)}<br><small>A browser with WebAssembly is required.</small>`;
  }
})();

/* ---------------------------------------------------------------- files */
const zone = $("#zone");
zone.onclick = () => $("#file").click();
zone.ondragover = e => {e.preventDefault(); zone.classList.add("over")};
zone.ondragleave = () => zone.classList.remove("over");
zone.ondrop = e => {e.preventDefault(); zone.classList.remove("over"); addFiles(e.dataTransfer.files)};
$("#file").onchange = e => addFiles(e.target.files);

function addFiles(list) {
  for (const f of list) if (!S.files.some(x => x.name === f.name)) S.files.push(f);
  $("#fileList").innerHTML = S.files.map(f =>
    `<span class="tile" style="padding:5px 9px"><b style="font-size:12.5px">${f.name}</b>
     <span>${(f.size / 1e6).toFixed(1)} MB</span></span>`).join("");
  $("#loadBtn").disabled = !S.files.length;
}

/* .gz is the norm for anything from GEO, and browsers can now inflate it natively */
async function bytesOf(file) {
  const buf = await file.arrayBuffer();
  if (!file.name.toLowerCase().endsWith(".gz")) return new Uint8Array(buf);
  const ds = new DecompressionStream("gzip");
  const out = new Response(new Blob([buf]).stream().pipeThrough(ds));
  return new Uint8Array(await out.arrayBuffer());
}
const textOf = async f => new TextDecoder().decode(await bytesOf(f));
const find = re => S.files.find(f => re.test(f.name.toLowerCase()));

$("#loadBtn").onclick = async () => {
  $("#loadBtn").disabled = true;
  try {
    const py = S.py, g = py.globals;
    const mtx = find(/\.mtx(\.gz)?$/), feat = find(/(features|genes)\.tsv(\.gz)?$/), bc = find(/barcodes\.tsv(\.gz)?$/);
    const h5ad = find(/\.h5ad$/);
    let qc;
    if (mtx && feat && bc) {
      g.set("_m", await bytesOf(mtx)); g.set("_f", await textOf(feat)); g.set("_b", await textOf(bc));
      qc = py.runPython("import pipeline; pipeline.load_mtx(_m.to_py(), _f, _b)").toJs({dict_converter: Object.fromEntries});
    } else if (h5ad) {
      g.set("_h", await bytesOf(h5ad));
      qc = py.runPython("import pipeline; pipeline.load_h5ad(_h.to_py())").toJs({dict_converter: Object.fromEntries});
    } else {
      const tbl = find(/\.(csv|tsv|txt)(\.gz)?$/);
      if (!tbl) throw new Error("Need a 10x triplet, an .h5ad, or a counts table.");
      g.set("_t", await textOf(tbl));
      qc = py.runPython("import pipeline; pipeline.load_dense(_t)").toJs({dict_converter: Object.fromEntries});
    }
    showQC(qc);
  } catch (e) {
    toast("Could not read that: " + String(e.message || e).slice(0, 120));
    $("#loadBtn").disabled = false;
  }
};

$("#demoBtn").onclick = async () => {
  // a small synthetic set, so the pipeline can be tried without finding a file first
  $("#demoBtn").disabled = true;
  const py = S.py;
  const qc = py.runPython(`
import numpy as np, scipy.sparse as sp, pipeline
rng = np.random.default_rng(0)
n, g, k = 1800, 1200, 5
lab = rng.integers(0, k, n)
base = rng.gamma(0.4, 1.2, (k, g))
for c in range(k):
    base[c, c*60:(c+1)*60] *= 14                          # each group gets its own marker block
# mitochondrial genes live at the end, deliberately clear of every marker block: when they overlapped
# one, that whole population read as dying and the QC step removed it
base[:, -12:] = rng.gamma(0.6, 1.0, (k, 12))
M = rng.poisson(base[lab] * rng.uniform(0.5, 1.8, (n, 1)))
names = np.array([f"GENE{i}" for i in range(g)], dtype=object)
names[-12:] = np.array([f"MT-{i}" for i in range(12)], dtype=object)
pipeline._finish_load(sp.csr_matrix(M.astype(np.float32)), names,
                      np.array([f"cell{i}" for i in range(n)], dtype=object))
`).toJs({dict_converter: Object.fromEntries});
  showQC(qc);
};

/* ---------------------------------------------------------------- QC panel */
function histSvg(el, counts, edges, cut) {
  const w = 300, h = 90, max = Math.max(...counts, 1);
  const bw = w / counts.length;
  const bars = counts.map((c, i) =>
    `<rect x="${(i * bw).toFixed(1)}" y="${(h - (c / max) * (h - 8)).toFixed(1)}" width="${(bw - .6).toFixed(1)}"
      height="${((c / max) * (h - 8)).toFixed(1)}" fill="var(--accent)" opacity=".75"/>`).join("");
  const lines = (cut || []).filter(v => v != null).map(v => {
    const t = (v - edges[0]) / (edges[edges.length - 1] - edges[0]);
    return t < 0 || t > 1 ? "" : `<line x1="${(t * w).toFixed(1)}" y1="0" x2="${(t * w).toFixed(1)}" y2="${h}"
      stroke="var(--up)" stroke-dasharray="3 3"/>`}).join("");
  el.setAttribute("viewBox", `0 0 ${w} ${h}`); el.innerHTML = bars + lines;
}

function showQC(qc) {
  S.qc = qc;
  $("#loadCard").classList.add("hidden");
  $("#qcCard").classList.remove("hidden");
  $("#qcSub").textContent = `${qc.n_cells.toLocaleString()} cells × ${qc.n_genes_total.toLocaleString()} genes read.`;
  $("#qcTiles").innerHTML = [
    [qc.n_cells.toLocaleString(), "cells"], [qc.n_genes_total.toLocaleString(), "genes"],
    [qc.median_genes.toLocaleString(), "median genes / cell"], [qc.median_umi.toLocaleString(), "median counts / cell"],
  ].map(([v, l]) => `<div class="tile"><b>${v}</b><span>${l}</span></div>`).join("");
  $("#minG").value = qc.auto.min_genes; $("#maxG").value = qc.auto.max_genes; $("#maxMt").value = qc.auto.max_mt;
  const draw = () => {
    histSvg($("#histGenes"), qc.hist.genes, qc.hist.genes_edges, [+$("#minG").value, +$("#maxG").value]);
    histSvg($("#histMt"), qc.hist.mt, qc.hist.mt_edges, [+$("#maxMt").value]);
  };
  ["#minG", "#maxG", "#maxMt"].forEach(s => $(s).oninput = draw);
  draw();
  $("#qcWhy").textContent = qc.has_mt
    ? `Cut-offs start from your data: median ± 5 MAD on genes per cell, median + 3 MAD on mitochondrial percent (clamped to 5–25%). Dashed lines mark them. Raise the mitochondrial limit rather than lose a population — dissociation-stressed tissue genuinely runs higher.`
    : `No mitochondrial genes were found in this matrix, so no cell is removed on that criterion. Gene cut-offs are median ± 5 MAD on your data.`;
}

/* ---------------------------------------------------------------- run */
const STEPS = [["Filtering cells and genes", 10], ["Normalising and finding variable genes", 28],
  ["Principal components", 48], ["Neighbour graph and Leiden clusters", 64],
  ["UMAP layout", 80], ["Marker genes", 92]];
function stepUI(i, pct, msg) {
  $("#runBar").style.width = pct + "%"; $("#runMsg").textContent = msg;
  $("#runSteps").innerHTML = STEPS.map(([l], j) =>
    `<li class="${j < i ? "done" : j === i ? "now" : ""}">${j < i ? "✓" : "·"} ${l}</li>`).join("");
}
const breathe = () => new Promise(r => setTimeout(r, 30));   // let the browser paint between steps

$("#runBtn").onclick = async () => {
  const py = S.py, g = py.globals;
  $("#qcCard").classList.add("hidden"); $("#runCard").classList.remove("hidden");
  try {
    stepUI(0, 6, "Filtering…"); await breathe();
    g.set("_mn", +$("#minG").value); g.set("_mx", +$("#maxG").value); g.set("_mt", +$("#maxMt").value);
    const f = py.runPython("pipeline.apply_qc(int(_mn), int(_mx), float(_mt))").toJs({dict_converter: Object.fromEntries});

    stepUI(1, 28, `${f.cells_kept.toLocaleString()} cells kept · normalising`); await breathe();
    g.set("_nh", +$("#nHvg").value);
    const hv = py.runPython("pipeline.normalize_hvg_pca(int(_nh))").toJs({dict_converter: Object.fromEntries});

    stepUI(3, 56, "Neighbours and clusters…"); await breathe();
    g.set("_np", hv.knee); g.set("_k", +$("#kNN").value); g.set("_r", +$("#res").value);
    const cl = py.runPython("pipeline.neighbors_and_cluster(int(_np), int(_k), float(_r))").toJs({dict_converter: Object.fromEntries});

    stepUI(4, 72, `${cl.n_clusters} clusters · laying out the UMAP`); await breathe();
    const n = f.cells_kept, k = cl.k;
    const idx = [], dst = [];
    for (let i = 0; i < n; i++) {
      idx.push(Array.from(cl.knn_idx.slice(i * k, (i + 1) * k)));
      dst.push(Array.from(cl.knn_dist.slice(i * k, (i + 1) * k)));
    }
    const pca = py.runPython(`pipeline.STATE["pca"][:, :int(_np)].tolist()`).toJs();
    const um = new UMAP.UMAP({nNeighbors: k, minDist: 0.5, nComponents: 2, spread: 1.0,
                              random: mulberry32(42)});
    if (um.setPrecomputedKNN) um.setPrecomputedKNN(idx, dst);  // reuse the graph the clustering used
    const xy = await um.fitAsync(pca, e => {
      if (e % 40 === 0) stepUI(4, 72 + Math.min(8, e / 50), "UMAP layout…");
    });

    stepUI(5, 88, "Marker genes…"); await breathe();
    g.set("_xy", xy.flat());
    py.runPython("pipeline.set_umap(_xy.to_py())");
    const mk = py.runPython("pipeline.markers()").toJs({dict_converter: Object.fromEntries});
    S.emb = JSON.parse(py.runPython("pipeline.embedding_payload()"));
    S.mk = mk;
    stepUI(6, 100, "Done");
    showResult(f, hv, cl, mk);
  } catch (e) {
    $("#runMsg").innerHTML = `<b style="color:var(--up)">Stopped:</b> ${String(e.message || e).slice(0, 400)}`;
    console.error(e);
  }
};
// a seeded PRNG so the same data lays out the same way twice
function mulberry32(a) {return function () {a |= 0; a = a + 0x6D2B79F5 | 0;
  let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
  return ((t ^ t >>> 14) >>> 0) / 4294967296}}

/* ---------------------------------------------------------------- result */
function showResult(f, hv, cl, mk) {
  $("#runCard").classList.add("hidden"); $("#resCard").classList.remove("hidden");
  $("#resTiles").innerHTML = [
    [f.cells_kept.toLocaleString(), "cells after QC"], [f.cells_removed.toLocaleString(), "cells removed"],
    [cl.n_clusters, "clusters"], [hv.n_hvg.toLocaleString(), "variable genes"], [hv.knee, "PCs used"],
  ].map(([v, l]) => `<div class="tile"><b>${v}</b><span>${l}</span></div>`).join("");
  $("#geneList").innerHTML = S.emb.genes.slice(0, 4000).map(g => `<option value="${g}">`).join("");
  $("#mkTable").querySelector("tbody").innerHTML = mk.table.map(r =>
    `<tr><td>${r[0]}</td><td><b>${r[1]}</b></td><td class="mono">${r[2]}</td>
     <td class="mono">${r[3]}%</td><td class="mono">${r[4]}%</td><td class="mono">${r[5]}</td></tr>`).join("");
  const K = cl.n_clusters;
  $("#legend").innerHTML = Array.from({length: K}, (_, i) =>
    `<button data-i="${i}"><i style="background:${PALETTE[i % PALETTE.length]}"></i>${i} · ${cl.sizes[i]}</button>`).join("");
  $$("#legend button").forEach(b => b.onclick = () => {
    const i = +b.dataset.i; S.hidden.has(i) ? S.hidden.delete(i) : S.hidden.add(i);
    b.classList.toggle("off"); draw();
  });
  S.mode = "cluster"; draw();
}

const cv = $("#cv"), ctx = cv.getContext("2d");
function draw() {
  if (!S.emb) return;
  const r = cv.getBoundingClientRect();
  cv.width = r.width * devicePixelRatio; cv.height = r.height * devicePixelRatio;
  ctx.setTransform(devicePixelRatio, 0, 0, devicePixelRatio, 0, 0);
  ctx.clearRect(0, 0, r.width, r.height);
  const {x, y, cluster} = S.emb, n = x.length;
  const x0 = Math.min(...x), x1 = Math.max(...x), y0 = Math.min(...y), y1 = Math.max(...y);
  const pad = 16, sc = Math.min((r.width - 2 * pad) / (x1 - x0 || 1), (r.height - 2 * pad) / (y1 - y0 || 1));
  const px = v => pad + (v - x0) * sc, py2 = v => r.height - pad - (v - y0) * sc;
  const size = Math.max(1.4, Math.min(3.4, 2.4 * Math.sqrt(4000 / n)));
  let order = [...Array(n).keys()];
  if (S.mode !== "cluster" && S.vals) order.sort((a, b) => S.vals[a] - S.vals[b]);  // bright cells last
  const vmax = S.vals ? Math.max(...S.vals) || 1 : 1;
  for (const i of order) {
    if (S.mode === "cluster" && S.hidden.has(cluster[i])) continue;
    let col;
    if (S.mode === "cluster") col = PALETTE[cluster[i] % PALETTE.length];
    else if (S.vals[i] <= 0) col = "#d7dde0";                    // undetected cells stay background grey
    else {const t = Math.min(1, S.vals[i] / vmax); col = `rgb(${255 - t * 40},${230 - t * 200},${230 - t * 200})`}
    ctx.fillStyle = col; ctx.fillRect(px(x[i]), py2(y[i]), size, size);
  }
}
new ResizeObserver(draw).observe(cv);
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", draw);

$("#colorBy").onchange = async () => {
  const m = $("#colorBy").value;
  if (m === "gene") {$("#geneQ").focus(); return}
  S.mode = m; S.vals = m === "cluster" ? null : S.emb[m]; $("#geneMsg").textContent = ""; draw();
};
$("#geneQ").onchange = () => {
  const g = $("#geneQ").value.trim(); if (!g) return;
  S.py.globals.set("_g", g);
  const r = S.py.runPython("pipeline.gene_values(_g)").toJs({dict_converter: Object.fromEntries});
  if (!r.found) {$("#geneMsg").textContent = `${g} is not in this dataset.`; return}
  S.mode = "gene"; S.vals = r.values; $("#colorBy").value = "gene";
  $("#geneMsg").innerHTML = `<b>${r.gene}</b> detected in ${r.pct}% of cells.`;
  draw();
};

const dl = (name, text) => {const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], {type: "text/csv"})); a.download = name; a.click()};
$("#dlCells").onclick = () => {
  const {x, y, cluster, n_genes, pct_mt} = S.emb;
  dl("cells.csv", "umap_1,umap_2,cluster,n_genes,pct_mt\n" +
    x.map((_, i) => [x[i], y[i], cluster[i], n_genes[i], pct_mt[i]].join(",")).join("\n"));
};
$("#dlMarkers").onclick = () => dl("markers.csv",
  "cluster,gene,log2FC,pct_in,pct_out,z\n" + S.mk.table.map(r => r.join(",")).join("\n"));
$("#againBtn").onclick = () => location.reload();
