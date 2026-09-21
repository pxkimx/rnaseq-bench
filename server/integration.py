"""Judging how hard a dataset is to integrate, and doing it with Harmony or scVI.

The choice between the two is not a matter of taste. Harmony corrects coordinates after PCA and is
fast; scVI models the raw counts with sequencing depth and batch as explicit terms and learns a
representation in which the batch was never a factor, at the cost of training a network. Which one a
dataset needs is something you can measure before running either, so the app measures it and says so
rather than leaving the user to guess.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .common import log_exc


def assess(adata, batch_key: str, n_pcs: int, k: int = 30) -> dict:
    """How separated are the batches, before any correction?

    Two numbers decide it:

    `excess` — of a cell's k nearest neighbours in PCA space, how many more come from its own batch
    than batch sizes alone would predict. 0 means the batches are already mixed and correction has
    nothing to do; 1 means every cell sits exclusively among its own batch.

    `depth_ratio` — the spread of median genes per cell across batches. This is the case Harmony is
    weakest on: it moves cells in PCA space, but a library sequenced a third as deeply differs in how
    much was detected, not in where it sits, and no rotation of the components fixes that. scVI carries
    a per-cell library-size term, so it treats depth as something to explain rather than to undo.
    """
    out = {"n_batches": 0, "excess": 0.0, "depth_ratio": 1.0, "smallest": 0}
    try:
        b = adata.obs[batch_key].astype(str).values
        levels, counts = np.unique(b, return_counts=True)
        out["n_batches"] = int(len(levels))
        out["smallest"] = int(counts.min())
        if len(levels) < 2:
            return out

        med = adata.obs.groupby(batch_key, observed=True)["n_genes_by_counts"].median()
        out["depth_ratio"] = float(med.max() / max(med.min(), 1))

        # sub-sample: the neighbour search is only a diagnostic and must not cost more than the analysis
        rng = np.random.default_rng(0)
        idx = rng.choice(adata.n_obs, size=min(adata.n_obs, 4000), replace=False)
        X = np.asarray(adata.obsm["X_pca"][idx, :n_pcs])
        from sklearn.neighbors import NearestNeighbors
        nn = NearestNeighbors(n_neighbors=min(k + 1, len(idx))).fit(X)
        _, ind = nn.kneighbors(X)
        bi = pd.Categorical(b[idx]).codes
        same = (bi[ind[:, 1:]] == bi[:, None]).mean(1)          # own-batch share among neighbours
        share = pd.Series(bi).value_counts(normalize=True).sort_index().values
        expected = share[bi]                                    # what chance alone would give
        out["excess"] = float(np.clip(((same - expected) / np.maximum(1 - expected, 1e-9)).mean(), 0, 1))
    except Exception:  # noqa: BLE001
        log_exc("integration assess")
    return out


def recommend(m: dict, n_cells: int, scvi_available: bool) -> tuple[str, str]:
    """(method, the reason in plain language). Method is 'none', 'harmony' or 'scvi'."""
    n, ex, dr = m["n_batches"], m["excess"], m["depth_ratio"]
    if n < 2:
        return "none", "There is only one batch, so there is nothing to integrate."
    mixed = f"{100*ex:.0f}% of a cell's neighbours are from its own batch beyond what batch sizes alone explain"
    depth = f"the deepest batch detects {dr:.1f}× as many genes per cell as the shallowest"

    if ex < 0.15 and dr < 1.6:
        return "none", (f"The batches are already well mixed ({mixed}) and sequenced to a similar depth "
                        f"({depth}), so correction would have little to do and risks removing real signal.")
    # Near-total separation is ambiguous: it is what a real batch effect looks like, and also what two
    # genuinely different tissues or conditions look like. Integrating the second case erases the thing
    # being studied, so say so rather than quietly merging them.
    caution = (" Note that the batches barely overlap at all. That is what a strong batch effect looks "
               "like, but it is also what two genuinely different samples look like — different tissues, "
               "timepoints or conditions. If these are not meant to contain the same cell types, "
               "integrating them will remove the difference you are studying." if ex >= 0.9 else "")
    hard = ex >= 0.45 or dr >= 2.0
    if hard and scvi_available:
        return "scvi", (f"These batches are strongly separated — {mixed}, and {depth}. Harmony moves cells "
                        "within the principal components, which cannot repair a difference in how much was "
                        "detected; scVI models counts with per-cell library size and batch as explicit "
                        "terms, so it handles this case better. It trains a network, so expect minutes "
                        "rather than seconds." + caution)
    if hard:
        return "harmony", (f"These batches are strongly separated — {mixed}, and {depth}. Harmony is being "
                           "used because scvi-tools is not installed; for a difference in depth this large "
                           "scVI would usually do better, since it models library size directly." + caution)
    return "harmony", (f"The batches overlap moderately ({mixed}) at a similar depth ({depth}). Harmony "
                       "handles this well and takes seconds; scVI would cost minutes for little gain.")


def run_scvi(adata, batch_key: str, n_latent: int = 30, max_epochs: int | None = None):
    """Train scVI on the raw counts of the variable genes and return the latent representation."""
    import scvi

    scvi.settings.seed = 0
    sub = adata[:, adata.var["highly_variable"]].copy() if "highly_variable" in adata.var else adata.copy()
    sub.X = sub.layers["counts"].copy()          # scVI wants raw counts, not the log-normalised matrix
    scvi.model.SCVI.setup_anndata(sub, batch_key=batch_key, layer=None)
    model = scvi.model.SCVI(sub, n_latent=n_latent)
    model.train(max_epochs=max_epochs, early_stopping=True, enable_progress_bar=False)
    return np.asarray(model.get_latent_representation())
