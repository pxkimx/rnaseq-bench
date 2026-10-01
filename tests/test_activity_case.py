"""Pathway/TF activity must score the whole network whichever case the gene symbols use.

The bug this guards against: mouse symbols were upper-cased before being scored against the mouse
PROGENy and CollecTRI networks, which are themselves mouse-case. In bulk that raised decoupler's
"No sources with more than tmin=5 targets" and the step reported itself as skipped; in single cell it
did not fail at all — it scored one or two of the fourteen pathways from the few genes that happened
to match, and presented them as the result.

Needs the cached networks (they are downloaded on first use):  python tests/test_activity_case.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from server import resources  # noqa: E402

fails = []


def check(name, got, want):
    ok = got == want
    if not ok:
        fails.append(f"{name}: got {got!r}, wanted {want!r}")
    print(f"{'ok  ' if ok else 'FAIL'} {name}")


def score(net, columns):
    """What the pipeline does: translate to the network's case, then run ULM."""
    import decoupler as dc
    rng = np.random.default_rng(0)
    df = pd.DataFrame(rng.normal(size=(4, len(columns))), index=[f"s{i}" for i in range(4)], columns=list(columns))
    df.columns = resources.match_case(df.columns, net)
    df = df.loc[:, ~df.columns.duplicated()]
    es, _ = dc.mt.ulm(data=df, net=net, tmin=5)
    return es.shape[1]


for organism in ("mouse", "human"):
    prog = resources.progeny(organism)
    if prog is None:
        print(f"-- {organism} PROGENy unavailable (no cached network and no internet) — skipping")
        continue
    native = sorted(set(prog["target"].astype(str)))
    total = prog["source"].nunique()
    print(f"\n-- {organism} PROGENy: {total} pathways, {len(native)} targets, "
          f"{sum(t == t.upper() for t in native)} of them already upper-case")

    # the data in the network's own spelling
    check(f"{organism}: all pathways scored, data in the network's case", score(prog, native), total)
    # the same data in the other case — a table that happens to use upper-case mouse symbols, or
    # title-case human ones, must work just as well
    check(f"{organism}: all pathways scored, data upper-cased", score(prog, [t.upper() for t in native]), total)
    check(f"{organism}: all pathways scored, data title-cased",
          score(prog, [t[:1].upper() + t[1:].lower() for t in native]), total)

    # and the old behaviour, so the regression is recorded rather than merely fixed: upper-casing the
    # data without translating the network loses most of mouse, and nothing of human
    import decoupler as dc
    rng = np.random.default_rng(0)
    raw = pd.DataFrame(rng.normal(size=(4, len(native))), index=[f"s{i}" for i in range(4)],
                       columns=[t.upper() for t in native])
    raw = raw.loc[:, ~raw.columns.duplicated()]
    try:
        old_n = dc.mt.ulm(data=raw, net=prog, tmin=5)[0].shape[1]
    except AssertionError:
        old_n = 0
    print(f"     (before the fix, upper-casing {organism} data scored {old_n} of {total})")
    if organism == "mouse":
        check("mouse: the old upper-casing really was broken", old_n < total / 2, True)

tri = resources.collectri("mouse")
if tri is not None:
    native = sorted(set(tri["target"].astype(str)))
    total = tri["source"].nunique()
    got = score(tri, native)
    print(f"\n-- mouse CollecTRI: {total} TFs in the network, {got} scored from {len(native)} targets")
    check("mouse CollecTRI: most TFs scored", got > total / 2, True)

# ---------------------------------------------------------------- the single-cell path, end to end
# This is the one that failed silently: it scored one or two pathways from the handful of genes that
# happened to match and presented them as the result, with nothing to say it was incomplete.
prog = resources.progeny("mouse")
if prog is not None:
    import shutil
    import tempfile

    import anndata as ad

    from server import sc_extra
    from server.common import Job

    targets = sorted(set(prog["target"].astype(str)))
    rng = np.random.default_rng(0)
    tmp = Path(tempfile.mkdtemp(prefix="rnaseq-bench-activity-"))

    def sc_run(genes, label):
        out = tmp / "job"
        shutil.rmtree(out, ignore_errors=True)
        out.mkdir(parents=True)
        job = Job(out, "sc", label)
        job.section("activity", "Step", "Pathway activity", "test")
        X = rng.lognormal(size=(400, len(genes))).astype("float32")
        a = ad.AnnData(X)
        a.var_names = list(genes)
        a.obs["cell_type_simple"] = pd.Categorical(rng.choice(["A", "B", "C"], 400))
        a.layers["counts"] = X.copy()
        ok = sc_extra.pathway_activity(job, a, "mouse")
        n = pd.read_csv(out / "pathway_activity_by_cluster.csv", index_col=0).shape[1] if ok else 0
        warns = [f["text"] for f in job.result["flags"] if f["level"] == "warn"]
        return n, warns

    print("\n-- single cell, mouse PROGENy")
    n, w = sc_run(targets[:3000], "mouse-case")
    check("single cell: mouse-case matrix scores all 14", (n, bool(w)), (14, False))
    n, w = sc_run([g.upper() for g in targets[:3000]], "upper-case")
    check("single cell: upper-case mouse matrix scores all 14", (n, bool(w)), (14, False))

    # and a matrix that really does have too few of the network's genes must say so rather than
    # quietly return a handful of pathways
    sparse = list(rng.choice(targets, 40, replace=False)) + [f"NOTAGENE{i}" for i in range(1500)]
    n, w = sc_run(sparse, "sparse")
    check("single cell: a partial result is flagged", n < 7 and bool(w), True)
    shutil.rmtree(tmp, ignore_errors=True)

print()
if fails:
    print(f"{len(fails)} failed:")
    for f in fails:
        print("  " + f)
    sys.exit(1)
print("all activity-case checks passed")
