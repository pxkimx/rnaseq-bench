"""Offline checks for the gene dossier (server/dossier.py). Run: python tests/test_dossier.py

The dossier's job is to never fail as a whole: every source that cannot be reached, or does not answer,
must turn into a note while the rest still shows. These checks simulate both — no network is used.
"""
import sys
import tempfile
import time
import urllib.error
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from server import dossier  # noqa: E402

dossier.DOSSIER_DIR = Path(tempfile.mkdtemp())      # never touch the real cache


def check(cond, what):
    print(("ok   " if cond else "FAIL ") + what)
    if not cond:
        sys.exit(1)


check(dossier.guess_species("LMNA") == "human", "upper-case symbol is read as human")
check(dossier.guess_species("Lmna") == "mouse", "capitalised symbol is read as mouse")
check(dossier.guess_species("C1orf112") == "human", "C1orf112 stays human")
try:
    dossier.build("<script>")
    check(False, "a non-symbol is refused")
except ValueError:
    check(True, "a non-symbol is refused")


def offline(*a, **k):
    raise urllib.error.URLError("offline (simulated)")


real = {n: getattr(dossier, n) for n in ("ncbi_gene", "uniprot", "pubmed", "kegg", "alphafold")}
for n in real:
    setattr(dossier, n, offline)
d = dossier.build("TP53", "human", refresh=True)
check(d["ncbi"] is None and d["uniprot"] is None and d["pubmed"] is None, "offline: every source is empty, none raises")
check({n["source"] for n in d["notes"]} == {"ncbi", "uniprot", "pubmed"}, "offline: one note per source tried")
check(not list(dossier.DOSSIER_DIR.rglob("*.json")), "offline: a partial dossier is not cached")

# a source that hangs is abandoned at the deadline instead of holding the drawer open
dossier.DEADLINE = 1
dossier.ncbi_gene = lambda *a: {"found": True, "gene_id": "7157", "symbol": "TP53", "summary": "x"}
dossier.uniprot = lambda *a: time.sleep(5)
dossier.pubmed = lambda *a: {"query": "TP53[tiab]", "count": 3, "relevant": [], "newest": [], "url": ""}
dossier.kegg = lambda *a: {"kegg_id": "hsa:7157", "pathways": [], "url": ""}
t = time.time()
d = dossier.build("TP53", "human", refresh=True)
check(time.time() - t < 3, f"hung source abandoned after the deadline ({time.time() - t:.1f} s)")
check(any(n["source"] == "uniprot" and "did not answer" in n["text"] for n in d["notes"]), "hung source gets a note")
check(d["pubmed"]["count"] == 3 and d["kegg"]["kegg_id"] == "hsa:7157", "the sources that answered still show")
# identical requests at the same time share one fetch (double-click, or the assistant asking meanwhile)
from concurrent.futures import ThreadPoolExecutor  # noqa: E402
calls = []
dossier.DEADLINE = 5
dossier.uniprot = lambda *a: (calls.append(1), time.sleep(0.5), {"found": False})[-1]
with ThreadPoolExecutor(6) as ex:
    outs = list(ex.map(lambda _: dossier.build("TP53", "human", refresh=True), range(6)))
check(len(calls) == 1 and all(o["symbol"] == "TP53" for o in outs), f"6 simultaneous requests → {len(calls)} fetch")

c = dossier.compact(d)
check("pubmed" in c and "url" not in c["pubmed"], "compact() drops links for the assistant")
print("all dossier checks passed")
