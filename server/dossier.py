"""Gene dossier: what the public databases say about one gene, fetched with Biopython.

Clicking a gene anywhere in a result opens this. It answers "what is this gene?" without leaving the
analysis: NCBI Gene (name, location, RefSeq summary), UniProt (function, localisation, domains, disease),
KEGG (pathways, linked to the map with the gene highlighted), AlphaFold DB (how much of the predicted
structure is confident) and PubMed (most relevant and newest papers).

Every source is independent and optional. Whatever cannot be reached becomes a note naming that source,
and the rest of the dossier still shows — the same rule as every other internet step in this app. A dossier
that fetched completely is cached for two weeks, so reopening a gene is instant and works offline; a partial
one is not cached, so the missing parts are tried again next time.
"""
from __future__ import annotations

import hashlib
import json
import re
import threading
import time
import urllib.error
import urllib.request
import warnings
from concurrent.futures import Future, ThreadPoolExecutor, wait

from .common import log_exc
from .geo import CACHE, load_settings

DOSSIER_DIR = CACHE / "dossier"
TTL = 14 * 86400
TIMEOUT = 20          # per request, where the client lets us set one
DEADLINE = 25         # per source: Entrez and alphafold_db open URLs with no timeout, so we stop waiting instead
_POOL = ThreadPoolExecutor(max_workers=8, thread_name_prefix="dossier")
STRUCTURE_BENCH = "http://localhost:8767"

# taxon id, KEGG organism code, name. KEGG gene ids for these organisms are the NCBI Gene ids.
SPECIES = {
    "human": (9606, "hsa", "Homo sapiens"),
    "mouse": (10090, "mmu", "Mus musculus"),
    "rat": (10116, "rno", "Rattus norvegicus"),
}


def guess_species(gene: str) -> str:
    """Mouse and rat symbols are capitalised (Lmna), human ones are upper case (LMNA) — except the human
    open-reading-frame names (C1orf112), whose lower-case 'orf' is not a mouse sign."""
    return "human" if gene == gene.upper() or re.search(r"orf\d", gene) else "mouse"


def _entrez():
    from Bio import Entrez
    cfg = load_settings()
    Entrez.tool = "RNAseqBench"
    # Biopython retries three times, 15 s apart — right for a batch script, far too long for someone waiting
    # on a drawer to open. _ecall does the pacing and retrying instead.
    Entrez.max_tries = 1
    Entrez.email = cfg.get("ncbi_email") or None
    Entrez.api_key = cfg.get("ncbi_api_key") or None
    return Entrez


def _read(handle):
    Entrez = _entrez()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")        # "Email address is not specified" — optional, set in Settings
        try:
            return Entrez.read(handle)
        finally:
            handle.close()


_ELOCK = threading.Lock()
_ELAST = [0.0]


def _ecall(fn, **kw):
    """One NCBI request, paced for every thread in the server together.

    NCBI allows 3 requests a second without an API key (10 with one) and answers HTTP 429 beyond that.
    Biopython paces calls with a module-level timestamp that is not locked, so the dossier's parallel NCBI and
    PubMed lookups — or two dossiers opened in quick succession — went over the limit and came back empty.
    The request is opened inside the lock; reading the response happens outside it.
    """
    Entrez = _entrez()
    gap = 0.12 if Entrez.api_key else 0.36
    for attempt in range(4):
        with _ELOCK:
            wait_s = _ELAST[0] + gap - time.time()
            if wait_s > 0:
                time.sleep(wait_s)
            try:
                return _call(fn, **kw)
            except urllib.error.HTTPError as e:
                if e.code != 429 or attempt == 3:
                    raise
            finally:
                _ELAST[0] = time.time()
        time.sleep(0.8 * (attempt + 1))       # over the limit anyway: back off, then queue again
    raise RuntimeError("unreachable")


def _call(fn, *a, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return fn(*a, **kw)


# ---------------------------------------------------------------- sources
def ncbi_gene(gene: str, taxid: int) -> dict:
    Entrez = _entrez()
    ids = _read(_ecall(Entrez.esearch, db="gene", term=f"{gene}[sym] AND {taxid}[taxid]"))["IdList"]
    via_alias = False
    if not ids:
        # an outdated symbol, an alias or an Ensembl id: search every field and say that we did
        ids = _read(_ecall(Entrez.esearch, db="gene",
                          term=f'"{gene}"[All Fields] AND {taxid}[taxid] AND alive[prop]'))["IdList"]
        via_alias = bool(ids)
    if not ids:
        return {"found": False}
    gid = ids[0]
    d = _read(_ecall(Entrez.esummary, db="gene", id=gid))["DocumentSummarySet"]["DocumentSummary"][0]
    aliases = [a.strip() for a in str(d.get("OtherAliases", "")).split(",") if a.strip()]
    # the official symbol; "Name" can differ (MT-ND4's is "ND4")
    symbol = str(d.get("NomenclatureSymbol") or d.get("Name") or gene)
    via_alias = via_alias or symbol.lower() != gene.lower()     # [sym] also matches aliases (HGPS → LMNA)
    return {
        "found": True, "gene_id": gid, "symbol": symbol, "name": str(d.get("Description", "")),
        "summary": str(d.get("Summary", "")).strip(), "chromosome": str(d.get("Chromosome", "")),
        "location": str(d.get("MapLocation", "")), "aliases": aliases[:12],
        "mim": [str(m) for m in d.get("Mim", [])][:3], "via_alias": via_alias,
        "url": f"https://www.ncbi.nlm.nih.gov/gene/{gid}",
    }


def pubmed(gene: str, context: str = "") -> dict:
    Entrez = _entrez()
    term = f"{gene}[tiab]" + (f" AND ({context})" if context else "")
    top = _read(_ecall(Entrez.esearch, db="pubmed", term=term, retmax=6, sort="relevance"))
    new = _read(_ecall(Entrez.esearch, db="pubmed", term=term, retmax=4, sort="pub_date"))
    ids = list(dict.fromkeys(list(top["IdList"]) + list(new["IdList"])))
    papers = {}
    if ids:
        for s in _read(_ecall(Entrez.esummary, db="pubmed", id=",".join(ids))):
            authors = list(s.get("AuthorList", []))
            papers[str(s["Id"])] = {
                "pmid": str(s["Id"]), "title": str(s.get("Title", "")).rstrip("."),
                "journal": str(s.get("Source", "")), "year": str(s.get("PubDate", ""))[:4],
                "authors": (authors[0] + (" et al." if len(authors) > 1 else "")) if authors else "",
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{s['Id']}/",
            }
    return {
        "query": term, "count": int(top["Count"]),
        "relevant": [papers[i] for i in top["IdList"] if i in papers],
        "newest": [papers[i] for i in new["IdList"] if i in papers and i not in top["IdList"]],
        "url": "https://pubmed.ncbi.nlm.nih.gov/?term=" + urllib.request.quote(term),
    }


_REF = re.compile(r"\s*\((?:PubMed|ECO|Ref\.)[^)]*\)")


def _texts(comment: dict) -> list[str]:
    return [_REF.sub("", t["value"]).strip() for t in comment.get("texts", []) if t.get("value")]


def uniprot(gene: str, taxid: int) -> dict:
    from Bio import UniProt
    fields = ["accession", "protein_name", "length", "gene_primary", "reviewed"]
    hits = _call(UniProt.search, f"gene_exact:{gene} AND organism_id:{taxid} AND reviewed:true", fields=fields)
    reviewed = True
    if not len(hits):
        hits = _call(UniProt.search, f"gene_exact:{gene} AND organism_id:{taxid}", fields=fields)
        reviewed = False
    if not len(hits):
        return {"found": False}
    acc = hits[0]["primaryAccession"]
    with urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=TIMEOUT) as r:
        e = json.load(r)
    function, locations, diseases = [], [], []
    for c in e.get("comments", []):
        t = c.get("commentType")
        if t == "FUNCTION":
            function += _texts(c)
        elif t == "SUBCELLULAR LOCATION":
            locations += [s["location"]["value"] for s in c.get("subcellularLocations", []) if "location" in s]
        elif t == "DISEASE" and "disease" in c:
            dz = c["disease"]
            diseases.append({"name": dz.get("diseaseId", ""), "acronym": dz.get("acronym", "")})
    keep = {"Domain", "Zinc finger", "Transmembrane", "Coiled coil", "Motif", "DNA binding", "Signal", "Repeat"}
    domains = [{"type": f["type"], "name": f.get("description") or f["type"],
                "start": f["location"]["start"].get("value"), "end": f["location"]["end"].get("value")}
               for f in e.get("features", []) if f["type"] in keep][:24]
    pdb = [x["id"] for x in e.get("uniProtKBCrossReferences", []) if x.get("database") == "PDB"]
    desc = e.get("proteinDescription", {})
    name = (desc.get("recommendedName") or (desc.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value", "")
    return {
        "found": True, "accession": acc, "reviewed": reviewed, "name": name,
        "length": e.get("sequence", {}).get("length"), "function": " ".join(function)[:1600],
        "locations": list(dict.fromkeys(locations))[:8], "diseases": diseases[:12], "domains": domains,
        "n_pdb": len(pdb), "pdb": pdb[:12], "url": f"https://www.uniprot.org/uniprotkb/{acc}/entry",
    }


def _kegg_names(org: str) -> dict:
    p = DOSSIER_DIR / f"kegg_pathways_{org}.tsv"
    if not p.exists() or time.time() - p.stat().st_mtime > 90 * 86400:
        from Bio.KEGG import REST
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(_call(REST.kegg_list, "pathway", org).read())
    out = {}
    for line in p.read_text().splitlines():
        if "\t" in line:
            pid, name = line.split("\t", 1)
            out[pid.replace("path:", "")] = re.sub(r" - [A-Z][a-z]+ [a-z]+ \(.*\)$", "", name)
    return out


def kegg(gene_id: str, org: str) -> dict:
    from Bio.KEGG import REST
    kid = f"{org}:{gene_id}"
    text = _call(REST.kegg_link, "pathway", kid).read()
    pids = [ln.split("\t")[1].replace("path:", "") for ln in text.splitlines() if "\t" in ln]
    names = _kegg_names(org) if pids else {}
    # the global overview maps (01100 "Metabolic pathways" …) list thousands of genes and say nothing specific
    rows = [{"id": p, "name": names.get(p, p), "url": f"https://www.kegg.jp/pathway/{p}+{kid}"}
            for p in pids if not p[len(org):].startswith("011") and not p[len(org):].startswith("012")]
    return {"kegg_id": kid, "pathways": rows, "url": f"https://www.kegg.jp/entry/{kid}"}


def alphafold(acc: str) -> dict:
    from Bio.PDB import alphafold_db
    try:
        preds = list(_call(alphafold_db.get_predictions, acc))
    except urllib.error.HTTPError as e:
        if e.code == 404:          # the database answered: it has no model for this entry
            return {"found": False}
        raise
    if not preds:
        return {"found": False}
    p = next((x for x in preds if x.get("entryId") == f"AF-{acc}-F1"), preds[0])
    fr = {k: p.get(f"fractionPlddt{k}") for k in ("VeryHigh", "Confident", "Low", "VeryLow")}
    return {
        "found": True, "entry": p.get("entryId"), "mean_plddt": p.get("globalMetricValue"),
        "fractions": fr, "version": p.get("latestVersion"), "isoforms": len(preds) - 1,
        "url": f"https://alphafold.ebi.ac.uk/entry/{acc}",
    }


# ---------------------------------------------------------------- assembly
def _explain_af(af: dict) -> str:
    fr = af["fractions"]
    if None in fr.values():
        return f"Mean pLDDT {af['mean_plddt']:.0f}."
    good = 100 * (fr["VeryHigh"] + fr["Confident"])
    poor = 100 * (fr["Low"] + fr["VeryLow"])
    s = f"Mean pLDDT {af['mean_plddt']:.0f}: {good:.0f}% of residues are modelled confidently (pLDDT ≥ 70)"
    if poor >= 30:
        s += (f", while {poor:.0f}% are low confidence — usually regions that are disordered or flexible on their "
              "own, and not a structure to read detail from")
    return s + "."


def build(gene: str, species: str = "auto", context: str = "", refresh: bool = False) -> dict:
    gene = gene.strip()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._\-/@]{0,40}", gene):
        raise ValueError(f"'{gene[:40]}' does not look like a gene symbol.")
    if species in (None, "", "auto"):
        species = guess_species(gene)
    if species not in SPECIES:
        raise ValueError(f"Species must be one of {', '.join(SPECIES)}.")
    context = (context or "").strip()[:120]
    taxid, org, sci = SPECIES[species]
    key = f"{gene}__{hashlib.sha1(context.encode()).hexdigest()[:8]}" if context else gene
    cp = DOSSIER_DIR / species / f"{key}.json"
    if not refresh and cp.exists() and time.time() - cp.stat().st_mtime < TTL:
        d = json.loads(cp.read_text())
        d["cached"] = True
        return d
    # a double-click, or the assistant asking while the drawer loads: share the fetch rather than repeat it
    with _FLIGHT_LOCK:
        fut = _FLIGHT.get(str(cp))
        owner = fut is None
        if owner:
            fut = _FLIGHT[str(cp)] = Future()
    if not owner:
        return fut.result(timeout=3 * DEADLINE + 10)
    try:
        d = _fetch(gene, species, context, taxid, org, sci, cp)
        fut.set_result(d)
        return d
    except BaseException as e:
        fut.set_exception(e)
        raise
    finally:
        with _FLIGHT_LOCK:
            _FLIGHT.pop(str(cp), None)


_FLIGHT: dict = {}
_FLIGHT_LOCK = threading.Lock()


def _fetch(gene, species, context, taxid, org, sci, cp) -> dict:
    t0 = time.time()
    notes, parts = [], {}

    def guarded(name, label, fn, *a):
        try:
            parts[name] = fn(*a)
        except Exception:  # noqa: BLE001 - one unreachable database must not blank the dossier
            msg = log_exc(f"dossier {name} {gene}")
            parts[name] = None
            notes.append({"source": name, "text": f"{label} could not be reached ({msg}). The rest is unaffected; "
                                                  "open the dossier again later to retry."})

    def run(*jobs):
        """Run (name, label, fn, *args) jobs side by side; a source that has not answered by DEADLINE is
        reported as not answering. Its thread is abandoned rather than awaited (the pool is shared, not
        closed), so a hung connection costs a worker, not the user's wait."""
        fs = {_POOL.submit(guarded, *j): j for j in jobs}
        done, pending = wait(fs, timeout=DEADLINE)
        for f in pending:
            name, label = fs[f][:2]
            parts[name] = None
            notes.append({"source": name, "text": f"{label} did not answer within {DEADLINE} s. The rest is "
                                                  "unaffected; open the dossier again later to retry."})

    # NCBI first: it resolves aliases, and its gene id is what KEGG uses. The others run alongside KEGG.
    run(("ncbi", "NCBI Gene", ncbi_gene, gene, taxid))
    ncbi = parts.get("ncbi") or {}
    symbol = ncbi.get("symbol") if ncbi.get("found") else gene
    jobs = [("uniprot", "UniProt", uniprot, symbol, taxid), ("pubmed", "PubMed", pubmed, symbol, context)]
    if ncbi.get("found"):
        jobs.append(("kegg", "KEGG", kegg, ncbi["gene_id"], org))
    run(*jobs)
    up = parts.get("uniprot") or {}
    if up.get("found"):
        run(("alphafold", "AlphaFold DB", alphafold, up["accession"]))

    # plain-language notes about what came back, in the voice of the rest of the app
    if parts.get("ncbi") is not None and not ncbi.get("found"):
        notes.append({"source": "ncbi", "text": f"NCBI Gene has no {species} gene called {gene}. If this is an "
                                                "Ensembl or transcript id, it was not mapped to a symbol."})
    elif ncbi.get("via_alias"):
        kind = ("a database id" if re.match(r"(ENS[A-Z]*G\d|[NX][MR]_\d)", gene) else
                "an alias or former symbol")
        notes.append({"source": "ncbi", "text": f"{gene} is {kind}, not the official {species} symbol; NCBI Gene "
                                                f"matches it to <b>{ncbi['symbol']}</b>, which is what is shown."})
    if ncbi.get("found") and not ncbi.get("summary"):
        notes.append({"source": "ncbi", "text": "NCBI has no RefSeq summary for this gene (common for mouse genes); "
                                                "the UniProt function below is usually the better description."})
    if parts.get("uniprot") is not None and not up.get("found"):
        notes.append({"source": "uniprot", "text": "No UniProt entry for this gene — it may be non-coding."})
    elif up.get("found") and not up.get("reviewed"):
        noncoding = re.search(r"transcript|non-protein coding|long intergenic|antisense|microRNA|small nucleolar",
                              ncbi.get("name", ""), re.I)
        notes.append({"source": "uniprot", "text": "Only an unreviewed (TrEMBL) UniProt entry exists, so the function "
                                                   "text is computer-annotated." + (
            f" NCBI describes {symbol} as <b>{ncbi['name']}</b> — a non-coding RNA — so this predicted protein is "
            "probably not real." if noncoding else "")})
    af = parts.get("alphafold") or {}
    if parts.get("alphafold") is not None and not af.get("found"):
        long_ = (up.get("length") or 0) > 2700
        notes.append({"source": "alphafold", "text": f"AlphaFold DB has no model for {up['accession']}. " + (
            f"At {up['length']:,} residues it is longer than the 2,700 AlphaFold DB models as one piece; it holds such "
            "proteins only as overlapping fragments in bulk downloads." if long_ else
            "It covers most reviewed proteins, so this is usual only for computer-predicted or very recent entries.")})
    if af.get("found"):
        af["explain"] = _explain_af(af)
    pm = parts.get("pubmed") or {}
    if pm and pm.get("count", 0) > 20000 and not context:
        notes.append({"source": "pubmed", "text": f"{pm['count']:,} papers mention {symbol}. Add a context term below "
                                                  "(your tissue, cell type or condition) to find the ones that bear on "
                                                  "this analysis — and if the symbol is also an ordinary word or "
                                                  "abbreviation, many of these are about something else."})

    d = {"gene": gene, "symbol": symbol, "species": species, "organism": sci, "taxid": taxid, "context": context,
         "fetched": time.time(), "elapsed": round(time.time() - t0, 1), "cached": False, "notes": notes, **parts,
         "structure_bench": f"{STRUCTURE_BENCH}/#" + (f"uniprot={up['accession']}" if up.get("found")
                                                      else f"gene={symbol}&species={species}")}
    if not any(n["text"].endswith("to retry.") for n in notes):
        cp.parent.mkdir(parents=True, exist_ok=True)
        cp.write_text(json.dumps(d))
    return d


def compact(d: dict) -> dict:
    """The dossier without links and long lists — what the assistant needs to answer from it."""
    out = {k: d.get(k) for k in ("symbol", "species", "notes")}
    if d.get("ncbi"):
        out["ncbi"] = {k: d["ncbi"].get(k) for k in ("name", "summary", "location", "aliases", "gene_id")}
    if d.get("uniprot"):
        out["uniprot"] = {k: d["uniprot"].get(k) for k in ("accession", "name", "function", "locations", "diseases",
                                                          "domains", "n_pdb")}
    if d.get("kegg"):
        out["kegg_pathways"] = [p["name"] for p in d["kegg"]["pathways"]]
    if d.get("alphafold"):
        out["alphafold"] = {k: d["alphafold"].get(k) for k in ("mean_plddt", "fractions", "explain")}
    if d.get("pubmed"):
        out["pubmed"] = {"query": d["pubmed"]["query"], "count": d["pubmed"]["count"],
                         "papers": [{k: p[k] for k in ("pmid", "title", "journal", "year")}
                                    for p in d["pubmed"]["relevant"] + d["pubmed"]["newest"]]}
    return out
