"""Declarative description of the parameters a pipeline actually used.

`result["params"]` carries the *values* (the agent, codegen and the re-run endpoint read it).
`result["params_panel"]` carries the same values wrapped in everything the UI needs to show a
Parameters panel: a label, the control to edit it with, where the value came from and — the point
of the whole thing — a plain-language reason it has that value.

Keep the two in sync: every editable field here must be a key the pipeline's `run(params=...)`
understands, otherwise the Re-analyze button will silently drop it.
"""
from __future__ import annotations


def F(key, label, value, why, type="text", source="default", **extra) -> dict:
    """One parameter row.

    source: "user"    — you set this by hand, so it is used as given
            "auto"    — derived from this dataset (the why says how)
            "data"    — read off the file/metadata (column names, levels)
            "default" — the built-in default, which suited this dataset
    """
    d = {"key": key, "label": label, "value": value, "why": why, "type": type, "source": source}
    d.update({k: v for k, v in extra.items() if v is not None})
    return d


def G(id, title, fields, note="") -> dict:
    """A group of parameter rows; renders as one titled block in the panel."""
    return {"id": id, "title": title, "note": note, "fields": [f for f in fields if f]}


def src(user_params: dict, key: str, auto=False) -> str:
    """Where did this value come from? `auto=True` when the pipeline derived it from the data."""
    v = (user_params or {}).get(key)
    if v not in (None, "", [], "auto"):
        return "user"
    return "auto" if auto else "default"


def onoff(user_params: dict, key: str, default=True) -> str:
    return "user" if (user_params or {}).get(key) not in (None, "") and (user_params or {}).get(key) != default else "default"


def skipped(flags, *needles) -> str:
    """The '… skipped (…)' note a pipeline logged for an optional analysis, if it logged one.

    An add-on toggle that reads "on" while the analysis never ran is a lie the panel must not tell.
    """
    for f in flags or []:
        t = str(f.get("text", ""))
        if "skipped" in t.lower() and any(n.lower() in t.lower() for n in needles):
            return " " + t.rstrip(".").rstrip() + "."
    return ""


# obs columns the pipeline itself adds — never offer them as a batch / sample / condition choice
DERIVED = {"qc_pass", "predicted_doublet", "doublet_score", "leiden", "cell_type", "cell_type_simple",
           "phase", "S_score", "G2M_score", "dpt_pseudotime", "celltypist", "celltypist_conf"}


def real_obs(cols) -> list:
    return [c for c in (cols or []) if c not in DERIVED and not str(c).startswith("leiden")]
