"""Build a publication-style PDF report from a finished job's result.json."""
from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

from PIL import Image
from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (Image as RLImage, KeepTogether, PageBreak, Paragraph, SimpleDocTemplate,
                                Spacer, Table, TableStyle)

INK = colors.HexColor("#13212b")
MUTED = colors.HexColor("#5b6b75")
ACCENT = colors.HexColor("#0f766e")
WARN = colors.HexColor("#b45309")
LINE = colors.HexColor("#d8e0e4")

S = {
    "title": ParagraphStyle("title", fontName="Helvetica-Bold", fontSize=24, leading=28, textColor=INK),
    "kicker": ParagraphStyle("kicker", fontName="Helvetica-Bold", fontSize=8, leading=10, textColor=ACCENT, spaceAfter=2),
    "h1": ParagraphStyle("h1", fontName="Helvetica-Bold", fontSize=16, leading=20, textColor=INK, spaceAfter=4),
    "h2": ParagraphStyle("h2", fontName="Helvetica-Bold", fontSize=11, leading=14, textColor=INK, spaceBefore=4, spaceAfter=2),
    "body": ParagraphStyle("body", fontName="Helvetica", fontSize=9, leading=12.5, textColor=INK, alignment=TA_LEFT),
    "muted": ParagraphStyle("muted", fontName="Helvetica", fontSize=8.5, leading=11.5, textColor=MUTED),
    "sub": ParagraphStyle("sub", fontName="Helvetica-Oblique", fontSize=7.5, leading=10, textColor=MUTED),
    "cell": ParagraphStyle("cell", fontName="Helvetica", fontSize=7.5, leading=9.5, textColor=INK),
    "tile_v": ParagraphStyle("tv", fontName="Helvetica-Bold", fontSize=15, leading=18, textColor=INK),
    "tile_l": ParagraphStyle("tl", fontName="Helvetica", fontSize=7.5, leading=9, textColor=MUTED),
}


def clean(html: str) -> str:
    html = html or ""
    html = re.sub(r"<code>(.*?)</code>", r"<font face='Courier'>\1</font>", html)
    html = html.replace("<i>", "<i>").replace("≥", "&gt;=").replace("≤", "&lt;=")
    html = re.sub(r"<(?!/?(b|i|font)\b)[^>]+>", "", html)
    return html.replace("−", "-").replace("γδ", "gd").replace("×", "x").replace("≈", "~").replace("²", "2")


def build_pdf(root: Path):
    res = json.loads((root / "result.json").read_text())
    W, H = A4
    margin = 16 * mm
    width = W - 2 * margin
    story = []

    story += [Paragraph("TRANSCRIPT LENS · " + ("SINGLE-CELL RNA-SEQ" if res["kind"] == "sc" else "BULK RNA-SEQ") + " REPORT", S["kicker"]),
              Paragraph(clean(res["name"]), S["title"]), Spacer(1, 3),
              Paragraph(datetime.now().strftime("Generated %d %B %Y, %H:%M"), S["muted"]), Spacer(1, 10)]

    tiles = res["tiles"]
    cells = [[Paragraph(clean(str(t["value"])), S["tile_v"]), Paragraph(clean(t["label"]), S["tile_l"])] for t in tiles]
    row = [[Table([[c[0]], [c[1]]], style=[("LEFTPADDING", (0, 0), (-1, -1), 0)]) for c in cells[i:i + 3]] for i in range(0, len(cells), 3)]
    if row:
        for r in row:
            while len(r) < 3:
                r.append("")
        t = Table(row, colWidths=[width / 3] * 3)
        t.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.6, LINE), ("INNERGRID", (0, 0), (-1, -1), 0.6, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 7), ("BOTTOMPADDING", (0, 0), (-1, -1), 7), ("LEFTPADDING", (0, 0), (-1, -1), 9),
                               ("VALIGN", (0, 0), (-1, -1), "TOP")]))
        story += [t, Spacer(1, 12)]

    story.append(Paragraph("Key findings", S["h1"]))
    fl = []
    for f in res["flags"]:
        tag = {"ok": ("OK", ACCENT), "warn": ("CHECK", WARN), "info": ("NOTE", MUTED)}[f["level"]]
        fl.append([Paragraph(f"<font color='#{tag[1].hexval()[2:]}'><b>{tag[0]}</b></font>", S["cell"]), Paragraph(clean(f["text"]), S["body"])])
    if fl:
        t = Table(fl, colWidths=[16 * mm, width - 16 * mm])
        t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LINEBELOW", (0, 0), (-1, -2), 0.4, LINE),
                               ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5), ("LEFTPADDING", (0, 0), (-1, -1), 0)]))
        story.append(t)

    for sec in res["sections"]:
        story += [PageBreak(), Paragraph(clean(sec["kicker"]).upper(), S["kicker"]), Paragraph(clean(sec["title"]), S["h1"]),
                  Paragraph(clean(sec["lede"]), S["muted"]), Spacer(1, 8)]
        for it in sec["items"]:
            if it["type"] == "figure":
                p = root / it["src"]
                if not p.exists():
                    continue
                with Image.open(p) as im:
                    w, h = im.size
                maxw = width if it.get("wide") else width * 0.72
                dw = min(maxw, w / 200 * 72 * 0.95)
                dh = dw * h / w
                if dh > 150 * mm:
                    dh = 150 * mm; dw = dh * w / h
                block = [Paragraph(clean(it["title"]), S["h2"])]
                if it.get("sub"):
                    block.append(Paragraph(clean(it["sub"]), S["sub"]))
                block += [Spacer(1, 3), RLImage(str(p), width=dw, height=dh, hAlign="LEFT"), Spacer(1, 4)]
                if it.get("how"):
                    block.append(Paragraph("<b>How to read it.</b> " + clean(it["how"]), S["muted"]))
                if it.get("yours"):
                    block.append(Paragraph("<b>In your data.</b> " + clean(it["yours"]), S["body"]))
                block.append(Spacer(1, 12))
                story.append(KeepTogether(block))
            elif it["type"] == "table":
                cols = it["columns"]
                data = [[Paragraph(f"<b>{clean(c)}</b>", S["cell"]) for c in cols]] + [[Paragraph(clean(str(v)), S["cell"]) for v in r] for r in it["rows"][:60]]
                t = Table(data, repeatRows=1, colWidths=None)
                t.setStyle(TableStyle([("LINEBELOW", (0, 0), (-1, 0), 0.8, INK), ("LINEBELOW", (0, 1), (-1, -1), 0.3, LINE),
                                       ("TOPPADDING", (0, 0), (-1, -1), 2.5), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5), ("VALIGN", (0, 0), (-1, -1), "TOP")]))
                story += [Paragraph(clean(it["title"]), S["h2"]), t]
                if it.get("note"):
                    story.append(Paragraph(clean(it["note"]), S["sub"]))
                story.append(Spacer(1, 12))
            elif it["type"] == "explorer":
                pass

    story += [PageBreak(), Paragraph("METHODS", S["kicker"]), Paragraph("How this analysis was computed", S["h1"])]
    for h, p in res["methods"]:
        story += [Paragraph(clean(h), S["h2"]), Paragraph(clean(p), S["body"])]
    if res.get("versions"):
        story += [Spacer(1, 8), Paragraph("Software versions", S["h2"]),
                  Paragraph(", ".join(f"{k} {v}" for k, v in res["versions"].items()), S["muted"])]
    story += [Spacer(1, 10), Paragraph("Interpretations are generated automatically from the numbers above and are meant to guide, not replace, expert review.", S["sub"])]

    def deco(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica", 7.5); canvas.setFillColor(MUTED)
        canvas.drawString(margin, 10 * mm, "RNAseq Bench")
        canvas.drawRightString(W - margin, 10 * mm, f"{doc.page}")
        canvas.setStrokeColor(ACCENT); canvas.setLineWidth(2); canvas.line(margin, H - 10 * mm, margin + 18 * mm, H - 10 * mm)
        canvas.restoreState()

    doc = SimpleDocTemplate(str(root / "report.pdf"), pagesize=A4, leftMargin=margin, rightMargin=margin, topMargin=16 * mm, bottomMargin=16 * mm,
                            title=f"RNAseq Bench report — {res['name']}", author="RNAseq Bench")
    doc.build(story, onFirstPage=deco, onLaterPages=deco)
