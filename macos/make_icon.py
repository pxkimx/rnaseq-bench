#!/usr/bin/env python3
"""Draw the app icon from the same geometry as the logo in web/index.html, and build AppIcon.icns.

    python macos/make_icon.py

Kept as code rather than a checked-in binary so the icon and the in-app mark cannot drift apart: the
coordinates below are the SVG's 32-unit viewBox, scaled onto Apple's 1024 canvas.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image, ImageDraw

HERE = Path(__file__).resolve().parent
RES = HERE / "RNAseq Bench.app" / "Contents" / "Resources"

SS = 4                      # supersampling; PIL has no antialiased shapes, so draw big and shrink
CANVAS = 1024
INSET, BOX = 100, 824       # Apple leaves a margin around the rounded square
PANEL = "#0d1b22"           # --ink in the light theme: reads on both a light and a dark dock
AXES = "#6a7b84"            # --muted

# (cx, cy, r, colour) in the SVG's 32-unit space
CLUSTERS = [
    (13.0, 11.0, 3.4, "#00e0cf"), (21.2, 17.0, 4.4, "#ffb020"),
    (14.6, 20.4, 2.4, "#ff5c8a"), (22.2, 9.0, 1.6, "#4ea3ff"),
    (24.6, 21.8, 1.3, "#b27bff"), (10.0, 17.4, 1.2, "#a3e635"),
]


def draw(size: int) -> Image.Image:
    s = size * SS
    k = (BOX / 32) * (s / CANVAS)           # 32 viewBox units -> pixels
    off = INSET * (s / CANVAS)
    to = lambda v: off + v * k              # noqa: E731

    im = Image.new("RGBA", (s, s), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    d.rounded_rectangle([to(0.5), to(0.5), to(31.5), to(31.5)], radius=7.5 * k, fill=PANEL)

    w = max(1, int(1.05 * k))
    cap = w / 2
    for (x0, y0, x1, y1) in [(7, 4.6, 7, 25.4), (6.6, 25, 27.4, 25)]:
        d.line([to(x0), to(y0), to(x1), to(y1)], fill=AXES, width=w)
        for (cx, cy) in ((x0, y0), (x1, y1)):      # round the ends, as stroke-linecap does
            d.ellipse([to(cx) - cap, to(cy) - cap, to(cx) + cap, to(cy) + cap], fill=AXES)

    for cx, cy, r, col in CLUSTERS:
        d.ellipse([to(cx - r), to(cy - r), to(cx + r), to(cy + r)], fill=col)
    return im.resize((size, size), Image.LANCZOS)


def main() -> int:
    if not RES.is_dir():
        print(f"no Resources folder at {RES}", file=sys.stderr)
        return 1
    draw(512).save(RES / "AppIcon.png")

    iconset = HERE / "AppIcon.iconset"
    shutil.rmtree(iconset, ignore_errors=True)
    iconset.mkdir()
    for px in (16, 32, 64, 128, 256, 512, 1024):
        draw(px).save(iconset / f"icon_{px}x{px}.png")
        if px > 16:                                  # @2x is the next size up, named as the half
            draw(px).save(iconset / f"icon_{px // 2}x{px // 2}@2x.png")
    r = subprocess.run(["iconutil", "-c", "icns", str(iconset), "-o", str(RES / "AppIcon.icns")])
    shutil.rmtree(iconset, ignore_errors=True)
    if r.returncode:
        return r.returncode
    print(f"wrote {RES/'AppIcon.icns'} and AppIcon.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
