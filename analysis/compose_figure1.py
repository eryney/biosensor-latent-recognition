"""Compose the four Figure 1 panels into a single figure.

Layout: (a) scaffold architecture top-left, (c) docked ciprofloxacin top-right,
(b) sensing-mechanism schematic bottom-left, (d) prediction-vs-crystal overlay
bottom-right. Panels are flattened onto white and autocropped before placement.

Reads the four panel PNGs from figures/panels/ and writes
figures/figure1_scaffold_structure.{png,pdf} at 300 dpi.

Run:  python analysis/compose_figure1.py
"""
import os, numpy as np
from PIL import Image, ImageDraw, ImageFont

FIGDIR = os.path.join(os.path.dirname(__file__), "..", "figures")
PANELDIR = os.path.join(FIGDIR, "panels")

def flatten_white(im):
    im = im.convert("RGBA")
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
    return Image.alpha_composite(bg, im).convert("RGB")

def autocrop(im, pad=20, bg=(255, 255, 255)):
    arr = np.asarray(im.convert("RGB"))
    mask = np.any(np.abs(arr.astype(int) - np.array(bg)) > 12, axis=2)
    if not mask.any():
        return im
    ys, xs = np.where(mask)
    y0, y1, x0, x1 = ys.min(), ys.max(), xs.min(), xs.max()
    return im.crop((max(0, x0 - pad), max(0, y0 - pad),
                    min(arr.shape[1], x1 + pad), min(arr.shape[0], y1 + pad)))

def load_font(sz):
    for p in ["/System/Library/Fonts/Helvetica.ttc",
              "/System/Library/Fonts/Supplemental/Arial Bold.ttf"]:
        if os.path.exists(p):
            try: return ImageFont.truetype(p, sz)
            except Exception: pass
    return ImageFont.load_default()

pa = autocrop(flatten_white(Image.open(f"{PANELDIR}/panel_a_7S7U.png")))
pc = autocrop(flatten_white(Image.open(f"{PANELDIR}/panel_c_7S7U_reference.png")))
pd_ = autocrop(flatten_white(Image.open(f"{PANELDIR}/panel_d_7S7U.png")))

F_letter, F_title, F_call = load_font(58), load_font(38), load_font(34)
CELL_W, CELL_H, PAD, TOP, GRID_GAP = 1500, 1080, 30, 70, 20

def fit(im, mw, mh):
    r = min(mw / im.width, mh / im.height)
    return im.resize((int(im.width * r), int(im.height * r)), Image.LANCZOS)

canvas = Image.new("RGB", (2 * CELL_W + 3 * GRID_GAP, 2 * CELL_H + 3 * GRID_GAP), "white")
draw = ImageDraw.Draw(canvas)

def place(im, letter, title, col, row, callout=None):
    cx = GRID_GAP + col * (CELL_W + GRID_GAP)
    cy = GRID_GAP + row * (CELL_H + GRID_GAP)
    draw.text((cx + 10, cy + 8), letter, fill="black", font=F_letter)
    if title:
        draw.text((cx + 80, cy + 22), title, fill="black", font=F_title)
    fitted = fit(im, CELL_W - 2 * PAD, CELL_H - TOP - 2 * PAD)
    ox = cx + (CELL_W - fitted.width) // 2
    oy = cy + TOP + (CELL_H - TOP - fitted.height) // 2
    canvas.paste(fitted, (ox, oy))
    if callout:
        tx, ty = ox + fitted.width - 20, oy + fitted.height - 30
        for dx in (-2, -1, 0, 1, 2):
            for dy in (-2, -1, 0, 1, 2):
                draw.text((tx + dx, ty + dy), callout, fill="white", font=F_call, anchor="rs")
        draw.text((tx, ty), callout, fill="black", font=F_call, anchor="rs")

place(pa, "a", "Scaffold architecture", 0, 0)
place(pc, "c", "Ciprofloxacin docked in pocket", 1, 0)
pb = autocrop(flatten_white(Image.open(f"{PANELDIR}/panel_b_mechanism.png")))
place(pb, "b", "Sensing mechanism", 0, 1)
place(pd_, "d", "ESMFold prediction vs crystal", 1, 1, callout="PBP RMSD 1.7 \u00c5")

canvas.save(f"{FIGDIR}/figure1_scaffold_structure.png", dpi=(300, 300))
canvas.save(f"{FIGDIR}/figure1_scaffold_structure.pdf", "PDF", resolution=300)
print("wrote figure1_scaffold_structure.{png,pdf}", canvas.size)
