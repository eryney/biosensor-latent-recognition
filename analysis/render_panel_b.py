"""Figure 1 panel (b) - schematic of the sensing mechanism.

Two states drawn side by side. In the unbound state the PBP lobes are open and
a linker glutamate rests against the cpGFP chromophore, quenching it. Ligand
binding closes the lobes, which withdraws the glutamate and relieves the
quenching, so the fluorophore brightens.

Colours match panel (a): PBP magenta, cpGFP green, connecting linker cyan.

Run:  python analysis/render_panel_b.py
Writes ../figures/panels/panel_b_mechanism.png
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

PANELDIR = os.path.join(os.path.dirname(__file__), "..", "figures", "panels")
OUT = os.path.join(PANELDIR, "panel_b_mechanism.png")

PBP = "#c020c0"
CPGFP_DIM = "#9ec49e"
CPGFP_BRIGHT = "#12b012"
LINKER = "#17b8c8"
LIGAND = "#f08000"
GLU = "#1c3a5e"
CHROMO_DIM = "#5a6b5a"
CHROMO_BRIGHT = "#eaff70"

LOBE_W, LOBE_H = 0.28, 0.24
BARREL_W, BARREL_H = 0.34, 0.36
LINK_W = 0.09
SPAN = LOBE_W + LINK_W + BARREL_W
MID = 0.52


def state(ax, x0, bound):
    """Draw one state; x0 is the left edge of the PBP lobes."""
    # Open lobes are splayed apart; ligand binding closes them.
    half = 0.050 if bound else 0.115
    for sign in (+1, -1):
        y = MID + half if sign > 0 else MID - half - LOBE_H
        ax.add_patch(FancyBboxPatch(
            (x0, y), LOBE_W, LOBE_H,
            boxstyle="round,pad=0.010,rounding_size=0.042",
            linewidth=0, facecolor=PBP, zorder=3))

    lx = x0 + LOBE_W
    ax.add_patch(FancyBboxPatch(
        (lx, MID - 0.028), LINK_W, 0.056,
        boxstyle="square,pad=0.0", linewidth=0, facecolor=LINKER, zorder=2))

    bx = lx + LINK_W
    bcx = bx + BARREL_W / 2
    if bound:
        for r, a in ((0.31, 0.08), (0.26, 0.12), (0.21, 0.17)):
            ax.add_patch(Circle((bcx, MID), r, facecolor=CPGFP_BRIGHT,
                                alpha=a, linewidth=0, zorder=1))
    ax.add_patch(FancyBboxPatch(
        (bx, MID - BARREL_H / 2), BARREL_W, BARREL_H,
        boxstyle="round,pad=0.010,rounding_size=0.055",
        linewidth=0, facecolor=CPGFP_BRIGHT if bound else CPGFP_DIM, zorder=3))

    # Chromophore, drawn inside the barrel so the quenching is legible.
    ax.add_patch(Circle((bcx, MID), 0.052,
                        facecolor=CHROMO_BRIGHT if bound else CHROMO_DIM,
                        edgecolor="white", linewidth=1.4, zorder=5))

    # Linker glutamate: against the chromophore when unbound, withdrawn into
    # the linker when the lobes close.
    gx = bcx - 0.098 if not bound else lx + 0.028
    ax.add_patch(Circle((gx, MID), 0.030, facecolor=GLU,
                        edgecolor="white", linewidth=1.2, zorder=6))

    # Labels sit outside the shapes so nothing is occluded.
    ax.text(x0 + LOBE_W / 2, MID - half - LOBE_H - 0.085, "PBP",
            ha="center", va="top", color=PBP, fontsize=13, fontweight="bold")
    ax.text(bcx, MID - BARREL_H / 2 - 0.055, "cpGFP", ha="center", va="top",
            color=CPGFP_BRIGHT if bound else "#6f8f6f", fontsize=13,
            fontweight="bold")
    ax.annotate("Glu", xy=(gx, MID + 0.030), xytext=(gx, MID + 0.255),
                ha="center", va="bottom", fontsize=11, color=GLU,
                style="italic", zorder=7,
                arrowprops=dict(arrowstyle="-", color=GLU, linewidth=1.0,
                                shrinkA=1, shrinkB=1))

    if bound:
        ax.add_patch(Circle((x0 + LOBE_W + 0.004, MID), 0.042,
                            facecolor=LIGAND, edgecolor="#a05000",
                            linewidth=1.4, zorder=7))
    else:
        ax.add_patch(Circle((x0 + LOBE_W / 2, MID + 0.375), 0.038,
                            facecolor=LIGAND, edgecolor="#a05000",
                            linewidth=1.4, zorder=7))
        ax.text(x0 + LOBE_W / 2 - 0.075, MID + 0.375, "ligand", ha="right",
                va="center", fontsize=11, color="#a05000", style="italic")


def main():
    os.makedirs(PANELDIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.6, 3.1), dpi=600)

    left_x, right_x = 0.10, 1.26
    ax.set_xlim(0.0, right_x + SPAN + 0.10)
    ax.set_ylim(0.00, 1.16)
    ax.set_aspect("equal")
    ax.axis("off")

    state(ax, left_x, bound=False)
    state(ax, right_x, bound=True)

    ax.text(left_x + SPAN / 2, 1.10, "Unbound \u2014 dim",
            ha="center", va="center", fontsize=14)
    ax.text(right_x + SPAN / 2, 1.10, "Ligand-bound \u2014 bright",
            ha="center", va="center", fontsize=14)

    ax.add_patch(FancyArrowPatch((left_x + SPAN + 0.045, MID),
                                 (right_x - 0.045, MID),
                                 arrowstyle="<|-|>", mutation_scale=16,
                                 linewidth=1.8, color="#404040",
                                 shrinkA=0, shrinkB=0))

    fig.savefig(OUT, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", os.path.normpath(OUT))


if __name__ == "__main__":
    main()
