"""Generate main-text Figures 2-4 for the latent-recognition biosensor screen.

Every value plotted is derived from the processed data tables in this
repository. Inputs and outputs:

  responses / dose-response curves / fits -> data/processed/
  sensor sequences                        -> data/metadata/sensor_sequences.csv
  figures (PNG at 600 dpi + vector PDF)   -> figures/
  per-figure source tables (CSV)          -> analysis/tables/

  Figure 2  full 63-ligand x 18-sensor screen, grouped by ligand scope class
  Figure 3  sequence-function cliffs across the OpuBC-cpGFP lineage
  Figure 4  dose-response curves for the strongest latent-recognition leads

Figure 1 is a composite of rendered structural panels and is assembled by
compose_figure1.py, not by this script.

Run:  python analysis/make_figures.py
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
METADATA = ROOT / "data" / "metadata"
SEQ_CSV = METADATA / "sensor_sequences.csv"
STRUCTURES = ROOT / "data" / "structures"
FIG = ROOT / "figures"
ANALYSIS = ROOT / "analysis" / "tables"

COLOR = {
    "blue": "#2B6CB0",
    "sky": "#63B3ED",
    "teal": "#319795",
    "green": "#2F855A",
    "orange": "#DD6B20",
    "red": "#C53030",
    "purple": "#6B46C1",
    "gray": "#4A5568",
    "light": "#EDF2F7",
}

SHORT_SENSOR_LABELS = {
    "iLevaphenolSnFR1.0": "iLevorphanolSnFR1.0",
}

SENSOR_DISPLAY = {
    "v4.6": "iNicSnFR1 (v4.6)",
    "V4.8.1.2": "Nic-1 int. (V4.8.1.2)",
    "V6": "Nic-2 int. (V6)",
    "cc93": "iNicSnFR3a (cc93)",
    "L194": "iSeroSnFR int. (L194)",
    "AK1": "AK1 int. (AK1)",
    "Fent2 436L": "iFentanylSnFR neg. ctrl.",
    "V7": "iNicSnFR3b (V7)",
    "iEscSnFR": "iEscSnFR",
    "Tap1.0": "iTapentadolSnFR (Tap1.0)",
    "V7.1": "ACh-1 int. (V7.1)",
    "V7.1.2": "ACh-2 int. (V7.1.2)",
    "V8": "ACh-3 int. (V8)",
    "V9": "iAChSnFR (V9)",
    "iCytSnFR": "iCytSnFR",
    "iCytBrEtSnFR": "iCyt_BrEt_SnFR",
    "iFloxSnFR": "iFluoxSnFR (iFloxSnFR)",
    "iLevaphenolSnFR1.0": "iLevorphanolSnFR1.0",
}

TRUNK_ORDER = [
    "v4.6",
    "V4.8.1.2",
    "V6",
    "cc93",
    "V7",
    "V7.1",
    "V7.1.2",
    "V8",
    "V9",
]

BRANCH_ORDER = [
    "L194",
    "AK1",
    "Fent2 436L",
    "iEscSnFR",
    "Tap1.0",
    "iCytSnFR",
    "iCytBrEtSnFR",
    "iFloxSnFR",
    "iLevaphenolSnFR1.0",
]

LIGAND_LABELS = {
    "ciprofloxacin hydrochloride hydrate": "ciprofloxacin",
    "L-(+)-ergothioneine": "ergothioneine",
    "thiamine (hydrochloride)": "thiamine",
    "mono(2‐ethylhexyl) phthalate (mehp)": "MEHP",
    "tryptamine hydrochloride": "tryptamine",
    "Histamine dihydrochloride": "histamine",
    "valproic acid sodium salt": "valproic acid",
    "atorvastatin (calcium salt hydrate)": "atorvastatin",
    "D-(+)-glucosamine hydrochloride": "glucosamine",
    "(-)-epinepherine (+)-bitartrate salt": "epinephrine",
    "(–)-Norepinephrine (bitartrate hydrate)": "norepinephrine",
    "17 beta estradiol": "estradiol",
}


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8.5,
            "axes.labelsize": 8.5,
            "axes.titlesize": 10,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.5,
            "ytick.major.width": 0.5,
            "xtick.major.size": 2.2,
            "ytick.major.size": 2.2,
            "figure.dpi": 150,
        }
    )


def save(fig: plt.Figure, name: str) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=600, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def panel(ax: plt.Axes, label: str) -> None:
    if hasattr(ax, "text2D"):
        ax.text2D(-0.11, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=11.0)
    else:
        ax.text(-0.11, 1.04, label, transform=ax.transAxes, fontweight="bold", fontsize=11.0)


def clean_ligand(name: str) -> str:
    return LIGAND_LABELS.get(name, name).replace("‐", "-")


def load_matrix() -> pd.DataFrame:
    return pd.read_csv(DATA / "response_matrix.csv", index_col=0)


def mol_from_smiles(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None and "." in smiles:
        parts = [Chem.MolFromSmiles(x) for x in smiles.split(".")]
        parts = [m for m in parts if m is not None and m.GetNumHeavyAtoms() > 1]
        mol = max(parts, key=lambda m: m.GetNumHeavyAtoms()) if parts else None
    return mol


def classify_ligands() -> pd.DataFrame:
    smiles = pd.read_csv(METADATA / "ligand_smiles.csv")
    categories = pd.read_csv(METADATA / "ligand_categories.csv")
    df = smiles.merge(categories, on="ligand", how="left")

    alkaloid_like = {
        "Nicotine",
        "Caffeine",
        "theophylline",
        "theobromine",
        "tryptamine hydrochloride",
        "Histamine dihydrochloride",
        "betahistine",
        "melatonin",
        "ritalinic acid",
        "nicotinamide",
    }

    rows = []
    for _, row in df.iterrows():
        mol = mol_from_smiles(row["SMILES"])
        formal_charge = Chem.GetFormalCharge(mol) if mol is not None else np.nan
        n_atoms = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7) if mol is not None else 0
        has_perm_cation = any(a.GetFormalCharge() > 0 for a in mol.GetAtoms()) if mol is not None else False
        has_anion = any(a.GetFormalCharge() < 0 for a in mol.GetAtoms()) if mol is not None else False
        heavy = mol.GetNumHeavyAtoms() if mol is not None else np.nan
        mw = Descriptors.MolWt(mol) if mol is not None else np.nan

        if row["ligand"] in alkaloid_like:
            scope = "alkaloid-like"
        elif has_perm_cation or has_anion or formal_charge != 0:
            scope = "non-alkaloid charged"
        else:
            scope = "non-alkaloid uncharged"

        rows.append(
            {
                "ligand": row["ligand"],
                "display": clean_ligand(row["ligand"]),
                "category": row["category"],
                "scope": scope,
                "formal_charge": formal_charge,
                "nitrogen_atoms": n_atoms,
                "has_permanent_cation": has_perm_cation,
                "has_anion": has_anion,
                "heavy_atoms": heavy,
                "mol_wt": mw,
                "SMILES": row["SMILES"],
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(ANALYSIS / "ligand_scope_classification.csv", index=False)
    return out


def figure2_screen_scope() -> None:
    """Figure 2 - full screen heatmap, per-ligand hit breadth, scope-class structures.

    Panel A displays values <= 0 as zero and caps the colour scale at
    ΔF/F0 = 4, so moderate hits stay visible alongside the much larger
    nicotine response (max 12.1).
    """
    matrix = load_matrix()
    lig = classify_ligands()
    lig["max_response"] = lig["ligand"].map(matrix.max(axis=1))
    lig["hit_breadth_03"] = lig["ligand"].map((matrix > 0.3).sum(axis=1))
    lig["hit_breadth_10"] = lig["ligand"].map((matrix > 1.0).sum(axis=1))
    lig["any_hit_03"] = lig["max_response"] > 0.3
    lig["any_hit_10"] = lig["max_response"] > 1.0

    scope_order = ["alkaloid-like", "non-alkaloid charged", "non-alkaloid uncharged"]
    sensor_order = TRUNK_ORDER + [s for s in BRANCH_ORDER if s in matrix.columns]
    lig = lig.sort_values(
        ["scope", "category", "hit_breadth_03", "max_response", "display"],
        ascending=[True, True, False, False, True],
        key=lambda s: s.map({v: i for i, v in enumerate(scope_order)}) if s.name == "scope" else s,
    )
    ordered = matrix.loc[lig["ligand"], sensor_order]
    display_values = ordered.clip(lower=0)

    fig = plt.figure(figsize=(7.4, 8.0), constrained_layout=True)
    gs = GridSpec(1, 3, figure=fig, width_ratios=[5.0, 1.15, 1.75])
    ax = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1], sharey=ax)
    ax_structures = fig.add_subplot(gs[0, 2])

    im = ax.imshow(
        display_values.values,
        aspect="auto",
        cmap="Blues",
        vmin=0,
        vmax=4.0,
        interpolation="nearest",
    )
    ax.set_xticks(range(len(sensor_order)))
    ax.set_xticklabels([SENSOR_DISPLAY.get(s, s) for s in sensor_order], rotation=90, ha="center", fontsize=6.1)
    ax.set_yticks(range(len(lig)))
    ax.set_yticklabels(lig["display"], fontsize=6.2)
    ax.tick_params(length=0)
    ax.set_title("Full single-concentration screen grouped by ligand scope")
    panel(ax, "a")
    cbar = fig.colorbar(im, ax=ax, fraction=0.026, pad=0.01)
    cbar.set_label("ΔF/F0\n(values <= 0 shown as 0;\nscale capped at 4)")

    starts = lig.groupby("scope", sort=False).head(1).index
    scope_two_line = {
        "alkaloid-like": "alkaloid-\nlike",
        "non-alkaloid charged": "non-alkaloid\ncharged",
        "non-alkaloid uncharged": "non-alkaloid\nuncharged",
    }
    boundary_positions = []
    start = 0
    for scope in scope_order:
        n = int((lig["scope"] == scope).sum())
        if n == 0:
            continue
        mid = start + n / 2 - 0.5
        ax.text(-8.4, mid, scope_two_line.get(scope, scope), ha="right", va="center",
                fontsize=7.8, fontweight="bold", linespacing=1.0)
        boundary_positions.append(start - 0.5)
        start += n
    for y in boundary_positions[1:] + [len(lig) - 0.5]:
        ax.axhline(y, color="black", lw=1.35, zorder=5)
        ax_b.axhline(y, color="black", lw=1.0, zorder=5)

    y = np.arange(len(lig))
    ax_b.barh(y - 0.19, lig["hit_breadth_03"], color=COLOR["sky"], height=0.34, label="ΔF/F0 > 0.3")
    ax_b.barh(y + 0.19, lig["hit_breadth_10"], color=COLOR["blue"], height=0.34, label="ΔF/F0 > 1.0")
    ax_b.tick_params(axis="y", labelleft=False, length=0)
    ax_b.set_xlim(0, 18)
    ax_b.set_xlabel("Number of sensors")
    ax_b.set_title("Sensors responding\nper ligand")
    ax_b.legend(frameon=False, loc="lower right", fontsize=6.9)
    ax_b.text(-0.32, 1.10, "b", transform=ax_b.transAxes, fontweight="bold", fontsize=11.0)

    ax_structures.axis("off")
    representatives = [
        ("Nicotine", "Nicotine (alkaloid-like)"),
        ("thiamine (hydrochloride)", "Thiamine (non-alkaloid charged)"),
        ("DEHP", "DEHP (non-alkaloid uncharged)"),
    ]
    smiles_by_ligand = lig.set_index("ligand")["SMILES"]
    representative_mols = [mol_from_smiles(smiles_by_ligand.loc[k]) for k, _ in representatives]
    structure_grid = Draw.MolsToGridImage(
        representative_mols,
        molsPerRow=1,
        subImgSize=(520, 300),
        legends=[label for _, label in representatives],
        useSVG=False,
    )
    ax_structures.imshow(structure_grid)
    ax_structures.set_title("Representative\nstructures", fontsize=9.4, pad=4, linespacing=1.05)
    ax_structures.text(-0.05, 1.15, "c", transform=ax_structures.transAxes, fontweight="bold", fontsize=11.0)

    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    lig.to_csv(ANALYSIS / "figure2_ligand_scope_summary.csv", index=False)
    save(fig, "figure2_screen_scope_heatmap")


def hamming_with_terminal_tolerance(a: str, b: str) -> int:
    a = a.upper()
    b = b.upper()
    n = min(len(a), len(b))
    diffs = sum(aa != bb for aa, bb in zip(a[:n], b[:n]))
    extra = abs(len(a) - len(b))
    if extra == 1 and (a.endswith("L") or b.endswith("L")):
        extra = 0
    return diffs + extra


def sensor_sequences() -> dict[str, str]:
    seqs = pd.read_csv(SEQ_CSV)
    return {str(row["sensor"]): str(row["sequence"]).upper() for _, row in seqs.iterrows()}


def figure3_activity_cliffs() -> None:
    """Figure 3 - sensor lineage, hit-status changes per parent-child step, V7/Tap1.0.

    Fent2 436L is the sequence outgroup and is excluded from all pairwise
    sequence statistics. Panel C uses a broken y-axis so the nicotine response
    does not compress the remaining ligands.
    """
    matrix = load_matrix()
    seqs = sensor_sequences()
    sensors = [s for s in matrix.columns if s != "Fent2 436L"]

    pairs = []
    for i, a in enumerate(sensors):
        for b in sensors[i + 1 :]:
            mut = hamming_with_terminal_tolerance(seqs[a], seqs[b])
            corr = float(np.corrcoef(matrix[a], matrix[b])[0, 1])
            pairs.append({"sensor_a": a, "sensor_b": b, "mutations": mut, "response_corr": corr})
    pairs_df = pd.DataFrame(pairs)
    pairs_df.to_csv(ANALYSIS / "figure3_sequence_activity_cliff_pairs.csv", index=False)

    tap_v7 = matrix[["Tap1.0", "V7"]].copy()
    tap_v7["difference"] = tap_v7["V7"] - tap_v7["Tap1.0"]
    top = tap_v7.reindex(tap_v7["difference"].abs().sort_values(ascending=False).head(14).index)

    fig = plt.figure(figsize=(7.4, 6.0), constrained_layout=True)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.05, 1.35], width_ratios=[1.0, 1.45], wspace=0.42)
    ax_tree = fig.add_subplot(gs[0, :])
    ax_rewire = fig.add_subplot(gs[1, 0])
    fp_gs = gs[1, 1].subgridspec(2, 1, height_ratios=[0.42, 1.0], hspace=0.05)
    ax_fp_high = fig.add_subplot(fp_gs[0, 0])
    ax_fp = fig.add_subplot(fp_gs[1, 0], sharex=ax_fp_high)

    ax_tree.axis("off")
    y0 = 0.76
    # Trunk node labels sit close under their nodes (0.045 rather than 0.10) so a
    # clear band remains between the legend row above and the sensor names.
    TRUNK_LABEL_DY = 0.045
    trunk_x = np.linspace(0.06, 0.90, len(TRUNK_ORDER))
    for i, s in enumerate(TRUNK_ORDER):
        ax_tree.scatter(trunk_x[i], y0, s=112, color=COLOR["blue"], edgecolor="white", lw=0.6, zorder=3)
        ax_tree.text(trunk_x[i], y0 + TRUNK_LABEL_DY, s, ha="center", va="bottom", fontsize=8.6)
        if i:
            ax_tree.plot([trunk_x[i - 1], trunk_x[i]], [y0, y0], color="#2D3748", lw=1.5, zorder=1)
    positions = {s: (trunk_x[i], y0) for i, s in enumerate(TRUNK_ORDER)}
    # (x, node_y, parent, label_dy, label_dx) — label_dy staggers long sibling labels
    # into two tiers so adjacent names never share a horizontal band.
    # label_dx shifts a label sideways off its own descending connector. AK1 is a
    # pass-through node (the edge to Fent2 436L continues straight down through
    # x=0.37), so a centred label would sit on that line and is nudged left.
    branch_positions = {
        "L194": (0.25, 0.40, "cc93", -0.075, 0.0),
        "AK1": (0.37, 0.40, "cc93", -0.075, -0.035),
        "Fent2 436L": (0.37, 0.13, "AK1", -0.075, 0.0),
        "iEscSnFR": (0.49, 0.40, "V7", -0.075, 0.0),
        "Tap1.0": (0.60, 0.40, "V7", -0.075, 0.0),
        "iCytSnFR": (0.70, 0.36, "V9", -0.075, 0.0),
        "iCytBrEtSnFR": (0.82, 0.24, "V9", -0.075, 0.0),
        "iFloxSnFR": (0.94, 0.36, "V9", -0.075, 0.0),
        "iLevaphenolSnFR1.0": (1.06, 0.24, "V9", -0.075, 0.0),
    }
    for s, (x, y, parent, label_dy, label_dx) in branch_positions.items():
        px, py = positions[parent]
        positions[s] = (x, y)
        color = COLOR["red"] if s == "Tap1.0" else COLOR["teal"]
        joint_y = 0.55 if parent in {"cc93", "V7"} else (0.46 if parent == "V9" else 0.26)
        ax_tree.plot([px, px, x, x], [py, joint_y, joint_y, y], color="#2D3748", lw=1.35, zorder=1)
        ax_tree.scatter(x, y, s=76, color=color, edgecolor="white", lw=0.5, zorder=3)
        ha = "center" if label_dx == 0 else ("right" if label_dx < 0 else "left")
        ax_tree.text(x + label_dx, y + label_dy, SHORT_SENSOR_LABELS.get(s, s),
                     ha=ha, va="top", fontsize=7.8)
    # W436A annotates the V7 -> Tap1.0 edge; placed just right of that vertical, clear of node labels
    ax_tree.text(0.615, 0.50, "W436A", color=COLOR["red"], fontsize=9.1, fontweight="bold", ha="left", va="center")
    ax_tree.set_xlim(0, 1.14)
    ax_tree.set_ylim(0.02, 0.98)
    ax_tree.scatter([], [], s=55, color=COLOR["blue"], label="trunk")
    ax_tree.scatter([], [], s=42, color=COLOR["teal"], label="branch")
    ax_tree.scatter([], [], s=42, color=COLOR["red"], label="Tap1.0")
    ax_tree.legend(frameon=False, loc="upper left", fontsize=7.5, handletextpad=0.35, borderpad=0.1, ncol=3)
    ax_tree.set_title("Screened OpuBC/cpGFP lineage", fontsize=11.0)
    panel(ax_tree, "a")

    lineage_edges = list(zip(TRUNK_ORDER[:-1], TRUNK_ORDER[1:])) + [
        ("cc93", "L194"),
        ("cc93", "AK1"),
        ("V7", "iEscSnFR"),
        ("V7", "Tap1.0"),
        ("V9", "iCytSnFR"),
        ("V9", "iCytBrEtSnFR"),
        ("V9", "iFloxSnFR"),
        ("V9", "iLevaphenolSnFR1.0"),
    ]
    transition_rows = []
    for parent_name, child_name in lineage_edges:
        mutations = hamming_with_terminal_tolerance(seqs[parent_name], seqs[child_name])
        if mutations > 3:
            continue
        parent_hits = matrix[parent_name] > 0.3
        child_hits = matrix[child_name] > 0.3
        gained = int((~parent_hits & child_hits).sum())
        lost = int((parent_hits & ~child_hits).sum())
        transition_rows.append(
            {
                "parent": parent_name,
                "child": child_name,
                "mutations": mutations,
                "gained": gained,
                "lost": lost,
                "changed": gained + lost,
            }
        )
    transitions = pd.DataFrame(transition_rows).sort_values(
        ["mutations", "changed", "parent", "child"],
        ascending=[True, False, True, True],
    )
    short_names = {
        "iEscSnFR": "iEsc",
        "iLevaphenolSnFR1.0": "iLevorphanol",
    }
    row_labels = [
        f"{short_names.get(row.parent, row.parent)} to {short_names.get(row.child, row.child)}  ({row.mutations} mut.)"
        for row in transitions.itertuples()
    ]
    y_positions = np.arange(len(transitions))
    bar_colors = [COLOR["red"] if child == "Tap1.0" else "#718096" for child in transitions["child"]]
    ax_rewire.barh(y_positions, transitions["changed"], color=bar_colors, height=0.64)
    for y_pos, changed in zip(y_positions, transitions["changed"]):
        ax_rewire.text(changed + 0.25, y_pos, str(changed), va="center", ha="left", fontsize=8.1)
    ax_rewire.set_yticks(y_positions)
    ax_rewire.set_yticklabels(row_labels, fontsize=7.5)
    ax_rewire.invert_yaxis()
    ax_rewire.set_xlim(0, 13.5)
    ax_rewire.set_xticks([0, 3, 6, 9, 12])
    ax_rewire.set_xlabel("Ligands with changed hit status (of 63)")
    ax_rewire.set_title("1-3 mutations rewire many responses", fontsize=9.4)
    ax_rewire.grid(axis="x", color="#E2E8F0", lw=0.6)
    ax_rewire.set_axisbelow(True)
    ax_rewire.text(-0.11, 1.14, "b", transform=ax_rewire.transAxes, fontweight="bold", fontsize=11.0)

    x = np.arange(len(top))
    width = 0.38
    for fp_axis in [ax_fp, ax_fp_high]:
        fp_axis.bar(x - width / 2, top["Tap1.0"], width=width, color=COLOR["orange"], label="Tap1.0")
        fp_axis.bar(x + width / 2, top["V7"], width=width, color=COLOR["blue"], label="V7")
    ax_fp.axhline(0.3, color=COLOR["gray"], lw=0.7, ls=":")
    ax_fp.axhline(1.0, color=COLOR["gray"], lw=0.7, ls="--")
    ax_fp.set_xticks(x)
    ax_fp.set_xticklabels([clean_ligand(v) for v in top.index], rotation=48, ha="right", fontsize=7.5)
    ax_fp.set_ylabel("ΔF/F0")
    ax_fp.set_ylim(-1.05, 3.15)
    ax_fp_high.set_ylim(9.6, 11.35)
    ax_fp_high.set_title("V7 vs Tap1.0: one mutation (W436A)\nswitches responses from near zero to >1", fontsize=9.0, pad=6)
    ax_fp_high.legend(frameon=False, ncol=2, loc="upper right", fontsize=7.0)
    ax_fp_high.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_fp_high.spines["bottom"].set_visible(False)
    ax_fp.spines["top"].set_visible(False)
    break_kwargs = dict(color="#2D3748", clip_on=False, lw=0.8)
    ax_fp_high.plot((-0.012, 0.012), (-0.02, 0.02), transform=ax_fp_high.transAxes, **break_kwargs)
    ax_fp_high.plot((0.988, 1.012), (-0.02, 0.02), transform=ax_fp_high.transAxes, **break_kwargs)
    ax_fp.plot((-0.012, 0.012), (0.98, 1.02), transform=ax_fp.transAxes, **break_kwargs)
    ax_fp.plot((0.988, 1.012), (0.98, 1.02), transform=ax_fp.transAxes, **break_kwargs)
    # Threshold key sits in the right margin over the low-response ligands, clear
    # of the long DEHP and thiamine bars.
    ax_fp.text(0.985, 0.97, "dotted: 0.3; dashed: 1.0", transform=ax_fp.transAxes,
               fontsize=7.2, va="top", ha="right")
    ax_fp_high.text(-0.13, 1.30, "c", transform=ax_fp_high.transAxes, fontweight="bold", fontsize=11.0)

    for a in [ax_rewire, ax_fp, ax_fp_high]:
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    ax_fp_high.spines["top"].set_visible(False)
    save(fig, "figure3_sequence_activity_cliffs")


def hill(x, baseline, df_max, ec50, n):
    x = np.asarray(x, dtype=float)
    return baseline + df_max / (1.0 + (ec50 / x) ** n)


def read_curve(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    conc_col = next((cols[c] for c in cols if "conc" in c or "um" in c or "µm" in c), df.columns[0])
    mean_col = next((c for c in df.columns if "mean" in c.lower() or "df" in c.lower()), df.columns[1])
    sem_col = next((c for c in df.columns if "sem" in c.lower()), None)
    out = pd.DataFrame({"conc_uM": df[conc_col], "mean": df[mean_col]})
    out["sem"] = df[sem_col].abs() if sem_col else 0.0
    return out.sort_values("conc_uM")


def figure4_dose_response_leads() -> None:
    """Figure 4 - Hill fits for the strongest leads across nine analyte-sensor curves.

    Solid curves are within-range fits (EC50 < 400 uM, dFmax < 20); dashed grey
    curves mark dose responses that never plateau in the measured range, whose
    fitted EC50 values are lower bounds rather than potency estimates.
    """
    fits = pd.read_csv(DATA / "dose_response_fits.csv")
    files = [
        "cipro_iEscSnFR.csv",
        "cipro_v7.1.csv",
        "cipro_v7.csv",
        "cipro_v9.csv",
        "dehp_v4.8.1.2.csv",
        "dehp_L194.csv",
        "thiamine_iEscSnFR.csv",
        "ergothioneine_v9.csv",
        "carnitine_v9.csv",
    ]
    titles = {
        "cipro_iEscSnFR.csv": "ciprofloxacin, iEscSnFR",
        "cipro_v7.1.csv": "ciprofloxacin, V7.1",
        "cipro_v7.csv": "ciprofloxacin, V7",
        "cipro_v9.csv": "ciprofloxacin, V9",
        "dehp_v4.8.1.2.csv": "DEHP, V4.8.1.2",
        "dehp_L194.csv": "DEHP, L194",
        "thiamine_iEscSnFR.csv": "thiamine, iEscSnFR",
        "ergothioneine_v9.csv": "ergothioneine, V9",
        "carnitine_v9.csv": "L-carnitine, V9",
    }
    fig, axes = plt.subplots(3, 3, figsize=(7.1, 6.4), constrained_layout=True)
    for ax, fname in zip(axes.flat, files):
        curve = read_curve(DATA / "dose_responses" / fname)
        row = fits[fits["file"] == fname].iloc[0]
        ax.errorbar(curve["conc_uM"], curve["mean"], yerr=curve["sem"], fmt="o", ms=3.2, color=COLOR["blue"], ecolor="#718096", capsize=1.8)
        xgrid = np.logspace(math.log10(max(curve["conc_uM"].min(), 1e-4)), math.log10(curve["conc_uM"].max() * 1.25), 300)
        yfit = hill(xgrid, row["baseline"], row["df_max"], row["ec50_uM"], row["hill_n"])
        reliable = row["ec50_uM"] < 400 and row["df_max"] < 20
        ax.plot(
            xgrid,
            yfit,
            color=COLOR["red"] if reliable else COLOR["gray"],
            lw=1.2,
            ls="-" if reliable else "--",
        )
        ax.set_xscale("log")
        ax.set_title(titles[fname], fontsize=9.4)
        if reliable:
            note = f"EC50={row['ec50_uM']:.1f} µM"
        else:
            note = "lower-bound fit"
        ax.text(0.04, 0.92, note, transform=ax.transAxes, fontsize=7.5, va="top")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    for ax in axes[-1, :]:
        ax.set_xlabel("Ligand concentration (µM)")
    for ax in axes[:, 0]:
        ax.set_ylabel("ΔF/F0")
    panel(axes[0, 0], "a")
    panel(axes[0, 1], "b")
    panel(axes[0, 2], "c")
    save(fig, "figure4_dose_response_leads")

    summary = fits[fits["file"].isin(files)].copy()
    summary.to_csv(ANALYSIS / "figure4_dose_response_lead_summary.csv", index=False)


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    style()
    figure2_screen_scope()
    figure3_activity_cliffs()
    figure4_dose_response_leads()


if __name__ == "__main__":
    main()
