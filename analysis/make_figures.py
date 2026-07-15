"""Generate all main-text figures from the processed data tables."""

from __future__ import annotations

import math
import json
import shlex
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patches
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from scipy.optimize import curve_fit
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
METADATA = ROOT / "data" / "metadata"
STRUCTURES = ROOT / "data" / "structures"
FIG = ROOT / "figures"
ANALYSIS = ROOT / "analysis"

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


def figure1_scaffold_prediction() -> None:
    """Intro schematic for the OpuBC/cpGFP scaffold and current structure runs."""

    fig = plt.figure(figsize=(7.1, 4.25))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.05, 1.55], height_ratios=[1.0, 0.85], hspace=0.85, wspace=0.38)
    ax_mech = fig.add_subplot(gs[0, 0])
    ax_pred = fig.add_subplot(gs[1, 0])
    ax_tree = fig.add_subplot(gs[:, 1])

    ax_mech.axis("off")
    ax_mech.set_xlim(0, 1)
    ax_mech.set_ylim(0, 1)
    for center, angle in [((0.34, 0.55), 22), ((0.66, 0.55), -22)]:
        ax_mech.add_patch(
            patches.Ellipse(center, 0.28, 0.58, angle=angle, facecolor=COLOR["sky"], edgecolor=COLOR["blue"], lw=1.1)
        )
    ax_mech.add_patch(patches.Ellipse((0.50, 0.24), 0.32, 0.22, facecolor="#68D391", edgecolor=COLOR["green"], lw=1.1))
    ax_mech.add_patch(patches.Circle((0.50, 0.55), 0.045, facecolor=COLOR["orange"], edgecolor="#7B341E", lw=0.8))
    ax_mech.plot([0.42, 0.47], [0.35, 0.29], color=COLOR["gray"], lw=1.0)
    ax_mech.plot([0.58, 0.53], [0.35, 0.29], color=COLOR["gray"], lw=1.0)
    ax_mech.annotate("", xy=(0.78, 0.56), xytext=(0.70, 0.56), arrowprops=dict(arrowstyle="->", lw=0.9, color=COLOR["gray"]))
    ax_mech.set_title("OpuBC/cpGFP sensor mechanism", pad=4)
    ax_mech.text(0.50, 0.62, "ligand", ha="center", va="bottom", fontsize=6.2)
    ax_mech.text(0.50, 0.24, "cpGFP", ha="center", va="center", fontsize=6.4, color="#1C4532")
    ax_mech.text(0.50, 0.05, "Binding closes the PBP and shifts fluorescence.", ha="center", fontsize=6.2)
    panel(ax_mech, "a")

    ax_tree.axis("off")
    ax_tree.set_xlim(0, 1)
    ax_tree.set_ylim(0, 1)
    ax_tree.set_title("Screened OpuBC/cpGFP family", pad=4)
    trunk = ["v4.6", "V4.8", "V6", "cc93", "V7", "V7.1", "V7.1.2", "V8", "V9"]
    xs = np.linspace(0.07, 0.93, len(trunk))
    y = 0.78
    for i, sensor in enumerate(trunk):
        color = COLOR["red"] if sensor == "cc93" else COLOR["blue"]
        ax_tree.scatter(xs[i], y, s=58, color=color, edgecolor="white", lw=0.5, zorder=3)
        ax_tree.text(xs[i], y - 0.055, sensor, ha="center", va="top", fontsize=5.7, rotation=35)
        if i:
            ax_tree.plot([xs[i - 1], xs[i]], [y, y], color=COLOR["gray"], lw=0.9)
    branch_rows = [
        ("L194", "cc93", 0.58, "W436T"),
        ("AK1", "cc93", 0.44, "Y357G"),
        ("Tap1.0", "V7", 0.44, "W436A"),
        ("iEscSnFR", "V7", 0.58, ""),
        ("iCytSnFR", "V9", 0.60, ""),
        ("iFloxSnFR", "V9", 0.46, "W436S"),
        ("iLevorphanol", "V9", 0.32, ""),
    ]
    parent_map = {"cc93": xs[3], "V7": xs[4], "V9": xs[8]}
    for sensor, parent, by, label in branch_rows:
        px = parent_map[parent]
        bx = px + (0.06 if parent != "V9" else -0.10)
        color = COLOR["orange"] if label else COLOR["teal"]
        ax_tree.plot([px, bx], [y, by], color=COLOR["gray"], lw=0.75)
        ax_tree.scatter(bx, by, s=50, color=color, edgecolor="white", lw=0.4, zorder=3)
        ax_tree.text(bx + (0.015 if parent != "V9" else -0.015), by, sensor, ha="left" if parent != "V9" else "right", va="center", fontsize=5.8)
        if label:
            ax_tree.text(bx + (0.015 if parent != "V9" else -0.015), by - 0.047, label, ha="left" if parent != "V9" else "right", va="center", fontsize=5.3, color=COLOR["orange"])
    ax_tree.text(0.04, 0.16, "cc93/iNicSnFR3a was chosen for structure work because it is in the screen and has PDB 7S7U.", fontsize=6.5)
    ax_tree.text(0.04, 0.08, "Y357 and W436 mark recurrent pocket changes used later to interpret activity cliffs.", fontsize=6.5)
    panel(ax_tree, "b")

    ax_pred.axis("off")
    ax_pred.text(0.08, 0.98, "Structure prediction run status", fontsize=8, va="top")
    manifest = STRUCTURES / "predictions" / "ee88e2c3-manifest.json"
    ptms = []
    aggregates = []
    if manifest.exists():
        data = json.loads(manifest.read_text())
        for item in data.get("chai1", []):
            scores = item.get("score_summary", {})
            ptms.append(scores.get("ptm", {}).get("mean", np.nan))
            aggregates.append(scores.get("aggregate_score", {}).get("mean", np.nan))
    rows = [
        ("Chai-1", f"5 models; mean pTM {np.nanmean(ptms):.2f}; mean aggregate {np.nanmean(aggregates):.3f}" if ptms else "completed"),
        ("ColabFold/AF2", "attempted; container failed at JAX/Haiku import"),
        ("RoseTTAFold", "not completed in this pass"),
    ]
    y0 = 0.66
    for i, (name, status) in enumerate(rows):
        yy = y0 - i * 0.25
        ax_pred.text(0.02, yy, name, fontsize=6.6, fontweight="bold", va="top")
        ax_pred.text(0.43, yy, status, fontsize=6.2, va="top")
        ax_pred.plot([0.02, 0.98], [yy - 0.065, yy - 0.065], color="#E2E8F0", lw=0.7)
    ax_pred.text(0.02, 0.02, "The Chai confidence is modest, so the model is treated as orientation rather than mechanism.", fontsize=6.1)
    panel(ax_pred, "c")

    save(fig, "figure1_scaffold_structure_prediction")


def parse_pdb_ca(path: Path, chain: str = "A") -> dict[int, np.ndarray]:
    coords: dict[int, np.ndarray] = {}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        if line[21].strip() != chain:
            continue
        try:
            resid = int(line[22:26])
            coords[resid] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
        except ValueError:
            continue
    return coords


def parse_cif_ca(path: Path) -> dict[int, np.ndarray]:
    lines = path.read_text(errors="ignore").splitlines()
    coords: dict[int, np.ndarray] = {}
    i = 0
    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue
        fields: list[str] = []
        i += 1
        while i < len(lines) and lines[i].startswith("_"):
            fields.append(lines[i].strip())
            i += 1
        if "_atom_site.label_atom_id" not in fields:
            continue
        idx = {field: fields.index(field) for field in fields}
        required = [
            "_atom_site.label_atom_id",
            "_atom_site.label_seq_id",
            "_atom_site.Cartn_x",
            "_atom_site.Cartn_y",
            "_atom_site.Cartn_z",
        ]
        if not all(field in idx for field in required):
            continue
        while i < len(lines) and lines[i].strip() and not lines[i].startswith(("loop_", "_", "#")):
            parts = shlex.split(lines[i])
            if len(parts) >= len(fields) and parts[idx["_atom_site.label_atom_id"]] == "CA":
                try:
                    resid = int(parts[idx["_atom_site.label_seq_id"]])
                    coords[resid] = np.array(
                        [
                            float(parts[idx["_atom_site.Cartn_x"]]),
                            float(parts[idx["_atom_site.Cartn_y"]]),
                            float(parts[idx["_atom_site.Cartn_z"]]),
                        ],
                        dtype=float,
                    )
                except ValueError:
                    pass
            i += 1
        break
    return coords


def ordered_coords(coords: dict[int, np.ndarray]) -> tuple[np.ndarray, list[int]]:
    residues = sorted(coords)
    return np.vstack([coords[resid] for resid in residues]), residues


def align_to_reference(mobile: dict[int, np.ndarray], reference: dict[int, np.ndarray]) -> tuple[dict[int, np.ndarray], float, int]:
    common = sorted(set(mobile) & set(reference))
    if len(common) < 3:
        raise ValueError("Need at least 3 shared C-alpha atoms for structural alignment")
    p = np.vstack([mobile[resid] for resid in common])
    q = np.vstack([reference[resid] for resid in common])
    p_mean = p.mean(axis=0)
    q_mean = q.mean(axis=0)
    p0 = p - p_mean
    q0 = q - q_mean
    u, _, vt = np.linalg.svd(p0.T @ q0)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    aligned = {resid: (coord - p_mean) @ r + q_mean for resid, coord in mobile.items()}
    aligned_common = np.vstack([aligned[resid] for resid in common])
    rmsd = float(np.sqrt(np.mean(np.sum((aligned_common - q) ** 2, axis=1))))
    return aligned, rmsd, len(common)


def plot_trace(ax: plt.Axes, coords: dict[int, np.ndarray], color: str, title: str, subtitle: str = "") -> None:
    arr, _ = ordered_coords(coords)
    ax.plot(arr[:, 0], arr[:, 1], arr[:, 2], color=color, lw=1.15)
    ax.scatter(arr[0, 0], arr[0, 1], arr[0, 2], s=9, color=color, edgecolor="white", linewidth=0.3)
    ax.set_title(title, pad=1, fontsize=7.2)
    if subtitle:
        ax.text2D(0.5, 0.02, subtitle, transform=ax.transAxes, ha="center", fontsize=5.8, color=COLOR["gray"])
    ax.view_init(elev=19, azim=-64)
    ax.set_axis_off()


def set_3d_bounds(axes: list[plt.Axes], coord_sets: list[dict[int, np.ndarray]]) -> None:
    all_coords = np.vstack([ordered_coords(coords)[0] for coords in coord_sets])
    mins = all_coords.min(axis=0)
    maxs = all_coords.max(axis=0)
    center = (mins + maxs) / 2
    radius = float(np.max(maxs - mins) / 2) * 1.08
    for ax in axes:
        ax.set_xlim(center[0] - radius, center[0] + radius)
        ax.set_ylim(center[1] - radius, center[1] + radius)
        ax.set_zlim(center[2] - radius, center[2] + radius)


def structure_projection(reference: dict[int, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    arr, _ = ordered_coords(reference)
    center = arr.mean(axis=0)
    _, _, vt = np.linalg.svd(arr - center, full_matrices=False)
    return center, vt[:2].T


def project_ca(coords: dict[int, np.ndarray], center: np.ndarray, basis: np.ndarray) -> dict[int, np.ndarray]:
    return {resid: (coord - center) @ basis for resid, coord in coords.items()}


def plot_trace_2d(ax: plt.Axes, coords: dict[int, np.ndarray], color: str, title: str, subtitle: str = "") -> None:
    arr, _ = ordered_coords(coords)
    ax.plot(arr[:, 0], arr[:, 1], color=color, lw=1.2)
    ax.scatter(arr[0, 0], arr[0, 1], s=9, color=color, edgecolor="white", linewidth=0.3, zorder=3)
    ax.set_title(title, pad=1, fontsize=7.2)
    if subtitle:
        ax.text(0.5, 0.02, subtitle, transform=ax.transAxes, ha="center", fontsize=5.8, color=COLOR["gray"])
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")


def set_2d_bounds(axes: list[plt.Axes], coord_sets: list[dict[int, np.ndarray]]) -> None:
    all_coords = np.vstack([ordered_coords(coords)[0] for coords in coord_sets])
    mins = all_coords.min(axis=0)
    maxs = all_coords.max(axis=0)
    center = (mins + maxs) / 2
    radius = float(np.max(maxs - mins) / 2) * 1.08
    for ax in axes:
        ax.set_xlim(center[0] - radius, center[0] + radius)
        ax.set_ylim(center[1] - radius, center[1] + radius)


def figure1_structure_model_comparison() -> None:
    """Crystal, Chai-1, ColabFold/AF2 status, and RoseTTAFold2 comparison."""

    structure_dir = STRUCTURES
    pred_dir = structure_dir / "predictions"
    crystal = parse_pdb_ca(structure_dir / "7S7U.pdb")
    chai = parse_cif_ca(pred_dir / "chai1" / "ee88e2c3-pred.model_idx_4.cif")
    rf2 = parse_pdb_ca(pred_dir / "rosettafold2" / "63d6adf3-rf2_seed7_00_pred.pdb")
    chai_aligned, chai_rmsd, chai_n = align_to_reference(chai, crystal)
    rf2_aligned, rf2_rmsd, rf2_n = align_to_reference(rf2, crystal)
    center, basis = structure_projection(crystal)
    crystal_plot = project_ca(crystal, center, basis)
    chai_plot = project_ca(chai_aligned, center, basis)
    rf2_plot = project_ca(rf2_aligned, center, basis)

    chai_manifest = json.loads((pred_dir / "ee88e2c3-manifest.json").read_text())
    chai_model4 = next(item for item in chai_manifest["chai1"] if item["model"] == "4")
    chai_scores = chai_model4["score_summary"]
    rf2_manifest = json.loads((pred_dir / "63d6adf3-rosettafold2-manifest.json").read_text())
    rf2_plddt = rf2_manifest["rosettafold2"]["mean_plddt"]

    fig = plt.figure(figsize=(7.1, 4.55), constrained_layout=True)
    gs = fig.add_gridspec(2, 4, height_ratios=[1.0, 1.25], hspace=0.08, wspace=0.08)
    ax_crystal = fig.add_subplot(gs[0, 0])
    ax_chai = fig.add_subplot(gs[0, 1])
    ax_af2 = fig.add_subplot(gs[0, 2])
    ax_rf2 = fig.add_subplot(gs[0, 3])
    ax_overlay = fig.add_subplot(gs[1, :3])
    ax_status = fig.add_subplot(gs[1, 3])

    panel(ax_crystal, "a")
    panel(ax_chai, "b")
    panel(ax_af2, "c")
    panel(ax_rf2, "d")
    panel(ax_overlay, "e")
    panel(ax_status, "f")

    plot_trace_2d(ax_crystal, crystal_plot, "#1A202C", "Crystal 7S7U", f"{len(crystal)} C-alpha atoms")
    plot_trace_2d(ax_chai, chai_plot, COLOR["blue"], "Chai-1", f"model 4, pTM {chai_scores['ptm']['mean']:.2f}")
    plot_trace_2d(ax_rf2, rf2_plot, COLOR["orange"], "RoseTTAFold2", f"mean pLDDT {rf2_plddt:.1f}")

    ax_af2.axis("off")
    ax_af2.set_title("ColabFold/AF2", pad=1, fontsize=7.2)
    ax_af2.add_patch(
        patches.FancyBboxPatch(
            (0.10, 0.23),
            0.80,
            0.48,
            boxstyle="round,pad=0.018,rounding_size=0.02",
            facecolor="#F7FAFC",
            edgecolor="#CBD5E0",
            lw=0.8,
        )
    )
    ax_af2.text(0.50, 0.55, "No model returned", ha="center", va="center", fontsize=7.0, fontweight="bold", color=COLOR["gray"])
    ax_af2.text(
        0.50,
        0.38,
        "PyPI run hit JAX/Haiku;\nofficial image canceled\nbefore output.",
        ha="center",
        va="center",
        fontsize=5.9,
        color=COLOR["gray"],
        linespacing=1.25,
    )

    for coords, color, label, lw in [
        (crystal_plot, "#1A202C", "7S7U crystal", 1.3),
        (chai_plot, COLOR["blue"], "Chai-1", 1.0),
        (rf2_plot, COLOR["orange"], "RoseTTAFold2", 1.0),
    ]:
        arr, _ = ordered_coords(coords)
        ax_overlay.plot(arr[:, 0], arr[:, 1], color=color, lw=lw, label=label, alpha=0.88)
    ax_overlay.set_title("Overlay on 7S7U crystal", pad=1, fontsize=7.2)
    ax_overlay.set_aspect("equal", adjustable="box")
    ax_overlay.axis("off")
    ax_overlay.legend(loc="lower left", bbox_to_anchor=(0.02, 0.02), frameon=False, fontsize=6)

    set_2d_bounds([ax_crystal, ax_chai, ax_rf2, ax_overlay], [crystal_plot, chai_plot, rf2_plot])

    ax_status.axis("off")
    ax_status.set_title("Run summary", pad=1, fontsize=7.2)
    rows = [
        ("Reference", "PDB 7S7U"),
        ("Chai-1", f"RMSD {chai_rmsd:.1f} A, n={chai_n}"),
        ("Chai confidence", f"pTM {chai_scores['ptm']['mean']:.2f}; aggregate {chai_scores['aggregate_score']['mean']:.3f}"),
        ("ColabFold/AF2", "attempted; no structure"),
        ("RoseTTAFold2", f"RMSD {rf2_rmsd:.1f} A, n={rf2_n}"),
        ("RF2 confidence", f"mean pLDDT {rf2_plddt:.1f}"),
    ]
    y = 0.92
    for key, value in rows:
        ax_status.text(0.02, y, key, fontsize=6.2, fontweight="bold", va="top")
        ax_status.text(0.02, y - 0.075, value, fontsize=5.9, va="top", color=COLOR["gray"])
        ax_status.plot([0.02, 0.98], [y - 0.11, y - 0.11], color="#E2E8F0", lw=0.65)
        y -= 0.15

    save(fig, "figure1_structure_model_comparison")


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
    # Sequences ship as a curated CSV under data/metadata.
    seqs = pd.read_csv(METADATA / "sensor_sequences.csv")
    return {str(row["sensor"]): str(row["sequence"]).upper() for _, row in seqs.iterrows()}


def figure3_activity_cliffs() -> None:
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
    trunk_x = np.linspace(0.06, 0.90, len(TRUNK_ORDER))
    for i, s in enumerate(TRUNK_ORDER):
        ax_tree.scatter(trunk_x[i], y0, s=112, color=COLOR["blue"], edgecolor="white", lw=0.6, zorder=3)
        ax_tree.text(trunk_x[i], y0 + 0.10, s, ha="center", va="bottom", fontsize=8.6)
        if i:
            ax_tree.plot([trunk_x[i - 1], trunk_x[i]], [y0, y0], color="#2D3748", lw=1.5, zorder=1)
    positions = {s: (trunk_x[i], y0) for i, s in enumerate(TRUNK_ORDER)}
    # (x, node_y, parent, label_dy) — label_dy staggers long sibling labels into
    # two tiers so adjacent names never share a horizontal band.
    branch_positions = {
        "L194": (0.25, 0.40, "cc93", -0.075),
        "AK1": (0.37, 0.40, "cc93", -0.075),
        "Fent2 436L": (0.37, 0.13, "AK1", -0.075),
        "iEscSnFR": (0.49, 0.40, "V7", -0.075),
        "Tap1.0": (0.60, 0.40, "V7", -0.075),
        "iCytSnFR": (0.70, 0.36, "V9", -0.075),
        "iCytBrEtSnFR": (0.82, 0.24, "V9", -0.075),
        "iFloxSnFR": (0.94, 0.36, "V9", -0.075),
        "iLevaphenolSnFR1.0": (1.06, 0.24, "V9", -0.075),
    }
    for s, (x, y, parent, label_dy) in branch_positions.items():
        px, py = positions[parent]
        positions[s] = (x, y)
        color = COLOR["red"] if s == "Tap1.0" else COLOR["teal"]
        joint_y = 0.55 if parent in {"cc93", "V7"} else (0.46 if parent == "V9" else 0.26)
        ax_tree.plot([px, px, x, x], [py, joint_y, joint_y, y], color="#2D3748", lw=1.35, zorder=1)
        ax_tree.scatter(x, y, s=76, color=color, edgecolor="white", lw=0.5, zorder=3)
        ax_tree.text(x, y + label_dy, SHORT_SENSOR_LABELS.get(s, s), ha="center", va="top", fontsize=7.8)
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
    ax_fp.text(0.02, 0.93, "dotted: 0.3; dashed: 1.0", transform=ax_fp.transAxes, fontsize=7.2, va="top")
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


def figure4_translational_relevance() -> None:
    smiles = pd.read_csv(METADATA / "ligand_smiles.csv").set_index("ligand")["SMILES"]
    one_mutation_ec50_gain = 53.656326 / 22.383613
    one_mutation_slope_gain = (2.751028 / 22.383613) / (2.182608 / 53.656326)
    compounds = [
        {
            "ligand": "ciprofloxacin hydrochloride hydrate",
            "name": "Ciprofloxacin",
            "use": "Clinical drug monitoring",
            "context": "EC50 22.4 µM vs plasma Cmax ~9 µM",
            "required_gain": 22.383613 / 9.0,
            "lower_bound": False,
            "gap_label": "2.5x",
        },
        {
            "ligand": "l-carnitine",
            "name": "L-carnitine",
            "use": "Metabolic health",
            "context": "no plateau by 200 µM; plasma to 75 µM",
            "required_gain": 200.0 / 75.0,
            "lower_bound": True,
            "gap_label": ">2.7x",
        },
        {
            "ligand": "L-(+)-ergothioneine",
            "name": "Ergothioneine",
            "use": "Nutrition / oxidative stress",
            "context": "no plateau by 200 µM; RBC to 160 µM",
            "required_gain": 200.0 / 160.0,
            "lower_bound": True,
            "gap_label": ">1.25x",
        },
        {
            "ligand": "thiamine (hydrochloride)",
            "name": "Thiamine",
            "use": "Vitamin status",
            "context": "no plateau by 200 µM; blood to 0.195 µM",
            "required_gain": 200.0 / 0.195,
            "lower_bound": True,
            "gap_label": ">1,000x",
        },
        {
            "ligand": "DEHP",
            "name": "DEHP",
            "use": "Environmental exposure",
            "context": "EC50 34.5 µM; background serum <0.005 µM",
            "required_gain": 34.521475 / 0.005,
            "lower_bound": True,
            "gap_label": ">6,900x",
        },
    ]

    fig = plt.figure(figsize=(7.4, 4.9), constrained_layout=True)
    gs = fig.add_gridspec(5, 3, width_ratios=[1.05, 1.55, 4.4], hspace=0.08, wspace=0.06)
    ax_gain = fig.add_subplot(gs[:, 2])

    y_positions = np.arange(len(compounds))[::-1]
    rows = []
    for row_index, (compound, y_pos) in enumerate(zip(compounds, y_positions)):
        ax_mol = fig.add_subplot(gs[row_index, 0])
        ax_info = fig.add_subplot(gs[row_index, 1])
        ax_mol.axis("off")
        ax_info.axis("off")

        mol = Chem.MolFromSmiles(smiles.loc[compound["ligand"]])
        fragments = Chem.GetMolFrags(mol, asMols=True)
        mol = max(fragments, key=lambda fragment: fragment.GetNumHeavyAtoms())
        image = Draw.MolToImage(mol, size=(520, 260), fitImage=True)
        ax_mol.imshow(image)

        ax_info.text(0.0, 0.72, compound["name"], fontsize=10.1, fontweight="bold", va="center")
        ax_info.text(0.0, 0.46, compound["use"], fontsize=8.3, color=COLOR["blue"], va="center")
        ax_info.text(0.0, 0.19, compound["context"], fontsize=7.4, color=COLOR["gray"], va="center", wrap=True)

        required_gain = compound["required_gain"]
        if compound["lower_bound"]:
            ax_gain.annotate(
                "",
                xy=(required_gain, y_pos),
                xytext=(1.0, y_pos),
                arrowprops=dict(arrowstyle="-|>", color="#718096", lw=4.0, shrinkA=0, shrinkB=0),
            )
        else:
            ax_gain.plot(
                [1.0, required_gain],
                [y_pos, y_pos],
                color=COLOR["red"],
                lw=4.5,
                solid_capstyle="butt",
            )
            ax_gain.scatter([required_gain], [y_pos], s=42, color=COLOR["red"], edgecolor="white", lw=0.5, zorder=4)

        # Push the label clear of the one-mutation band so text never sits on the
        # dotted band verticals or the red measured-requirement marker.
        band_right = 53.656326 / 22.383613  # upper edge of the blue band (~2.4x)
        label_x = max(required_gain, band_right) * 1.22
        ax_gain.text(
            label_x,
            y_pos,
            compound["gap_label"],
            fontsize=8.3,
            ha="left",
            va="center",
            color="#2D3748",
        )
        rows.append(
            {
                "compound": compound["name"],
                "use_case": compound["use"],
                "relevant_context": compound["context"],
                "minimum_required_gain": required_gain,
                "is_lower_bound": compound["lower_bound"],
                "gap_label": compound["gap_label"],
            }
        )

    for boundary in np.arange(0.5, len(compounds) - 0.5, 1.0):
        ax_gain.axhline(boundary, color="#E2E8F0", lw=0.7)

    ax_gain.axvspan(one_mutation_ec50_gain, one_mutation_slope_gain, color=COLOR["sky"], alpha=0.24, zorder=0)
    ax_gain.axvline(one_mutation_ec50_gain, color=COLOR["blue"], lw=0.8, ls=":")
    ax_gain.axvline(one_mutation_slope_gain, color=COLOR["blue"], lw=0.8, ls=":")
    benchmark_midpoint = math.sqrt(one_mutation_ec50_gain * one_mutation_slope_gain)
    ax_gain.text(
        25.0,
        0.86,
        "Blue band: 1-mutation gain observed\nin this lineage (V7→iEscSnFR, T360S)\n2.4x EC50 gain; 3.0x sensitivity gain",
        transform=ax_gain.get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=7.0,
        color=COLOR["blue"],
    )

    ax_gain.set_xscale("log")
    ax_gain.set_xlim(0.8, 20000)
    ax_gain.set_ylim(-0.55, len(compounds) - 0.45)
    ax_gain.set_yticks([])
    ax_gain.set_xticks([1, 3, 10, 100, 1000, 10000])
    ax_gain.set_xticklabels(["1", "3", "10", "100", "1,000", "10,000"])
    ax_gain.set_xlabel("Sensitivity improvement required (fold, log scale)")
    ax_gain.grid(axis="x", which="major", color="#EDF2F7", lw=0.7)
    ax_gain.set_axisbelow(True)
    ax_gain.spines["top"].set_visible(False)
    ax_gain.spines["right"].set_visible(False)
    ax_gain.spines["left"].set_visible(False)
    ax_gain.legend(
        handles=[
            Line2D([0], [0], color=COLOR["red"], marker="o", lw=4, markersize=4, label="measured requirement"),
            Line2D([0], [0], color="#718096", marker=">", lw=4, markersize=5, label="minimum requirement"),
            patches.Patch(facecolor=COLOR["sky"], alpha=0.24, label="one-mutation gain observed"),
        ],
        frameon=False,
        ncol=3,
        loc="upper center",
        bbox_to_anchor=(0.57, 1.085),
        fontsize=7.3,
        handlelength=2.0,
        columnspacing=1.0,
    )

    fig.suptitle("Small sequence changes can close near-term sensitivity gaps", fontsize=11.0)
    rows.append(
        {
            "compound": "one-mutation benchmark",
            "use_case": "V7 to iEscSnFR, T360S",
            "relevant_context": "ciprofloxacin",
            "minimum_required_gain": one_mutation_ec50_gain,
            "is_lower_bound": False,
            "gap_label": f"{one_mutation_ec50_gain:.2f}x EC50; {one_mutation_slope_gain:.2f}x delta-slope",
        }
    )
    pd.DataFrame(rows).to_csv(ANALYSIS / "figure4_translational_context.csv", index=False)
    save(fig, "figure5_translational_relevance")


def main() -> None:
    ANALYSIS.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    style()
    figure2_screen_scope()
    figure3_activity_cliffs()
    figure4_dose_response_leads()
    figure4_translational_relevance()


if __name__ == "__main__":
    main()
