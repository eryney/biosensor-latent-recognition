# Figure 1 panels

Figure 1 is a four-panel composite built from structural renders rather than from
the screening data. Panels (a), (c) and (d) are rendered in PyMOL; panel (b) is a
schematic. `analysis/compose_figure1.py` assembles them into the final figure.

## Structure used

All structural panels are built on **7S7U**, the crystal structure of iNicSnFR3a
(cc93) in the closed, nicotine-bound conformation.

| PDB  | State  | Y65-F391 Cα distance | Ligand      |
|------|--------|----------------------|-------------|
| 7S7V | open   | 17.7 Å               | —           |
| 7S7U | closed | 13.4 Å               | nicotine    |
| 7S7T | closed | 13.5 Å               | varenicline |

The pocket is formed only in the closed state, so 7S7U is the appropriate target
for the docking panel. 7S7U and 7S7T are closed to within 0.1 Å of each other by this
measure; 7S7U was selected as the nicotine-bound structure of the sensor from
which this family descends.

## Panel (a) — scaffold architecture

`render_panel_a.py` → `panels/panel_a_7S7U.png`

7S7U chain A as cartoon: PBP magenta, cpGFP barrel green, connecting linker cyan.
The chromophore (CRO 240) is shown in limon, the four aromatic-cage residues (Y65,
Y357, F391, W436) as yellow sticks, and the 31 positions that vary across the
17-member scaffold family as orange Cα spheres.

The reporter is a circularly permuted GFP, so the barrel is discontinuous in
crystal numbering. Alignment of chain A against avGFP places the avGFP N-half
(1-145) at crystal 176-319 and the C-half (145-238) at crystal 85-168; the barrel
therefore spans residues 85-319, and the PBP contributes 3-84 together with
329-523. The internal permutation linker (169-175) is unmodelled in 7S7U, and the
connecting linker is 320-328 (FPPPSSTDP).

## Panel (b) — sensing mechanism

`render_panel_b.py` → `panels/panel_b_mechanism.png`

Schematic of the conformational change that couples ligand binding to
fluorescence output. In the unbound state the PBP lobes are open and a linker
glutamate rests against the cpGFP chromophore, quenching it; ligand binding
closes the lobes, withdrawing the glutamate and relieving the quenching. Colours
match panel (a): PBP magenta, cpGFP green, connecting linker cyan.

This panel is a diagram, not a rendering of a structure, and no part of it is
derived from the screening data.

## Panel (c) — ciprofloxacin in the binding pocket

`build_panel_c.py` → `panels/panel_c_7S7U_reference.png`, `panel_c_7S7U_forlabels.pse`

The ligand pose is the rank-1 result of a blind DiffDock-L dock (v1.1.3) into
7S7U, with contacts to W436 (1.65 Å), F391 (2.0 Å) and Y65 (3.8 Å). All ten ranked
poses are provided under `data/structures/diffdock_7S7U_ranked_poses/`.

This is a predicted pose with unrelaxed geometry. It illustrates the pocket
environment and is not experimental evidence of the binding mode; no experimental
structure of a ciprofloxacin complex was determined in this work.

The script writes both a reference render and an editable PyMOL session in which
the aromatic-cage residues are exposed as named selections (`cage_Y65` …
`cage_Y460`). Labels are placed at each Cβ with per-residue offsets, since Y357
and W436 project close together in this view and their labels would otherwise
overlap.

## Panel (d) — ESMFold prediction versus crystal

`render_panel_d.py` → `panels/panel_d_7S7U.png`

Full-length ESMFold prediction of cc93 (mean pLDDT ≈ 70) superposed on 7S7U over
the PBP domain only (residues 3-84 + 329-522). The crystal is wheat, the
prediction blue.

The PBP-domain Cα RMSD is 1.74 Å, obtained independently with Biopython's
`Superimposer` and with PyMOL `align`. After that superposition the cpGFP module
is displaced by approximately 21 Å RMSD. The panel makes this contrast explicit:
the recognition domain is predicted accurately while the placement of the reporter
module relative to it is not.

## Inputs

| File | Contents |
|------|----------|
| `data/structures/7S7U.pdb` | crystal structure |
| `data/structures/7S7U_receptor.pdb` | chain A, waters removed; docking target |
| `data/structures/cipro_docked_7S7U_rank1.sdf` | DiffDock-L rank-1 pose |
| `data/structures/diffdock_7S7U_ranked_poses/` | all ten ranked poses |
| `data/structures/esmfold_cc93_fulllength.pdb` | ESMFold prediction of cc93 |
