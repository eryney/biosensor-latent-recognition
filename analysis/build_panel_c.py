"""Figure 1 panel (c) - ciprofloxacin docked in the 7S7U binding pocket.

Writes both a reference render and an editable PyMOL session. The aromatic-cage
residues are exposed as named selections (cage_Y65 ... cage_Y460) with labels
pre-placed at each CB so they can be repositioned without rebuilding the scene.

The pose is the rank-1 result of a blind DiffDock-L dock into 7S7U (the closed,
ligand-bound conformation) and is an illustrative predicted pose with unrelaxed
geometry, not an experimental structure.

Run in the `pymol-render` environment; writes ../figures/panels/panel_c_7S7U_*.
"""
import pymol
pymol.finish_launching(['pymol', '-qc'])
from pymol import cmd

REC = "../data/structures/7S7U_receptor.pdb"
LIG = "../data/structures/cipro_docked_7S7U_rank1.sdf"
PSE = "../figures/panels/panel_c_7S7U_forlabels.pse"
PNG = "../figures/panels/panel_c_7S7U_reference.png"

CAGE = {65:"Y65", 357:"Y357", 360:"T360", 391:"F391", 436:"W436", 460:"Y460"}

cmd.load(REC, "rec")
cmd.load(LIG, "cipro")
cmd.remove("solvent")
cmd.hide("everything"); cmd.bg_color("white")

# receptor as thin grey cartoon context
cmd.show("cartoon", "rec")
cmd.set("cartoon_transparency", 0.55, "rec")
cmd.color("grey80", "rec")

# chromophore for spatial reference
cmd.show("sticks", "rec and resn CRO")
cmd.color("palegreen", "rec and resn CRO and elem C")

# aromatic cage side chains as sticks, named selections so labels snap
for resi, lab in CAGE.items():
    sel = f"rec and resi {resi} and not (name N+C+O)"
    cmd.show("sticks", sel)
    cmd.color("slate", sel + " and elem C")
    cmd.select(f"cage_{lab}", f"rec and resi {resi}")
# de-select
cmd.deselect()

# ciprofloxacin — yellow sticks, prominent
cmd.show("sticks", "cipro")
cmd.set("stick_radius", 0.22, "cipro")
cmd.color("yellow", "cipro and elem C")
cmd.set("stick_radius", 0.20, "rec and (" + " or ".join(f"resi {r}" for r in CAGE) + ")")

cmd.set("ray_shadows", 0); cmd.set("ray_opaque_background", 0); cmd.set("antialias", 2)
cmd.set("cartoon_side_chain_helper", 1)

# zoom to pocket = cage residues + ligand
pocket_sel = "cipro or (rec and resi " + "+".join(map(str, CAGE)) + ")"
cmd.orient(pocket_sel)
cmd.zoom(pocket_sel, buffer=3.5)

# Labels at the CB of each cage residue. Y357 and W436 sit close together in
# this view, so their labels are offset along opposite directions to keep them
# legible; the offsets are in Angstroms in camera space.
LABEL_OFFSET = {
    65: [-1.6, 1.5, 0.0],
    357: [1.7, 2.0, 0.0],
    360: [1.5, 0.9, 0.0],
    391: [-2.2, -1.6, 0.0],
    436: [2.0, -1.7, 0.0],
    460: [0.6, 2.1, 0.0],
}

cmd.set("label_size", 20)
cmd.set("label_color", "black")
cmd.set("label_font_id", 7)
# An outline in the background colour keeps labels readable where they cross
# the pale cartoon.
cmd.set("label_outline_color", "white")
cmd.set("label_bg_color", "white")
cmd.set("label_bg_transparency", 0.25)
for resi, lab in CAGE.items():
    sel = f"rec and resi {resi} and name CB"
    cmd.label(sel, f'"{lab}"')
    cmd.set("label_position", LABEL_OFFSET[resi], sel)
if cmd.count_atoms("cipro and name N1"):
    cmd.label("cipro and name N1", '"ciprofloxacin"')
    cmd.set("label_position", [0.0, -2.4, 0.0], "cipro and name N1")

cmd.save(PSE)
cmd.ray(2000, 1600); cmd.png(PNG, dpi=300)
print("wrote", PSE, "and", PNG)
