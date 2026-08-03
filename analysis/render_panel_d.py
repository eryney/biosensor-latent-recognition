"""Figure 1 panel (d) - ESMFold prediction of cc93 overlaid on the 7S7U crystal.

The prediction is superposed on the PBP domain only (resi 3-84 + 329-522), which
is the comparison the panel makes: the PBP lobes track the crystal closely
(Ca RMSD 1.7 A) while the cpGFP module is displaced by ~21 A. Crystal is wheat,
prediction is blue.

Run in the `pymol-render` environment; writes ../figures/panel_d_7S7U.png.
"""
import pymol
pymol.finish_launching(['pymol', '-qc'])
from pymol import cmd

CRY = "../data/structures/7S7U.pdb"
ESM = "../data/structures/esmfold_cc93_fulllength.pdb"
OUT = "../figures/panels/panel_d_7S7U.png"

cmd.load(CRY, "cry")
cmd.load(ESM, "esm")
cmd.remove("solvent")
cmd.remove("cry and not chain A")
cmd.remove("not (alt ''+A)"); cmd.alter("all", "alt=''")
cmd.hide("everything"); cmd.bg_color("white")
cmd.show("cartoon", "cry or esm")

# Superpose the prediction onto the crystal over the PBP domain only.
pbp_sel = "resi 3-84+329-522 and name CA"
rms = cmd.align("esm and " + pbp_sel, "cry and " + pbp_sel, cycles=0)[0]
print("PBP-domain align RMSD:", round(rms, 3))

# Crystal wheat, prediction a single uniform blue: the divergence should read
# from the geometry after PBP superposition rather than from a second hue.
cmd.color("wheat", "cry")
cmd.color("skyblue", "esm")

cmd.set("cartoon_transparency", 0.0, "cry")
cmd.set("cartoon_transparency", 0.15, "esm")
cmd.set("ray_shadows", 0); cmd.set("ray_opaque_background", 0); cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1)

cmd.orient("cry")
cmd.zoom("cry or esm", buffer=4)
cmd.ray(2200, 1700); cmd.png(OUT, dpi=300)
print("wrote", OUT)
