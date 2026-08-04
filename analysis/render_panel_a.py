"""Figure 1 panel (a) - domain architecture of the OpuBC-cpGFP scaffold.

Rendered from 7S7U chain A (iNicSnFR3a / cc93, closed conformation). PBP is
magenta, the cpGFP barrel green, the connecting linker cyan; the 31 positions
that vary across the scaffold family are orange CA spheres and the four
aromatic-cage residues are yellow sticks.

The reporter is a circularly permuted GFP, so the barrel is discontinuous in
crystal numbering: it spans resi 85-319, with the PBP contributing 3-84 and
329-523. Module boundaries were assigned by aligning chain A against avGFP.

Run in the `pymol-render` environment; writes ../figures/panels/panel_a_7S7U.png.
"""
import pymol, numpy as np
pymol.finish_launching(['pymol', '-qc'])
from pymol import cmd

PDB = "../data/structures/7S7U.pdb"
OUT = "../figures/panels/panel_a_7S7U.png"

# Module boundaries (7S7U chain A, verified by avGFP alignment)
CPGFP = "85-319"            # circularly permuted GFP barrel
PBP_N = "3-84"             # PBP N-portion before the GFP insertion
PBP_C = "329-523"          # PBP C-portion after the connecting linker
LINK_C = "320-328"         # modeled connecting linker (FPPPSSTDP)
# cp-internal linker 169-175 is unmodeled; N-junction 84/85 is abrupt

MUT = [10,11,12,15,22,43,44,63,66,68,76,84,100,102,196,324,325,355,356,357,358,
       360,389,391,395,418,421,436,455,490,500]
CAGE = {65:"Y65",357:"Y357",391:"F391",436:"W436"}

cmd.load(PDB, "m")
cmd.remove("solvent")
cmd.remove("not (alt ''+A)"); cmd.alter("all", "alt=''")
cmd.remove("m and not chain A")
cmd.hide("everything"); cmd.bg_color("white")
cmd.show("cartoon", "m and polymer.protein")

# --- corrected coloring ---
cmd.color("magenta", "m and polymer.protein")          # PBP default
cmd.color("green",   f"m and resi {CPGFP}")            # cpGFP barrel (85-319)
cmd.color("cyan",    f"m and resi {LINK_C}")           # connecting linker
# NO cartoon_transparency (was 0.10) -> nothing bleeds through the barrel

# chromophore
cmd.show("sticks", "m and resn CRO")
cmd.color("limon", "m and resn CRO and elem C")

# mutation sites
sel_mut = "m and name CA and resi " + "+".join(map(str, MUT))
cmd.show("spheres", sel_mut); cmd.set("sphere_scale", 0.55, sel_mut)
cmd.color("orange", sel_mut)

# aromatic cage
for resi in CAGE:
    s = f"m and resi {resi}"
    cmd.show("sticks", s + " and not name N+C+O")
    cmd.color("yellow", s + " and elem C"); cmd.set("stick_radius", 0.30, s)

cmd.set("ray_shadows", 0); cmd.set("ray_opaque_background", 0); cmd.set("antialias", 2)
cmd.set("cartoon_fancy_helices", 1); cmd.set("sphere_quality", 3)

# --- deterministic camera: PBP left, cpGFP right ---
def centroid(sel):
    mdl = cmd.get_model(sel)
    return np.mean([a.coord for a in mdl.atom], axis=0)
pbp = (centroid(f"m and resi {PBP_N} and name CA") +
       centroid(f"m and resi {PBP_C} and name CA")) / 2
gfp = centroid(f"m and resi {CPGFP} and name CA")
xcam = gfp - pbp; xcam = xcam / np.linalg.norm(xcam)   # camera-right = PBP->GFP => PBP left, GFP right
worldz = np.array([0, 0, 1.0])
ycam = worldz - xcam * np.dot(worldz, xcam)
if np.linalg.norm(ycam) < 1e-3:
    ycam = np.array([0, 1.0, 0]) - xcam * xcam[1]
ycam = ycam / np.linalg.norm(ycam)
zcam = np.cross(xcam, ycam); zcam = zcam / np.linalg.norm(zcam)
R = np.vstack([xcam, ycam, zcam])
view = list(cmd.get_view())
view[0:9] = [float(x) for x in R.flatten()]
cmd.set_view(view)
cmd.zoom("m and polymer.protein", buffer=5)
# roll about the PBP->GFP axis so the barrel is viewed side-on (staves visible),
# not end-on. Small turn about camera-x (the inter-domain axis).
cmd.turn("x", -35)
cmd.zoom("m and polymer.protein", buffer=5)
# Ray tracing is pinned to one thread. PyMOL distributes ray tracing across
# threads by default, which makes the antialiased output differ by a few
# pixels between runs; single-threaded rendering is reproducible.
cmd.set("max_threads", 1)
cmd.ray(2200, 1700); cmd.png(OUT, dpi=300)
print("wrote", OUT)
