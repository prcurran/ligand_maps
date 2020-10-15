
try:
    import tkinter as tk      
except ImportError:
    import Tkinter as tk
from os.path import join
import tempfile

import zipfile
import math
from pymol import cmd, finish_launching, plugins
from pymol.cgo import *

finish_launching()

cmd.load("ligands.mol2", "ligands")
donor_projected_point_0 = [COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 0.8] + [SPHERE, float(21.0), float(0.0), float(14.5), float(1)]

cmd.load_cgo(donor_projected_point_0, "donor_projected_point_0_obj", 1)
cmd.pseudoatom(object="donor_projected_point_0_score", pos=(21.0, 0.0, 14.5), color=(1, 1, 1), label=27.1)

cmd.pseudoatom(object="donor_projected_point_0_proj_score", pos=(21.5, -0.25, 12.25), color=(1, 1, 1), label=14.7)

donor_projected_projection_0 = [COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 0.8] + [SPHERE, float(21.5), float(-0.25), float(12.25), float(1)]

cmd.load_cgo(donor_projected_projection_0, "donor_projected_projection_0_obj", 1)
cmd.pseudoatom(object="donor_projected_line_0pa1", pos=(21.0, 0.0, 14.5), color=(1, 1, 1))

cmd.pseudoatom(object="donor_projected_line_0pa2", pos=(21.5, -0.25, 12.25), color=(1, 1, 1))

cmd.distance(name="donor_projected_line_0", selection1="donor_projected_line_0pa1", selection2="donor_projected_line_0pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_projected_line_0")
cmd.set("dash_width", 4.0)
cmd.delete("donor_projected_line_0pa1")
cmd.delete("donor_projected_line_0pa2")
cmd.group("donor_projected_0", members="donor_projected_point_0_obj")
cmd.group("donor_projected_0", members="donor_projected_point_0_score")
cmd.group("donor_projected_0", members="donor_projected_projection_0")
cmd.group("donor_projected_0", members="donor_projected_line_0")
cmd.group("donor_projected_0", members="donor_projected_point_0_proj_score")
acceptor_projected_point_1 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 0.8] + [SPHERE, float(22.0), float(-1.5), float(19.0), float(1)]

cmd.load_cgo(acceptor_projected_point_1, "acceptor_projected_point_1_obj", 1)
cmd.pseudoatom(object="acceptor_projected_point_1_score", pos=(22.0, -1.5, 19.0), color=(1, 1, 1), label=21.9)

cmd.pseudoatom(object="acceptor_projected_point_1_proj_score", pos=(24.0, -2.0, 21.0), color=(1, 1, 1), label=18.3)

acceptor_projected_projection_1 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 0.8] + [SPHERE, float(24.0), float(-2.0), float(21.0), float(1)]

cmd.load_cgo(acceptor_projected_projection_1, "acceptor_projected_projection_1_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_1pa1", pos=(22.0, -1.5, 19.0), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_1pa2", pos=(24.0, -2.0, 21.0), color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_1", selection1="acceptor_projected_line_1pa1", selection2="acceptor_projected_line_1pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_1")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_1pa1")
cmd.delete("acceptor_projected_line_1pa2")
cmd.group("acceptor_projected_1", members="acceptor_projected_point_1_obj")
cmd.group("acceptor_projected_1", members="acceptor_projected_point_1_score")
cmd.group("acceptor_projected_1", members="acceptor_projected_projection_1")
cmd.group("acceptor_projected_1", members="acceptor_projected_line_1")
cmd.group("acceptor_projected_1", members="acceptor_projected_point_1_proj_score")
acceptor_projected_point_2 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 0.8] + [SPHERE, float(22.0), float(-1.5), float(19.0), float(1)]

cmd.load_cgo(acceptor_projected_point_2, "acceptor_projected_point_2_obj", 1)
cmd.pseudoatom(object="acceptor_projected_point_2_score", pos=(22.0, -1.5, 19.0), color=(1, 1, 1), label=21.9)

cmd.pseudoatom(object="acceptor_projected_point_2_proj_score", pos=(20.25, -3.25, 18.25), color=(1, 1, 1), label=17.7)

acceptor_projected_projection_2 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 0.8] + [SPHERE, float(20.25), float(-3.25), float(18.25), float(1)]

cmd.load_cgo(acceptor_projected_projection_2, "acceptor_projected_projection_2_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_2pa1", pos=(22.0, -1.5, 19.0), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_2pa2", pos=(20.25, -3.25, 18.25), color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_2", selection1="acceptor_projected_line_2pa1", selection2="acceptor_projected_line_2pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_2")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_2pa1")
cmd.delete("acceptor_projected_line_2pa2")
cmd.group("acceptor_projected_2", members="acceptor_projected_point_2_obj")
cmd.group("acceptor_projected_2", members="acceptor_projected_point_2_score")
cmd.group("acceptor_projected_2", members="acceptor_projected_projection_2")
cmd.group("acceptor_projected_2", members="acceptor_projected_line_2")
cmd.group("acceptor_projected_2", members="acceptor_projected_point_2_proj_score")
acceptor_projected_point_3 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 0.8] + [SPHERE, float(19.5), float(-4.5), float(16.5), float(1)]

cmd.load_cgo(acceptor_projected_point_3, "acceptor_projected_point_3_obj", 1)
cmd.pseudoatom(object="acceptor_projected_point_3_score", pos=(19.5, -4.5, 16.5), color=(1, 1, 1), label=18.4)

cmd.pseudoatom(object="acceptor_projected_point_3_proj_score", pos=(21.25, -3.25, 18.25), color=(1, 1, 1), label=16.3)

acceptor_projected_projection_3 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 0.8] + [SPHERE, float(21.25), float(-3.25), float(18.25), float(1)]

cmd.load_cgo(acceptor_projected_projection_3, "acceptor_projected_projection_3_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_3pa1", pos=(19.5, -4.5, 16.5), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_3pa2", pos=(21.25, -3.25, 18.25), color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_3", selection1="acceptor_projected_line_3pa1", selection2="acceptor_projected_line_3pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_3")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_3pa1")
cmd.delete("acceptor_projected_line_3pa2")
cmd.group("acceptor_projected_3", members="acceptor_projected_point_3_obj")
cmd.group("acceptor_projected_3", members="acceptor_projected_point_3_score")
cmd.group("acceptor_projected_3", members="acceptor_projected_projection_3")
cmd.group("acceptor_projected_3", members="acceptor_projected_line_3")
cmd.group("acceptor_projected_3", members="acceptor_projected_point_3_proj_score")
acceptor_projected_point_4 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 0.8] + [SPHERE, float(19.5), float(-4.5), float(16.5), float(1)]

cmd.load_cgo(acceptor_projected_point_4, "acceptor_projected_point_4_obj", 1)
cmd.pseudoatom(object="acceptor_projected_point_4_score", pos=(19.5, -4.5, 16.5), color=(1, 1, 1), label=18.4)

cmd.pseudoatom(object="acceptor_projected_point_4_proj_score", pos=(17.5, -5.75, 18.0), color=(1, 1, 1), label=14.7)

acceptor_projected_projection_4 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 0.8] + [SPHERE, float(17.5), float(-5.75), float(18.0), float(1)]

cmd.load_cgo(acceptor_projected_projection_4, "acceptor_projected_projection_4_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_4pa1", pos=(19.5, -4.5, 16.5), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_4pa2", pos=(17.5, -5.75, 18.0), color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_4", selection1="acceptor_projected_line_4pa1", selection2="acceptor_projected_line_4pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_4")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_4pa1")
cmd.delete("acceptor_projected_line_4pa2")
cmd.group("acceptor_projected_4", members="acceptor_projected_point_4_obj")
cmd.group("acceptor_projected_4", members="acceptor_projected_point_4_score")
cmd.group("acceptor_projected_4", members="acceptor_projected_projection_4")
cmd.group("acceptor_projected_4", members="acceptor_projected_line_4")
cmd.group("acceptor_projected_4", members="acceptor_projected_point_4_proj_score")
acceptor_projected_point_5 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 0.8] + [SPHERE, float(21.0), float(0.0), float(14.5), float(1)]

cmd.load_cgo(acceptor_projected_point_5, "acceptor_projected_point_5_obj", 1)
cmd.pseudoatom(object="acceptor_projected_point_5_score", pos=(21.0, 0.0, 14.5), color=(1, 1, 1), label=17.8)

cmd.pseudoatom(object="acceptor_projected_point_5_proj_score", pos=(22.5, 2.25, 15.75), color=(1, 1, 1), label=10.8)

acceptor_projected_projection_5 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 0.8] + [SPHERE, float(22.5), float(2.25), float(15.75), float(1)]

cmd.load_cgo(acceptor_projected_projection_5, "acceptor_projected_projection_5_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_5pa1", pos=(21.0, 0.0, 14.5), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_5pa2", pos=(22.5, 2.25, 15.75), color=(1, 1, 1))

cmd.distance(name="acceptor_projected_line_5", selection1="acceptor_projected_line_5pa1", selection2="acceptor_projected_line_5pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (1.0, 0.6, 0.6), selection="acceptor_projected_line_5")
cmd.set("dash_width", 4.0)
cmd.delete("acceptor_projected_line_5pa1")
cmd.delete("acceptor_projected_line_5pa2")
cmd.group("acceptor_projected_5", members="acceptor_projected_point_5_obj")
cmd.group("acceptor_projected_5", members="acceptor_projected_point_5_score")
cmd.group("acceptor_projected_5", members="acceptor_projected_projection_5")
cmd.group("acceptor_projected_5", members="acceptor_projected_line_5")
cmd.group("acceptor_projected_5", members="acceptor_projected_point_5_proj_score")
cmd.group("donor_projected_pts", members="donor_projected_0")
cmd.group("acceptor_projected_pts", members="acceptor_projected_1")
cmd.group("acceptor_projected_pts", members="acceptor_projected_2")
cmd.group("acceptor_projected_pts", members="acceptor_projected_3")
cmd.group("acceptor_projected_pts", members="acceptor_projected_4")
cmd.group("acceptor_projected_pts", members="acceptor_projected_5")
cmd.group("ligand_pharmacophore", members="donor_projected_pts")
cmd.group("ligand_pharmacophore", members="guanine_pts")
cmd.group("ligand_pharmacophore", members="uracil_pts")
cmd.group("ligand_pharmacophore", members="acceptor_pts")
cmd.group("ligand_pharmacophore", members="deoxyribose_pts")
cmd.group("ligand_pharmacophore", members="donor_ch_projected_pts")
cmd.group("ligand_pharmacophore", members="acceptor_projected_pts")
cmd.group("ligand_pharmacophore", members="cytosine_pts")
cmd.group("ligand_pharmacophore", members="ribose_pts")
cmd.group("ligand_pharmacophore", members="purine_pts")
cmd.group("ligand_pharmacophore", members="adenine_pts")
cmd.group("ligand_pharmacophore", members="ring_planar_projected_pts")
cmd.group("ligand_pharmacophore", members="pyrimidine_pts")
cmd.group("ligand_pharmacophore", members="thymine_pts")
cmd.group("ligand_pharmacophore", members="ring_pts")
cmd.group("ligand_pharmacophore", members="ring_non_planar_pts")
cmd.group("ligand_pharmacophore", members="hydrophobe_pts")
cmd.group("ligand_pharmacophore", members="heavy_atom_pts")
cmd.group("ligand_pharmacophore", members="bromine_pts")
cmd.group("ligand_pharmacophore", members="exit_vector_pts")
cmd.group("ligand_pharmacophore", members="halogen_pts")
cmd.group("ligand_pharmacophore", members="water_pts")
cmd.group("ligand_pharmacophore", members="fluorine_pts")
cmd.group("ligand_pharmacophore", members="metal_pts")
cmd.group("ligand_pharmacophore", members="chlorine_pts")
cmd.group("ligand_pharmacophore", members="GLN_pts")
cmd.group("ligand_pharmacophore", members="PRO_pts")
cmd.group("ligand_pharmacophore", members="ALA_pts")
cmd.group("ligand_pharmacophore", members="TYR_pts")
cmd.group("ligand_pharmacophore", members="SER_pts")
cmd.group("ligand_pharmacophore", members="GLY_pts")
cmd.group("ligand_pharmacophore", members="ASP_pts")
cmd.group("ligand_pharmacophore", members="CYS_pts")
cmd.group("ligand_pharmacophore", members="PHE_pts")
cmd.group("ligand_pharmacophore", members="LEU_pts")
cmd.group("ligand_pharmacophore", members="LYS_pts")
cmd.group("ligand_pharmacophore", members="ILE_pts")
cmd.group("ligand_pharmacophore", members="THR_pts")
cmd.group("ligand_pharmacophore", members="TRP_pts")
cmd.group("ligand_pharmacophore", members="ASN_pts")
cmd.group("ligand_pharmacophore", members="VAL_pts")
cmd.group("ligand_pharmacophore", members="HIS_pts")
cmd.group("ligand_pharmacophore", members="MET_pts")
cmd.group("ligand_pharmacophore", members="ARG_pts")
cmd.group("ligand_pharmacophore", members="GLU_pts")
cmd.set_color("ring_color", (0.33, 1.0, 0.0))
cmd.load("ring.grd", "ring_grid")
cmd.isosurface(name="surface_ring", map="ring_grid", level="1")

cmd.color("ring_color", "surface_ring")
cmd.set_color("acceptor_projected_color", (1.0, 0.0, 0.0))
cmd.load("acceptor_projected.grd", "acceptor_projected_grid")
cmd.isosurface(name="surface_acceptor_projected", map="acceptor_projected_grid", level="1")

cmd.color("acceptor_projected_color", "surface_acceptor_projected")
cmd.set_color("donor_projected_color", (0.0, 0.0, 1.0))
cmd.load("donor_projected.grd", "donor_projected_grid")
cmd.isosurface(name="surface_donor_projected", map="donor_projected_grid", level="1")

cmd.color("donor_projected_color", "surface_donor_projected")
cmd.group("feature_grids", members="ring")
cmd.group("feature_grids", members="acceptor_projected")
cmd.group("feature_grids", members="donor_projected")
cmd.group("feature_grids", members="surface_ring")
cmd.group("feature_grids", members="surface_acceptor_projected")
cmd.group("feature_grids", members="surface_donor_projected")


class IsoLevel(tk.Variable):
    def __init__(self, master, name, level):
        tk.Variable.__init__(self, master, value=level)
        self.name = name
        self.trace('w', self.callback)

    def callback(self, *args):
        cmd.isolevel(self.name, self.get())

    def increment(self, event=None, delta=0.1):
        self.set(round(float(self.get()) + delta, 2))

    def decrement(self, event=None):
        self.increment(None, -0.1)


surface_list = {None: {'feature_grids': ['surface_ring', 'surface_acceptor_projected', 'surface_donor_projected']}}
surface_max_list = {None: {'feature_grids': 27.1}}

top = tk.Toplevel(plugins.get_tk_root())

master = tk.Frame(top, padx=10, pady=10)
master.pack(fill="both", expand=1)

for child in list(master.children.values()):
    child.destroy()


row_counter = 0
for identifier, component_dic in surface_list.items():
    # add calculation identifier
    tk.Label(master, text=identifier).grid(row=row_counter, column=0, sticky="w")
    row_counter += 1
    
    for component_id, surfaces in component_dic.items():
        # add collection label, e.g. superstar or hotspot etc.
        tk.Label(master, text=component_id).grid(row=row_counter, column=1, sticky='w')
        row_counter += 1
        
        for i, surface in enumerate(surfaces):
            # add grid type label
            probe = surface.split("_")[-2]
            tk.Label(master, text=probe).grid(row=row_counter, column=2, sticky="w")
            
            # slider code 
            v = IsoLevel(master, surface, 5)
            e = tk.Scale(master, orient=tk.HORIZONTAL, from_=0, to=surface_max_list[identifier][component_id],
                         resolution=0.1, showvalue=0, variable=v)
            e.grid(row=row_counter, column=3, sticky="ew")

            e = tk.Entry(master, textvariable=v, width=4)
            e.grid(row=row_counter, column=4, sticky="e")
            master.columnconfigure(3, weight=1)
            row_counter += 1



