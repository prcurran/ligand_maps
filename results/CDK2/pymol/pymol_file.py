
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
ring_point_0 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 0.8] + [SPHERE, float(27.5), float(5.0), float(64.5), float(2.0)]

cmd.load_cgo(ring_point_0, "ring_point_0_obj", 1)
cmd.pseudoatom(object="ring_point_0_score", pos=(27.5, 5.0, 64.5), color=(1, 1, 1), label=191.4)

cmd.group("ring_0", members="ring_point_0_obj")
cmd.group("ring_0", members="ring_point_0_score")
acceptor_projected_point_1 = [COLOR, 1.0, 0.0, 0.0] +  [ALPHA, 0.8] + [SPHERE, float(28.5), float(4.5), float(65.5), float(1)]

cmd.load_cgo(acceptor_projected_point_1, "acceptor_projected_point_1_obj", 1)
cmd.pseudoatom(object="acceptor_projected_point_1_score", pos=(28.5, 4.5, 65.5), color=(1, 1, 1), label=190.6)

cmd.pseudoatom(object="acceptor_projected_point_1_proj_score", pos=(30.0, 2.75, 67.0), color=(1, 1, 1), label=135.7)

acceptor_projected_projection_1 = [COLOR, 1.0, 0.6, 0.6] +  [ALPHA, 0.8] + [SPHERE, float(30.0), float(2.75), float(67.0), float(1)]

cmd.load_cgo(acceptor_projected_projection_1, "acceptor_projected_projection_1_obj", 1)
cmd.pseudoatom(object="acceptor_projected_line_1pa1", pos=(28.5, 4.5, 65.5), color=(1, 1, 1))

cmd.pseudoatom(object="acceptor_projected_line_1pa2", pos=(30.0, 2.75, 67.0), color=(1, 1, 1))

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
donor_projected_point_2 = [COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 0.8] + [SPHERE, float(30.0), float(6.5), float(65.5), float(1)]

cmd.load_cgo(donor_projected_point_2, "donor_projected_point_2_obj", 1)
cmd.pseudoatom(object="donor_projected_point_2_score", pos=(30.0, 6.5, 65.5), color=(1, 1, 1), label=137.3)

cmd.pseudoatom(object="donor_projected_point_2_proj_score", pos=(30.75, 4.5, 67.25), color=(1, 1, 1), label=108.6)

donor_projected_projection_2 = [COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 0.8] + [SPHERE, float(30.75), float(4.5), float(67.25), float(1)]

cmd.load_cgo(donor_projected_projection_2, "donor_projected_projection_2_obj", 1)
cmd.pseudoatom(object="donor_projected_line_2pa1", pos=(30.0, 6.5, 65.5), color=(1, 1, 1))

cmd.pseudoatom(object="donor_projected_line_2pa2", pos=(30.75, 4.5, 67.25), color=(1, 1, 1))

cmd.distance(name="donor_projected_line_2", selection1="donor_projected_line_2pa1", selection2="donor_projected_line_2pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_projected_line_2")
cmd.set("dash_width", 4.0)
cmd.delete("donor_projected_line_2pa1")
cmd.delete("donor_projected_line_2pa2")
cmd.group("donor_projected_2", members="donor_projected_point_2_obj")
cmd.group("donor_projected_2", members="donor_projected_point_2_score")
cmd.group("donor_projected_2", members="donor_projected_projection_2")
cmd.group("donor_projected_2", members="donor_projected_line_2")
cmd.group("donor_projected_2", members="donor_projected_point_2_proj_score")
ring_point_3 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 0.8] + [SPHERE, float(31.0), float(8.5), float(65.5), float(2.0)]

cmd.load_cgo(ring_point_3, "ring_point_3_obj", 1)
cmd.pseudoatom(object="ring_point_3_score", pos=(31.0, 8.5, 65.5), color=(1, 1, 1), label=115.1)

cmd.group("ring_3", members="ring_point_3_obj")
cmd.group("ring_3", members="ring_point_3_score")
donor_projected_point_4 = [COLOR, 0.0, 0.0, 1.0] +  [ALPHA, 0.8] + [SPHERE, float(27.0), float(3.0), float(64.5), float(1)]

cmd.load_cgo(donor_projected_point_4, "donor_projected_point_4_obj", 1)
cmd.pseudoatom(object="donor_projected_point_4_score", pos=(27.0, 3.0, 64.5), color=(1, 1, 1), label=102.9)

cmd.pseudoatom(object="donor_projected_point_4_proj_score", pos=(27.25, 0.75, 65.75), color=(1, 1, 1), label=91.2)

donor_projected_projection_4 = [COLOR, 0.6, 0.6, 1.0] +  [ALPHA, 0.8] + [SPHERE, float(27.25), float(0.75), float(65.75), float(1)]

cmd.load_cgo(donor_projected_projection_4, "donor_projected_projection_4_obj", 1)
cmd.pseudoatom(object="donor_projected_line_4pa1", pos=(27.0, 3.0, 64.5), color=(1, 1, 1))

cmd.pseudoatom(object="donor_projected_line_4pa2", pos=(27.25, 0.75, 65.75), color=(1, 1, 1))

cmd.distance(name="donor_projected_line_4", selection1="donor_projected_line_4pa1", selection2="donor_projected_line_4pa2", width=0.5, gap=0.2, label=0, state=1)

cmd.set("dash_color", (0.6, 0.6, 1.0), selection="donor_projected_line_4")
cmd.set("dash_width", 4.0)
cmd.delete("donor_projected_line_4pa1")
cmd.delete("donor_projected_line_4pa2")
cmd.group("donor_projected_4", members="donor_projected_point_4_obj")
cmd.group("donor_projected_4", members="donor_projected_point_4_score")
cmd.group("donor_projected_4", members="donor_projected_projection_4")
cmd.group("donor_projected_4", members="donor_projected_line_4")
cmd.group("donor_projected_4", members="donor_projected_point_4_proj_score")
ring_point_5 = [COLOR, 0.33, 1.0, 0.0] +  [ALPHA, 0.8] + [SPHERE, float(25.0), float(7.0), float(61.5), float(2.0)]

cmd.load_cgo(ring_point_5, "ring_point_5_obj", 1)
cmd.pseudoatom(object="ring_point_5_score", pos=(25.0, 7.0, 61.5), color=(1, 1, 1), label=62.2)

cmd.group("ring_5", members="ring_point_5_obj")
cmd.group("ring_5", members="ring_point_5_score")
cmd.group("donor_projected_pts", members="donor_projected_2")
cmd.group("donor_projected_pts", members="donor_projected_4")
cmd.group("acceptor_projected_pts", members="acceptor_projected_1")
cmd.group("ring_pts", members="ring_0")
cmd.group("ring_pts", members="ring_3")
cmd.group("ring_pts", members="ring_5")
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
cmd.set_color("donor_projected_color", (0.0, 0.0, 1.0))
cmd.load("donor_projected.grd", "donor_projected_grid")
cmd.isosurface(name="surface_donor_projected", map="donor_projected_grid", level="1")

cmd.color("donor_projected_color", "surface_donor_projected")
cmd.set_color("acceptor_projected_color", (1.0, 0.0, 0.0))
cmd.load("acceptor_projected.grd", "acceptor_projected_grid")
cmd.isosurface(name="surface_acceptor_projected", map="acceptor_projected_grid", level="1")

cmd.color("acceptor_projected_color", "surface_acceptor_projected")
cmd.set_color("ring_color", (0.33, 1.0, 0.0))
cmd.load("ring.grd", "ring_grid")
cmd.isosurface(name="surface_ring", map="ring_grid", level="1")

cmd.color("ring_color", "surface_ring")
cmd.group("feature_grids", members="donor_projected")
cmd.group("feature_grids", members="acceptor_projected")
cmd.group("feature_grids", members="ring")
cmd.group("feature_grids", members="surface_donor_projected")
cmd.group("feature_grids", members="surface_acceptor_projected")
cmd.group("feature_grids", members="surface_ring")


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


surface_list = {None: {'feature_grids': ['surface_donor_projected', 'surface_acceptor_projected', 'surface_ring']}}
surface_max_list = {None: {'feature_grids': 191.4}}

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



