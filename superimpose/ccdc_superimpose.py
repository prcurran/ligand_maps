
import os
from ccdc.protein import Protein
from ccdc.io import MoleculeWriter, MoleculeReader
import pandas as pd


def align_proteins(pdbs, chains, hetids, input_dir, output_dir):
    # read all PDB's
    prots = [Protein.from_file(os.path.join(input_dir, f"{pdb}.pdb")) for pdb in pdbs]

    for p in prots:
        p.detect_ligand_bonds()

    ref_prot = Protein.from_file(os.path.join(input_dir, "2VTA.pdb"))
    ref_prot.detect_ligand_bonds()
    ref_mol = [l for l in ref_prot.ligands if l.identifier.split(":")[1][:3] == "LZ1"][0]

    with MoleculeWriter(os.path.join(output_dir, f"{ref_prot.identifier}_aligned.pdb")) as w:
        w.write(ref_prot)

    with MoleculeWriter(os.path.join(output_dir, f"{ref_prot.identifier}_{ref_mol.identifier.split(':')[1][:3]}_aligned.mol2")) as w:
        w.write(ref_mol)

    # create ref binding site
    ref_bind_site = Protein.BindingSiteFromMolecule(protein=ref_prot,
                                                    molecule=ref_mol,
                                                    distance=12)

    for p, h in zip(prots, hetids):
        try:
            print(p.identifier)
            # align prot
            chain_superposition = Protein.ChainSuperposition()
            (rmsd, transformation) = chain_superposition.superpose(ref_prot.chains[0],
                                                                   p.chains[0],
                                                                   ref_bind_site)
            p_name = p.identifier
            # get aligned mol
            p.detect_ligand_bonds()

            a_mol = [lig for lig in p.ligands if lig.identifier.split(":")[1][:3] == h[:3]][0]

            # write protein
            with MoleculeWriter(os.path.join(output_dir, f"prot_{p_name}_aligned.pdb")) as w:
                w.write(p)

            # write ligand
            with MoleculeWriter(os.path.join(output_dir,
                                             f"mol_{p_name}_{a_mol.identifier.split(':')[1][:3]}_aligned.mol2")) as w:
                w.write(a_mol)
        except:
            pass


if __name__ == "__main__":
    # test
    parent = os.path.dirname(os.path.dirname(__file__))
    df = pd.read_csv(os.path.join(parent, "cluster_ligands", "selected.csv"))

    align_proteins(list(df.pdb),
                   list(df.chain),
                   list(df.hetid),
                   os.path.join(parent, "download", "pdb"),
                   "pdb")
