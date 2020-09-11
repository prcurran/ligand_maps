from hotspots.pharmacophore_extension import LigandPharmacophoreModel, create_consensus
from ccdc import io
import os


if __name__ == "__main__":
    parent = os.path.dirname(os.path.dirname(__file__))
    mol_dir = os.path.join(parent, "superimpose", "pdb")

    mol_paths = [f for root, d, fs in os.walk(mol_dir) for f in fs if f.split(".")[1] == "mol2"]
    mols = [io.CrystalReader(os.path.join(mol_dir, f))[0] for f in mol_paths]

    lig_pharms = []
    for mol in mols:
        ligand_pharmacophore = LigandPharmacophoreModel()
        ligand_pharmacophore.feature_definitions = ["ring_planar_projected", "acceptor_projected", "donor_projected"]

        ligand_pharmacophore.detect_from_ligand(mol)
        for feat in ligand_pharmacophore.detected_features:
            ligand_pharmacophore.add_feature(feat)

        lig_pharms.append(ligand_pharmacophore)

    feats, feat_point_grds = create_consensus(lig_pharms)

    ensemble_pharm = LigandPharmacophoreModel()
    ensemble_pharm.detected_features = feats
    ensemble_pharm.feature_point_grids = feat_point_grds
    ensemble_pharm.detected_features = ensemble_pharm.top_features(num=5)
    ensemble_pharm.ligands = mols
    ensemble_pharm.pymol_visulisation("")
