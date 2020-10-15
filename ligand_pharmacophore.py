import argparse
import os
from ccdc.descriptors import GeometricDescriptors
from ccdc import io
from hotspots.pharmacophore_extension import LigandPharmacophoreModel, create_consensus
from hotspots.result import Results
from hotspots.hs_io import HotspotWriter
from ligand_overlay import ftp_download
import tempfile
from ccdc.protein import Protein


class Runner(argparse.ArgumentParser):
    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)

        self.add_argument(
            'ligand_overlay',
            help='Ligand overlay'
        )
        self.add_argument(
            'ref_pdb',
            help='Reference PDB'
        )
        self.add_argument(
            '-o', '--output-directory', default='.',
            help='Where output will be stored'
        )

        self.args = self.parse_args()

    @staticmethod
    def generate_pharmacophore(ligands, ref_pdb, out_dir):
        lig_pharms = []
        for ligand in ligands:
            ligand_pharmacophore = LigandPharmacophoreModel()
            ligand_pharmacophore.feature_definitions = ["ring",
                                                        "acceptor_projected",
                                                        "donor_projected"]

            ligand_pharmacophore.detect_from_ligand(ligand)

            for feat in ligand_pharmacophore.detected_features:
                ligand_pharmacophore.add_feature(feat)

            lig_pharms.append(ligand_pharmacophore)

        # 20 %
        cutoff = len(ligands) * 0.2
        feats, feat_point_grds = create_consensus(lig_pharms, cutoff=cutoff)
        print(feats)
        for feat in feats:
            if feat.identifier == "ring":
                p = feat.spheres[0].centre
                feat.spheres = (GeometricDescriptors.Sphere((p[0], p[1], p[2]), 2.0),)
                feat.point = feat.spheres[0]

        ensemble_pharm = LigandPharmacophoreModel()
        ensemble_pharm.detected_features = feats
        ensemble_pharm.feature_point_grids = feat_point_grds
        ensemble_pharm.ligands = ligands
        ensemble_pharm.detected_features = ensemble_pharm.top_features(num=6)
        pymol_o = os.path.join(out_dir, "pymol")
        if not os.path.exists(pymol_o):
            os.mkdir(pymol_o)
        ensemble_pharm.pymol_visulisation(pymol_o)

        #  enable rescoring
        tmp = tempfile.mkdtemp()
        ftp_download([ref_pdb, tmp])
        hr = Results(super_grids={"apolar": feat_point_grds["ring"],
                                  "donor": feat_point_grds["donor_projected"],
                                  "acceptor": feat_point_grds["acceptor_projected"]},
                     protein=Protein.from_file(os.path.join(tmp, f"{ref_pdb}.pdb")))

        hr_out = os.path.join(out_dir, "hr")
        if not os.path.exists(hr_out):
            os.mkdir(hr_out)
        with HotspotWriter(hr_out) as w:
            w.write(hr)

        p_out = os.path.join(out_dir, "ligand_pharmacophores")
        if not os.path.exists(p_out):
            os.mkdir(p_out)

        for n in [6, 5, 4, 3]:
            lp = LigandPharmacophoreModel()
            lp.detected_features = feats
            lp.detected_features = lp.top_features(num=n)
            for feat in lp.detected_features:
                lp.add_feature(feat)

            lp.intra_only = True

            lp.write(os.path.join(p_out, f"{n}.cm"))

    def run(self):
        ligands = io.CrystalReader(self.args.ligand_overlay)
        self.generate_pharmacophore(ligands,
                                    self.args.ref_pdb,
                                    self.args.output_directory)


if __name__ == "__main__":
    runner = Runner()
    runner.run()