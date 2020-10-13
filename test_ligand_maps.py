import unittest
from ccdc import io
from ligand_overlay import *
from rdkit.Chem import Mol
import os


class TestMethods(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = "testdata"

    def testpdb_search(self):
        uniprot = "P24941"
        max_resolution = 2.5
        pdbs, entities = pdb_search(uniprot, max_resolution)
        self.assertIsInstance(pdbs, list)

        pdb = pdbs[0]
        self.assertIsInstance(pdb, str)
        self.assertEqual(4, len(pdb))

        entity = entities[0]
        self.assertIsInstance(entity, int)

    def testftp_download(self):
        ins = ["3COV", self.tmp]
        outs = ftp_download(ins)
        self.assertTrue(os.path.exists(outs))

    def testdownload(self):
        ins = ["2VTA", "1HCL"]

        download(ins, out_dir=self.tmp)

        fs = [os.path.join(self.tmp, f"{i}.pdb") for i in ins]
        for f in fs:
            self.assertTrue(os.path.exists(f))

    def testprepare_protein(self):
        pdb = "3COV"
        entity = 0
        out_dir = self.tmp
        prot = prepare_protein(pdb, entity, out_dir)

        print(prot)

    def testalign(self):
        pdbs = ["2VTA", "1HCL"]
        entities = [0, 0]
        download(pdbs, out_dir=self.tmp)

        for p, e in zip(pdbs, entities):
            prepare_protein(p, e, self.tmp)

        args = ["2VTA", 'LZ11301', "1HCL", self.tmp]
        other, rms = align(args)

        with io.MoleculeWriter(os.path.join(self.tmp, f"aligned_{other.identifier}.pdb")) as w:
            w.write(other)

    # def testalign(self):





    # def testccdc_to_rdkit(self):
    #     mol = io.EntryReader('CSD').molecule('IBPRAC')
    #     rdmol = ccdc_to_rdkit(mol)
    #     self.assertIsInstance(rdmol, Mol)
