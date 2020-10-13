import argparse
import json
import os
import tempfile
from multiprocessing import Pool
from urllib.parse import quote_plus
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from ccdc import io
from ccdc.molecule import Molecule
from ccdc.protein import Protein
from hotspots.pharmacophore_extension import LigandPharmacophoreModel, create_consensus
from hotspots.data import common_solvents
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from tqdm import tqdm


def pdb_search(uniprot, max_resolution):
    """"""
    search_query_path = "template.json"

    with open(search_query_path, "r") as r:
        query_str = r.read()
        # add uniprot code to query
        query_json = json.loads(query_str)
        query_json['query']['nodes'][0]['parameters']['value'] = uniprot
        query_json['query']['nodes'][4]['parameters']['value'] = max_resolution
        query_str = json.dumps(query_json)

    results = http_get_query(query_str)

    # limit for testing
    lim = 10000
    # [pdb, entity - 1]
    pdbs = [i["identifier"].split("_")[0] for i in results["result_set"][:lim]]
    entities = [int(i["identifier"].split("_")[1]) - 1 for i in results["result_set"][:lim]]

    return pdbs, entities


def http_get_query(query):
    """
    """
    # encode the query
    url = f'http://search.rcsb.org/rcsbsearch/v1/query?json={quote_plus(query)}'
    # call the url
    response = urlopen(Request(url))

    if response.status == 200:
        # response received in bytes
        decoded = response.read().decode("utf-8")
        # covert json str to dict
        return json.loads(decoded)
    else:
        print(f"Status:{response.status}\n{response.msg}")


def ftp_download(args):
    """"""
    pdb, out_dir = args

    url = f"https://files.rcsb.org/download/{pdb}.pdb"
    response = urlopen(Request(url))
    f = response.read().decode("utf-8")

    # write out decoded file
    out_path = os.path.join(out_dir, f"{pdb}.pdb")

    with open(out_path, "w") as w:
        w.write(f)

    return out_path


def download(pdbs, out_dir=None, processes=6):
    """"""
    args = zip(pdbs, [out_dir] * len(pdbs))

    with Pool(processes=processes) as pool:
        list(tqdm(pool.imap(ftp_download, args), total=len(pdbs)))


def prepare_protein(pdb, entity, out_dir):
    """"""
    # protein prep
    prot_file = os.path.join(out_dir, f"{pdb}.pdb")
    prot = Protein.from_file(prot_file)
    prot.remove_all_waters()
    prot.remove_all_metals()
    prot.detect_ligand_bonds()
    discard_chains = {c.identifier for c in prot.chains} - {prot.chains[entity].identifier}

    for chain_id in discard_chains:
        prot.remove_chain(chain_id)

    for ligand in prot.ligands:
        if ligand.identifier.split(":")[0] in discard_chains:
            prot.remove_ligand(ligand.identifier)

    # overwrite protein
    with io.MoleculeWriter(prot_file) as w:
        w.write(prot)


def align(inputs):
    """"""
    try:
        ref_pdb, ref_mol, other_pdb, input_dir = inputs

        ref = Protein.from_file(os.path.join(input_dir, f"{ref_pdb}.pdb"))

        other_path = os.path.join(input_dir, f"{other_pdb}.pdb")
        other = Protein.from_file(other_path)

        if ref_mol:
            ref_mol_obj = [lig for lig in ref.ligands if lig.identifier.split(":")[1] == ref_mol][0]
            ref_bind_site = Protein.BindingSiteFromMolecule(protein=ref, molecule=ref_mol_obj, distance=12)
        else:
            ref_bind_site = None

        chain_superposition = Protein.ChainSuperposition()
        # other chains already striped
        rms, X = chain_superposition.superpose(ref.chains[0], other.chains[0], binding_site1=ref_bind_site)

        with io.MoleculeWriter(other_path) as w:
            w.write(other)

        return rms

    except:
        return 999


def ccdc_to_rdkit(mol):
    """
    Convert CCDC mol to RDKit mol
    :param mol:
    :return:
    """
    rd_mol = Chem.MolFromMolBlock(mol.to_string("mol"))
    rd_mol.SetProp("_Name", mol.identifier)
    return rd_mol


def rdkit_to_ccdc(mol):
    """
    Convert RDKit mol to CCDC mol
    :param mol:
    :return:
    """
    ccdc_mol = Molecule.from_string(Chem.MolToMolBlock(mol), format='mol')
    ccdc_mol.identifier = mol.GetProp("_Name")
    return ccdc_mol


def ligands_from_pdb(pdb, input_dir):
    """
    """
    ligands = []
    prot = Protein.from_file(os.path.join(input_dir, f"{pdb}.pdb"))
    prot.remove_all_metals()
    prot.detect_ligand_bonds()

    extracted = [l for l in prot.ligands]

    # may be more than one ligand per chain
    for extract in extracted:
        # change identifier so ligand can be tracked
        # identifier {pdb}:{chain}:{hetid}{resnum}
        extract.identifier = f"{pdb}:{extract.identifier}"
        ligands.append(extract)

    return ligands


def bs_filter(ligands, ref_mol, cutoff=6):
    all_ligands = []
    for lig in ligands:
        try:
            a = lig.centre_of_geometry()
            b = ref_mol.centre_of_geometry()
            dist = np.linalg.norm(np.array([a.x, a.y, a.z]) - np.array([b.x, b.y, b.z]))
            if dist < cutoff:
                all_ligands.append(lig)
        except:
            pass

    return all_ligands


def remove_solvents(ligands, min_atoms=5):
    """
    """
    return [lig for lig in ligands if all((lig.identifier.split(":")[2][:3] not in common_solvents(),
                                           int(len(lig.heavy_atoms)) >= min_atoms))]


def cluster_ligands(ligands, cutoff=0.2):
    """"""
    rdkit_ligands = []
    for lig in ligands:
        try:
            rdkit_ligands.append(ccdc_to_rdkit(lig))
        except:
            pass

    # from RDKit Cookbook
    fps = [AllChem.GetMorganFingerprintAsBitVect(lig, 2, 1024) for lig in rdkit_ligands]
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])

    # now cluster the data:
    clusters = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    all_ligands = []
    for cluster in clusters:
        try:
            all_ligands.append(rdkit_to_ccdc(rdkit_ligands[cluster[0]]))
        except:
            pass
    return all_ligands


class Runner(argparse.ArgumentParser):
    def __init__(self):
        super(self.__class__, self).__init__(description=__doc__)
        self.add_argument(
            'uniprot',
            help='A UNIPROT accession identifier for the overlay'
        )
        self.add_argument(
            'max_resolution', type=float,
            help='Any PDB structure with resolution greater than this value will be discarded'
        )
        self.add_argument(
            'ref_pdb',
            help='Reference PDB structure for alignment'
        )
        self.add_argument(
            'ref_chain',
            help='Reference Chain for alignment'
        )
        self.add_argument(
            '-o', '--output-directory', default='.',
            help='Where output will be stored'
        )
        self.add_argument(
            '-b', '--binding-site-ligand', default=None,
            help='Reference Binding-Site ligand'
        )
        self.add_argument(
            '-s', '--no-solvents', default=True,
            help='Include solvent in overlay?'
        )
        self.add_argument(
            '-c', '--cluster', default=True,
            help='Cluster ligands on structural similarity?'
        )
        self.add_argument(
            '-p', '--processes', type=int, default=6,
            help='Allocation of CPUs'
        )
        self.args = self.parse_args()

    def run(self):
        rms_cutoff = 1

        tmp = tempfile.mkdtemp()

        print("Searching the PDB")
        pdbs, entities = pdb_search(self.args.uniprot, self.args.max_resolution)

        download(pdbs, out_dir=tmp, processes=self.args.processes)

        if self.args.ref_pdb not in pdbs:
            print("Downloading reference data")
            download([self.args.ref_pdb], out_dir=tmp, processes=1)

        # chain to entity
        ref_prot = Protein.from_file(os.path.join(tmp, f"{self.args.ref_pdb}.pdb"))
        entity = [c.identifier for c in ref_prot.chains].index(self.args.ref_chain)
        prepare_protein(pdb=self.args.ref_pdb, entity=entity, out_dir=tmp)

        print("Prepare protein")
        for pdb, entity in tqdm(zip(pdbs, entities), total=len(pdbs)):
            # other chains removed, max number chains = 1
            try:
                # list(filter(None.__ne__, L))
                prepare_protein(pdb=pdb, entity=entity, out_dir=tmp)
            except:
                print("ERROR", pdb, entity)
                pdbs.remove(pdb)
                entities.remove(entity)

        print("Aligning PDBs")
        # align
        args = zip([self.args.ref_pdb] * len(pdbs),
                   [self.args.binding_site_ligand] * len(pdbs),
                   pdbs,
                   [tmp] * len(pdbs))

        with Pool(processes=self.args.processes) as pool:
            rms_list = list(tqdm(pool.imap(align, args), total=len(pdbs)))

        print("Extracting Ligands...")
        all_ligands = []
        # add the reference structure back to PDB list
        for pdb, rms in zip(pdbs + [self.args.ref_pdb], rms_list + [0]):
            if rms < rms_cutoff:
                all_ligands.extend(ligands_from_pdb(pdb, tmp))
        print(f"    {len(all_ligands)} detected")

        if self.args.binding_site_ligand:
            ref_prot = Protein.from_file(os.path.join(tmp, f"{self.args.ref_pdb}.pdb"))
            ref_mol = [lig for lig in ref_prot.ligands if
                       lig.identifier == f"{self.args.ref_chain}:{self.args.binding_site_ligand}"][0]
            all_ligands = bs_filter(all_ligands, ref_mol)
            print(f"    Binding site filter, {len(all_ligands)} ligands remaining")

        if self.args.no_solvents:
            all_ligands = remove_solvents(all_ligands)
            print(f"    Removed solvents and metals, {len(all_ligands)} ligands remaining")

        if self.args.cluster:
            all_ligands = cluster_ligands(all_ligands)
            print(f"    Clustered ligands, {len(all_ligands)} ligands remaining")

        pdbs, chains, hetids, resnum = zip(*[[lig.identifier.split(":")[0],
                                              lig.identifier.split(":")[1],
                                              lig.identifier.split(":")[2][:3],
                                              lig.identifier.split(":")[2][3:]]
                                             for lig in all_ligands])

        df = pd.DataFrame({"pdbs": pdbs, "chains": chains, "hetids": hetids, "resnum": resnum})
        df.to_csv(os.path.join(self.args.output_directory, "ligand_overlay.csv"))

        with io.MoleculeWriter(os.path.join(self.args.output_directory, "ligand_overlay.mol2")) as w:
            for lig in all_ligands:
                lig.add_hydrogens()
                w.write(lig)

if __name__ == '__main__':
    r = Runner()
    r.run()
