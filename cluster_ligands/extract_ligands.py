import os

import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ccdc.protein import Protein
from hotspots.data import common_solvents
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.manifold import TSNE
from ccdc import io
import tempfile


def fingerprint_array(ligands):
    X = []
    for l in ligands:
        arr = np.zeros((0,))
        DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(l, 2), arr)
        X.append(arr)
    return X


def tanimoto_dist(a, b):
    dotprod = np.dot(a, b)
    tc = dotprod / (np.sum(a) + np.sum(b) - dotprod)
    return 1.0 - tc


def cluster_ligands(ligands, t):
    cluster_dic = {}

    # generate fingerprint array
    X = fingerprint_array(ligands)
    if len(X) < 2:
        X = fingerprint_array(ligands)
        if len(X) < 2:
            raise ValueError("Fingerprint array must contain more than 1 entry")

    # dimensionality reduction
    tsne_X = TSNE(n_components=2, metric=tanimoto_dist).fit_transform(np.array(X, dtype=np.float32))

    # clustering
    cluster_tsne = hdbscan.HDBSCAN(min_cluster_size=2, gen_min_span_tree=True)
    cluster_tsne.fit(tsne_X)

    for i, label in enumerate(cluster_tsne.labels_):
        if label == -1:
            continue
        else:
            if label in cluster_dic:
                cluster_dic[label].append(ligands[i])
            else:
                cluster_dic.update({label: [ligands[i]]})

    x = [tsne_X.T[0][j] for j, l in enumerate(cluster_tsne.labels_) if l != -1]
    y = [tsne_X.T[1][j] for j, l in enumerate(cluster_tsne.labels_) if l != -1]
    hue = [l for j, l in enumerate(cluster_tsne.labels_) if l != -1]

    seen = [-1]
    sx = []
    sy = []
    for k, l in enumerate(cluster_tsne.labels_):
        if l in seen:
            continue
        else:
            sx.append(tsne_X.T[0][k])
            sy.append(tsne_X.T[1][k])
            seen.append(l)

    plt.scatter(x, y, c=hue, cmap='RdBu', alpha=0.7)
    plt.scatter(sx, sy, c="black", marker="x")

    plt.title("{} clusters".format(t))
    plt.savefig("{}.png".format(t))
    plt.close()

    if len(cluster_dic) == 0:
        print("NO CLUSTERS FOUND")
        try:
            unique = {}
            for l in ligands:
                hetid = l.chemical_id.split("_")[0]
                if not hetid in unique:
                    unique.update({hetid: l})

            ligands = unique.values()
            cluster_dic = {i: [ligands[i]] for i in range(0, len(ligands))}
        except:
            cluster_dic = {i: [ligands[i]] for i in range(0, len(ligands))}

    return cluster_dic


def ccdc_to_rdkit(mol):
    print(mol.identifier)
    temp = tempfile.mkdtemp()
    temp_path = os.path.join(temp, "mol.mol2")
    with io.MoleculeWriter(temp_path) as w:
        w.write(mol)

    rd_mol = Chem.MolFromMol2File(temp_path)
    rd_mol.SetProp("_Name", mol.identifier)
    return rd_mol


def main():

    pdb_dir = "/home/pcurran/github_packages/ligand_maps/download/pdb"
    df = pd.read_csv(os.path.join(pdb_dir, "chain_info.csv"))

    ligands = []

    for pdb in df.pdb:
        prot = Protein.from_file(os.path.join(pdb_dir, f"{pdb}.pdb"))
        prot.detect_ligand_bonds()

        chain = prot.chains[int(df.loc[df.pdb == pdb].chain)].identifier
        selected = [l for l in prot.ligands if l.identifier.split(":")[0] == chain and
                    l.identifier.split(":")[1][:3] not in common_solvents()]

        if len(selected) == 1:
            selected[0].identifier = f"{pdb}:{selected[0].identifier}"
            try:
                ligands.append(ccdc_to_rdkit(selected[0]))
            except AttributeError:
                pass

    cd = cluster_ligands(ligands, "CDK2")
    selected = [rd_mol[0].GetProp("_Name").split(":") for rd_mol in cd.values()]

    pdbs, chains, hetids = zip(*selected)

    other_df = pd.DataFrame({"pdb":pdbs, "chain": chains, "hetid":hetids})
    other_df.to_csv("selected.csv")


if __name__ == "__main__":
    main()
