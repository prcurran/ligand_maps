import os
from multiprocessing import Pool

from tqdm import tqdm

from pdb_superimposer.data_downloader import pdb_search_query, ftp_download

import pandas as pd


def wrap_ftp_download(inputs):
    #  simple wrapper to manage flow of args
    pdb, out_dir = inputs
    return ftp_download(pdb, out_dir)


def check_dir(d):
    if not os.path.exists(d):
        os.mkdir(d)
    return d


def main():
    ########################################################
    search_query_path = "example_search.json"
    out_dir = "pdb"
    ref_id = "2VTA"
    polymer_entity = 0  # PDB speak for chain
    processes = 6
    ########################################################

    check_dir(out_dir)

    # Step 1: Run the UniProt Search
    print("Searching PDB")
    with open(search_query_path, "r") as r:
        query = r.read()  # keep query as str

    results = pdb_search_query(query)

    pdb_entities = {i["identifier"].split("_")[0]: int(i["identifier"].split("_")[1]) - 1
                    for i in results["result_set"]}
    print(pdb_entities)

    # for testing
    if ref_id not in pdb_entities:
        pdb_entities.update({ref_id: polymer_entity})

    # Step 2: Download PDBs
    print("Downloading data")
    args = ((a, out_dir) for a in pdb_entities.keys())
    with Pool(processes=processes) as pool:
        list(tqdm(pool.imap_unordered(wrap_ftp_download, args), total=len(pdb_entities)))

    df = pd.DataFrame({"pdb": list(pdb_entities.keys()), "chain": list(pdb_entities.values())})
    df.to_csv(os.path.join(out_dir, "chain_info.csv"))


if __name__ == "__main__":
    main()

