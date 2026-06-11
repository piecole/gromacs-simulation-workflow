"""
Download a PDB file from the RCSB PDB and convert it to a PDB file.
"""

import os
import argparse

import gemmi
import requests

all_characters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")

def download_pdb(pdb_id):
    """
    Download a PDB file from the RCSB PDB and return the path to the downloaded file.
    """
    url = f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif.gz"
    response = requests.get(url, timeout=50)
    with open(f"{pdb_id}-assembly1.cif.gz", "wb") as f:
        f.write(response.content)
    return f"{pdb_id}-assembly1.cif.gz"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert CIF file to PDB file")
    parser.add_argument("pdb_id", type=str, help="PDB ID")
    args = parser.parse_args()

    download_pdb(args.pdb_id)
    cif = gemmi.read_structure(f"{args.pdb_id}-assembly1.cif.gz")

    for model in cif:
        chains = []
        for chain in model:
            chains.append(chain)
        if any(len(chain) > 1 for chain in chains):
            available_characters = [char for char in all_characters if char not in chains]
            for chain in chains:
                if len(chain) > 1:
                    chain.name = available_characters.pop(0)

    pdb = cif.write_pdb(f"{args.pdb_id}-assembly.pdb")

    os.remove(f"{args.pdb_id}-assembly1.cif.gz")
    print(f"Downloaded {args.pdb_id}")