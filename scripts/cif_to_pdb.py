import gemmi
import argparse

all_characters = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert CIF file to PDB file")
    parser.add_argument("cif_input", type=str, help="Path to the CIF file")
    parser.add_argument("pdb_output", type=str, help="Path to the PDB file")
    args = parser.parse_args()

    cif = gemmi.read_cif(args.cif_input)

    for model in cif.models:
        chains = []
        for chain in model.chains:
            chains.append(chain.name)
        if any(len(chain) > 1 for chain in chains):
            available_characters = [char for char in all_characters if char not in chains]
            for chain in chains:
                if len(chain) > 1:
                    chain.name = available_characters.pop(0)

    pdb = gemmi.write_pdb(args.pdb_output)
    print(f"Converted {args.cif_input} to {args.pdb_output}")