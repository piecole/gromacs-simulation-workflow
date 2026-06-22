"""
Insert TER records at chain breaks in a PDB file.

GROMACS .gro files carry no chain information, so a .gro -> .pdb round trip
fuses all chains into one. Re-running pdb2gmx on such a file fails because the
C-terminal carboxylate oxygens (OT1/OT2/OXT) of an internal chain look like
illegal extra atoms on a mid-chain residue, e.g.

    Atom OT1 in residue GLY 213 was not found in rtp entry GLY ...

This script restores the chain breaks. Every chain produced by a prior
CHARMM pdb2gmx run is capped COO-, so each chain end is the residue carrying
the terminal oxygens. We drop any existing TER lines and emit a fresh TER
after each residue that contains a terminal oxygen.

Usage
-----
    python insert_ter.py input.pdb output.pdb

Then feed output.pdb to pdb2gmx (with -ter -ignh -ss, plus -merge all for an
inter-chain disulfide).
"""

import argparse

TERMINAL_OXYGENS = {"OT1", "OT2", "OXT", "O1", "O2"}


def residue_key(line):
    # chain id (22), resSeq (23-26), iCode (27) in 1-based PDB columns
    return line[21:27]


def atom_name(line):
    return line[12:16].strip()


def insert_ter(in_path, out_path):
    with open(in_path) as handle:
        lines = handle.readlines()

    out = []
    current_key = None
    current_is_terminal = False
    chain_breaks = 0

    def flush_ter():
        nonlocal chain_breaks
        if current_is_terminal:
            out.append("TER\n")
            chain_breaks += 1

    for line in lines:
        record = line[:6].strip()

        if record == "TER":
            # Drop pre-existing TER records; we regenerate them.
            continue

        if record in ("ATOM", "HETATM"):
            key = residue_key(line)
            if current_key is None:
                current_key = key
            elif key != current_key:
                # Residue boundary: close the previous residue.
                flush_ter()
                current_key = key
                current_is_terminal = False

            if atom_name(line) in TERMINAL_OXYGENS:
                current_is_terminal = True
            out.append(line)
        else:
            # Non-atom record (END, CRYST1, etc.): close any open chain first.
            if record == "END":
                flush_ter()
                current_key = None
                current_is_terminal = False
            out.append(line)

    # Close a trailing chain if the file ended on atoms with no END.
    if current_key is not None:
        flush_ter()

    with open(out_path, "w") as handle:
        handle.writelines(out)

    print(f"Inserted {chain_breaks} TER record(s) -> {out_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="input PDB file")
    parser.add_argument("output", help="output PDB file with TER records")
    args = parser.parse_args()
    insert_ter(args.input, args.output)


if __name__ == "__main__":
    main()
