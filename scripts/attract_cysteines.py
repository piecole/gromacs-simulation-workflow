"""
Move two cysteine side chains so their SG atoms sit a target distance apart
(default 2.1 A, the typical S-S disulfide distance).

Strategy
--------
1. Rotation first: rotate each SG around its own CA-CB axis (a chi1 rotation)
   to bring the two SG atoms as close to the target distance as rotation alone
   allows. The CA-CB bond defines the rotation axis.
2. Bond stretch only for the residual: if rotation cannot reach the target, the
   remaining gap is closed by sliding each SG along the SG-SG line toward the
   other, with CB sharing part of the shift so the elongation is spread across
   the CA-CB and CB-SG bonds.

Usage
-----
    python attract_cysteines.py input.gro output.gro ca1 cb1 sg1 ca2 cb2 sg2

Atom numbers are 1-based GROMACS atom numbers, e.g.

    python attract_cysteines.py input.gro output.gro 110 111 112 211 212 113

Notes
-----
.gro coordinates are in nm, so 2.1 A = 0.21 nm. Periodic boundary conditions
are ignored; this is intended for two cysteines already in the same image.
"""

import argparse

import numpy as np


def read_gro_simple(gro_file):
    with open(gro_file) as f:
        lines = f.readlines()

    title = lines[0].rstrip("\n")
    natoms = int(lines[1].strip())
    atom_lines = lines[2:2 + natoms]
    box_line = lines[2 + natoms].rstrip("\n")

    atoms = []

    for line in atom_lines:
        atoms.append({
            "resid": line[0:5],
            "resname": line[5:10],
            "atomname": line[10:15],
            "atomnr": line[15:20],
            "coord": np.array([
                float(line[20:28]),
                float(line[28:36]),
                float(line[36:44]),
            ], dtype=float),
            "rest": line[44:].rstrip("\n"),
        })

    return title, atoms, box_line


def write_gro_simple(out_file, title, atoms, box_line):
    with open(out_file, "w") as f:
        f.write(title + "\n")
        f.write(f"{len(atoms):5d}\n")

        for atom in atoms:
            x, y, z = atom["coord"]
            f.write(
                f"{atom['resid']}"
                f"{atom['resname']}"
                f"{atom['atomname']}"
                f"{atom['atomnr']}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{atom['rest']}\n"
            )

        f.write(box_line + "\n")


def rotate_point_around_axis(point, axis_point_1, axis_point_2, angle_rad):
    """
    Rotate ``point`` around the axis defined by axis_point_1 -> axis_point_2
    using Rodrigues' rotation formula.
    """
    point = np.asarray(point, dtype=float)
    a = np.asarray(axis_point_1, dtype=float)
    b = np.asarray(axis_point_2, dtype=float)

    axis = b - a
    norm = np.linalg.norm(axis)
    if norm == 0:
        raise ValueError("Rotation axis has zero length (CA and CB overlap).")
    axis = axis / norm

    p = point - a

    p_rot = (
        p * np.cos(angle_rad)
        + np.cross(axis, p) * np.sin(angle_rad)
        + axis * np.dot(axis, p) * (1 - np.cos(angle_rad))
    )

    return a + p_rot


def _best_rotation(sg1_start, ca1, cb1, sg2_start, ca2, cb2, target_nm,
                   coarse_step_deg):
    """
    Find the (angle1, angle2) chi1 rotations that bring SG1-SG2 as close to
    ``target_nm`` as possible. Objective is |distance - target|, so the search
    reaches the target exactly when rotation allows, and otherwise lands at the
    closest approach. A coarse grid is followed by coordinate-descent refinement.

    Returns (score, distance, angle1, angle2, sg1_pos, sg2_pos).
    """
    def evaluate(a1, a2):
        s1 = rotate_point_around_axis(sg1_start, ca1, cb1, a1)
        s2 = rotate_point_around_axis(sg2_start, ca2, cb2, a2)
        d = np.linalg.norm(s2 - s1)
        return abs(d - target_nm), d, a1, a2, s1, s2

    angles = np.deg2rad(np.arange(0.0, 360.0, coarse_step_deg))

    best = None
    for a1 in angles:
        for a2 in angles:
            cand = evaluate(a1, a2)
            if best is None or cand[0] < best[0]:
                best = cand

    # Coordinate-descent refinement with a shrinking step.
    step = np.deg2rad(coarse_step_deg)
    min_step = np.deg2rad(1e-4)
    while step > min_step:
        improved = False
        for delta in (step, -step):
            for which in (1, 2):
                a1 = best[2] + (delta if which == 1 else 0.0)
                a2 = best[3] + (delta if which == 2 else 0.0)
                cand = evaluate(a1, a2)
                if cand[0] < best[0] - 1e-15:
                    best = cand
                    improved = True
        if not improved:
            step *= 0.5

    return best


def attract_cysteines(
    gro_file,
    out_file,
    ca1,
    cb1,
    sg1,
    ca2,
    cb2,
    sg2,
    target_angstrom=2.1,
    angle_step_degrees=2.0,
    cb_fraction=0.5,
):
    """
    Bring two cysteine SG atoms to ``target_angstrom`` apart, prioritising chi1
    rotation around each CA-CB axis and stretching the side-chain bonds only for
    the residual distance.

    Parameters
    ----------
    gro_file, out_file : str
        Input and output .gro paths.
    ca1, cb1, sg1, ca2, cb2, sg2 : int
        1-based GROMACS atom numbers for the two cysteines.
    target_angstrom : float
        Target SG-SG distance in Angstrom (default 2.1).
    angle_step_degrees : float
        Coarse grid step for the rotation search; refinement runs afterwards.
    cb_fraction : float
        Fraction of each SG's stretch shift that CB also moves (CA stays fixed).
        SG-SG distance depends only on the SG atoms, so this does not change the
        final separation; it shares the residual elongation between the CA-CB
        and CB-SG bonds instead of dumping it all on CB-SG. 0.0 keeps CB fixed
        (all distortion on CB-SG); 0.5 (default) splits it roughly evenly.
    """
    title, atoms, box_line = read_gro_simple(gro_file)

    # Convert 1-based GROMACS atom numbers to 0-based list positions.
    indices = {
        "ca1": ca1 - 1, "cb1": cb1 - 1, "sg1": sg1 - 1,
        "ca2": ca2 - 1, "cb2": cb2 - 1, "sg2": sg2 - 1,
    }
    for name, idx in indices.items():
        if idx < 0 or idx >= len(atoms):
            raise IndexError(
                f"Atom number for {name} is out of range "
                f"(got {idx + 1}, file has {len(atoms)} atoms)."
            )

    target_nm = target_angstrom / 10.0

    ca1_pos = atoms[indices["ca1"]]["coord"].copy()
    cb1_pos = atoms[indices["cb1"]]["coord"].copy()
    sg1_start = atoms[indices["sg1"]]["coord"].copy()

    ca2_pos = atoms[indices["ca2"]]["coord"].copy()
    cb2_pos = atoms[indices["cb2"]]["coord"].copy()
    sg2_start = atoms[indices["sg2"]]["coord"].copy()

    initial_dist = np.linalg.norm(sg2_start - sg1_start)
    cbsg1_start = np.linalg.norm(sg1_start - cb1_pos)
    cbsg2_start = np.linalg.norm(sg2_start - cb2_pos)
    cacb1_start = np.linalg.norm(cb1_pos - ca1_pos)
    cacb2_start = np.linalg.norm(cb2_pos - ca2_pos)

    # --- Step 1: rotation -------------------------------------------------
    _, dist_after_rotation, angle1, angle2, sg1_rot, sg2_rot = _best_rotation(
        sg1_start, ca1_pos, cb1_pos,
        sg2_start, ca2_pos, cb2_pos,
        target_nm, angle_step_degrees,
    )

    atoms[indices["sg1"]]["coord"] = sg1_rot.copy()
    atoms[indices["sg2"]]["coord"] = sg2_rot.copy()

    # --- Step 2: bond stretch for the residual ----------------------------
    # Slide each SG along the SG-SG line so the final separation is exactly the
    # target, and move CB by ``cb_fraction`` of that shift so the elongation is
    # shared between the CA-CB and CB-SG bonds rather than dumped on CB-SG. CA is
    # fixed. SG-SG distance depends only on the SG atoms, so moving CB does not
    # affect the final separation. When rotation already reached the target this
    # move is negligible.
    vec = sg2_rot - sg1_rot
    dist = np.linalg.norm(vec)
    if dist == 0:
        raise ValueError(
            "SG atoms overlap after rotation; cannot define a stretch direction."
        )
    unit = vec / dist

    excess = dist - target_nm
    stretch_applied = abs(excess) > 1e-9

    sg1_shift = 0.5 * excess * unit
    sg2_shift = -0.5 * excess * unit

    atoms[indices["cb1"]]["coord"] = cb1_pos + cb_fraction * sg1_shift
    atoms[indices["cb2"]]["coord"] = cb2_pos + cb_fraction * sg2_shift
    atoms[indices["sg1"]]["coord"] = sg1_rot + sg1_shift
    atoms[indices["sg2"]]["coord"] = sg2_rot + sg2_shift

    ca1_final = atoms[indices["ca1"]]["coord"]
    cb1_final = atoms[indices["cb1"]]["coord"]
    sg1_final = atoms[indices["sg1"]]["coord"]
    ca2_final = atoms[indices["ca2"]]["coord"]
    cb2_final = atoms[indices["cb2"]]["coord"]
    sg2_final = atoms[indices["sg2"]]["coord"]

    final_dist = np.linalg.norm(sg2_final - sg1_final)
    cbsg1_final = np.linalg.norm(sg1_final - cb1_final)
    cbsg2_final = np.linalg.norm(sg2_final - cb2_final)
    cacb1_final = np.linalg.norm(cb1_final - ca1_final)
    cacb2_final = np.linalg.norm(cb2_final - ca2_final)

    write_gro_simple(out_file, title, atoms, box_line)

    print(f"Target SG-SG distance:         {target_nm:.4f} nm = {target_angstrom:.2f} A")
    print(f"Initial SG-SG distance:        {initial_dist:.4f} nm = {initial_dist * 10:.2f} A")
    print(f"After rotation SG-SG distance: {dist_after_rotation:.4f} nm = {dist_after_rotation * 10:.2f} A")
    print(f"Final SG-SG distance:          {final_dist:.4f} nm = {final_dist * 10:.2f} A")
    print(f"Chi1 rotation cysteine 1:      {np.rad2deg(angle1) % 360:.2f} degrees")
    print(f"Chi1 rotation cysteine 2:      {np.rad2deg(angle2) % 360:.2f} degrees")
    if stretch_applied:
        print(f"Bond stretch needed (rotation alone could not reach target), cb_fraction={cb_fraction}:")
        print(f"  CA-CB bond 1: {cacb1_start * 10:.3f} A -> {cacb1_final * 10:.3f} A")
        print(f"  CB-SG bond 1: {cbsg1_start * 10:.3f} A -> {cbsg1_final * 10:.3f} A")
        print(f"  CA-CB bond 2: {cacb2_start * 10:.3f} A -> {cacb2_final * 10:.3f} A")
        print(f"  CB-SG bond 2: {cbsg2_start * 10:.3f} A -> {cbsg2_final * 10:.3f} A")
    else:
        print("Rotation alone reached the target; no bond stretch needed.")
    print(f"Wrote: {out_file}")

    return final_dist


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Move two cysteine side chains so their SG atoms are a target "
            "distance apart, prioritising chi1 rotation then bond stretching."
        )
    )
    parser.add_argument("input", help="Input .gro file")
    parser.add_argument("output", help="Output .gro file")
    parser.add_argument("ca1", type=int, help="1-based atom number of CA, cysteine 1")
    parser.add_argument("cb1", type=int, help="1-based atom number of CB, cysteine 1")
    parser.add_argument("sg1", type=int, help="1-based atom number of SG, cysteine 1")
    parser.add_argument("ca2", type=int, help="1-based atom number of CA, cysteine 2")
    parser.add_argument("cb2", type=int, help="1-based atom number of CB, cysteine 2")
    parser.add_argument("sg2", type=int, help="1-based atom number of SG, cysteine 2")
    parser.add_argument(
        "--target", type=float, default=2.1,
        help="Target SG-SG distance in Angstrom (default: 2.1)",
    )
    parser.add_argument(
        "--angle-step", type=float, default=2.0,
        help="Coarse rotation grid step in degrees (default: 2.0)",
    )
    parser.add_argument(
        "--cb-fraction", type=float, default=0.5,
        help=(
            "Fraction of the SG stretch that CB also moves, sharing the "
            "elongation between CA-CB and CB-SG bonds (default: 0.5; "
            "0.0 keeps CB fixed)"
        ),
    )

    args = parser.parse_args(argv)

    attract_cysteines(
        args.input,
        args.output,
        args.ca1,
        args.cb1,
        args.sg1,
        args.ca2,
        args.cb2,
        args.sg2,
        target_angstrom=args.target,
        angle_step_degrees=args.angle_step,
        cb_fraction=args.cb_fraction,
    )


if __name__ == "__main__":
    main()
