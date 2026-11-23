#!/bin/bash

# Cluster provides the MPI build as gmx_mpi (not thread-MPI gmx).
# Ensure the gromacs module is loaded before running this script.
GMX=gmx_mpi
if ! command -v "$GMX" >/dev/null 2>&1; then
    echo "Error: '$GMX' not found. Load the gromacs module or adjust GMX."
    exit 1
fi

# Loop through directories matching the pattern
for dir in run_pepd_p53_restraint_*/; do
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"
        cd "$dir"

        # Run your GMX command here
        "$GMX" distance -f md_0_1.xtc -s md_0_1.tpr \
            -select "atomnr 903 plus atomnr 10001; atomnr 8466 plus atomnr 2438" \
            -oall distance.xvg

        # Return to parent directory
        cd ..
    fi
done
