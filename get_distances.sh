#!/bin/bash

# Loop through directories matching the pattern
for dir in run_pepd_p53_restraint_*/; do
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"
        cd "$dir"

        # Run your GMX command here
        gmx distance -f md_0_1.xtc -s md_0_1.tpr \
            -select "atomnr 903 plus atomnr 10001; atomnr 8466 plus atomnr 2438" \
            -oall distance.xvg

        # Return to parent directory
        cd ..
    fi
done
