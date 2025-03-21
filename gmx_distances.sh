#!/bin/bash

folders="run_pepd_p53*restraint*"
module load gromacs
for folder in $folders; do
    echo "Processing $folder"
    cd $folder
    gmx distance -f md_0_1.xtc -s md_0_1.tpr -select "atomnr 903 plus atomnr 10001; atomnr 8466 plus atomnr 2438" -oall distance.xvg
    # Shorten the folder string to remove the run_pepd_p53 part
    folder_short=$(echo $folder | sed 's/run_pepd_p53_//')
    (head -n 30 distance.xvg; echo "..."; tail -n 10 distance.xvg) > distance_short_${folder_short}.xvg
    cd ..
done
