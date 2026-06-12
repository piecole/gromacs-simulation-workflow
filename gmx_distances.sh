#!/bin/bash

# Measure pairwise distances between all *custom* (user-defined) index groups
# for every trajectory found in the subdirectories.
#
# Standard GROMACS groups (System, Protein, Water, ions, ...) are ignored so
# that only the groups you added yourself with `gmx make_ndx` are measured.

index_file="index.ndx"

usage() {
    echo "Usage: gmx_distances.sh [-i index_file]"
    echo "  -i index_file    Name of the index (.ndx) file to use (default: index.ndx)"
}

while getopts "i:h" opt; do
    case $opt in
        i)
            index_file=$OPTARG
            ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            usage
            exit 1
            ;;
    esac
done

# Standard GROMACS groups that should NOT be treated as custom selections.
default_groups=(
    "System" "Protein" "Protein-H" "C-alpha" "Backbone" "MainChain"
    "MainChain+Cb" "MainChain+H" "SideChain" "SideChain-H" "Prot-Masses"
    "non-Protein" "Other" "Water" "SOL" "non-Water" "Ion" "Water_and_ions"
    "DNA" "RNA" "NA" "CL" "K" "MG" "CA" "ZN" "NA+" "CL-"
)

module load gromacs/2021.5-gcc-11.4.0-cuda-11.8.0

# Find all trajectories in subdirectories (actual files, not directories).
trajectories=( */*.xtc )

if [ ${#trajectories[@]} -eq 0 ]; then
    echo "No *.xtc files found in subdirectories."
    exit 1
fi

for traj in "${trajectories[@]}"; do
    echo "----------------------------------------"
    echo "Processing $traj"

    all_groups=($(grep -E '^\[.*\]' "$index_file" \
                | sed -E 's/^\[[[:space:]]*//; s/[[:space:]]*\][[:space:]]*$//'))

    custom_groups=()
    for group in "${all_groups[@]}"; do
        skip=false

        for remove in "${default_groups[@]}"; do
            if [ "$group" == "$remove" ]; then
                skip=true
                break
            fi
        done

        if [ "$skip" = false ]; then
            custom_groups+=("$group")
        fi
    done

    if [ ${#custom_groups[@]} -lt 2 ]; then
        echo " Skipping $traj: need at least 2 custom groups, found ${#custom_groups[@]} (${custom_groups[*]})"
        continue
    fi

    echo " Custom groups: ${custom_groups[*]}"

    # Build a selection for every unique pair of custom groups, measuring the
    # distance between their centres of mass.
    select_str=""
    for ((a = 0; a < ${#custom_groups[@]}; a++)); do
        for ((b = a + 1; b < ${#custom_groups[@]}; b++)); do
            pair="com of group \"${custom_groups[a]}\" plus com of group \"${custom_groups[b]}\""
            if [ -z "$select_str" ]; then
                select_str="$pair"
            else
                select_str="$select_str; $pair"
            fi
        done
    done

    
    out_base=$(echo "${traj%.xtc}" | tr '/' '_')
    gmx distance -f "$traj" -s "${traj%.xtc}.tpr" -n "$index_file" \
        -select "$select_str" -oall distance_${out_base}_${index_file%.ndx}.xvg
    
done
