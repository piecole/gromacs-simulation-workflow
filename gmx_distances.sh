#!/bin/bash

# Measure pairwise distances between all *custom* (user-defined) index groups
# for every trajectory found in the subdirectories.
#
# Standard GROMACS groups (System, Protein, Water, ions, ...) are ignored so
# that only the groups you added yourself with `gmx make_ndx` are measured.

index_file="index.ndx"
max_frames=1000

usage() {
    echo "Usage: gmx_distances.sh [-i index_file] [-m max_frames]"
    echo "  -i index_file    Name of the index (.ndx) file to use (default: index.ndx)"
    echo "  -m max_frames    Subsample each trajectory down to ~this many frames"
    echo "                   to speed up the distance calculation (default: 1000;"
    echo "                   use 0 to disable subsampling and use every frame)"
}

while getopts "i:m:h" opt; do
    case $opt in
        i)
            index_file=$OPTARG
            ;;
        m)
            max_frames=$OPTARG
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

# Make a results directory
mkdir -p results

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
    echo "Pairs:"
    for ((a = 0; a < ${#custom_groups[@]}; a++)); do
        for ((b = a + 1; b < ${#custom_groups[@]}; b++)); do
            pair="com of group \"${custom_groups[a]}\" plus com of group \"${custom_groups[b]}\""
            echo " $pair"
            if [ -z "$select_str" ]; then
                select_str="$pair"
            else
                select_str="$select_str; $pair"
            fi
        done
    done

    out_base=$(echo "${traj%.xtc}" | tr '/' '_')
    tpr="${traj%.xtc}.tpr"

    if [ ! -f "$tpr" ]; then
        echo " Skipping $traj: no matching $tpr"
        continue
    fi

    # Subsample on the fly: gmx distance processes every frame, so we tell it to
    # only read one frame every dt_new ps, targeting ~max_frames frames total.
    # (The trajectory-analysis tools support -b/-e/-dt but not -skip.)
    #
    # We read the timing from the .tpr, which is instant, rather than from
    # `gmx check -f`, which has to stream the whole trajectory to count frames.
    dt_opt=()
    if [ "$max_frames" -gt 0 ]; then
        # Pull dt (ps/step), nsteps, and the xtc write interval from the .tpr.
        read -r dt nsteps nstout < <(gmx dump -s "$tpr" 2>/dev/null | awk '
            /^[[:space:]]*dt[[:space:]]*=/                 {dt=$NF}
            /^[[:space:]]*nsteps[[:space:]]*=/             {ns=$NF}
            /(nstxout-compressed|nstxtcout)[[:space:]]*=/  {no=$NF}
            (dt!="" && ns!="" && no!="")                   {print dt, ns, no; exit}')

        if [ -z "$dt" ] || [ -z "$nstout" ] || [ "$nstout" -eq 0 ]; then
            echo " Could not read timing from $tpr; using every frame."
        else
            frame_dt=$(awk -v d="$dt" -v n="$nstout" 'BEGIN { printf "%g", d * n }')
            total_frames=$(( nsteps / nstout + 1 ))
            if [ "$total_frames" -le "$max_frames" ]; then
                echo " ~$total_frames frames (<= $max_frames); using every frame."
            else
                skip=$(( (total_frames + max_frames - 1) / max_frames ))
                dt_new=$(awk -v f="$frame_dt" -v s="$skip" 'BEGIN { printf "%g", f * s }')
                echo " ~$total_frames frames -> -dt $dt_new ps (~$(( total_frames / skip )) frames)"
                dt_opt=(-dt "$dt_new")
            fi
        fi
    fi

    echo "Running gmx distance..."
    gmx distance -f "$traj" -s "$tpr" -n "$index_file" "${dt_opt[@]}" \
        -select "$select_str" -oall results/distance_${out_base}_${index_file%.ndx}.xvg
done
