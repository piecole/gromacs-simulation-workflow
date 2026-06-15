#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --job-name=gromacs_distances
#SBATCH --output=gromacs-distances-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL

module load gromacs/2021.5-gcc-11.4.0-cuda-11.8.0

# Measure pairwise distances between all *custom* (user-defined) index groups
# for every trajectory found in the subdirectories.
#
# Standard GROMACS groups (System, Protein, Water, ions, ...) are ignored so
# that only the groups you added yourself with `gmx make_ndx` are measured.

index_file="index.ndx"
max_frames=1000
# Group number to cluster + center on. Group 1 is Protein in GROMACS' default
# numbering (group 0 is always System), so the default just works for protein
# systems. Using a number avoids ambiguity when the .ndx has duplicate names.
center_group="1"

usage() {
    echo "Usage: gmx_distances.sh [-i index_file] [-m max_frames] [-c center_group]"
    echo "  -i index_file    Name of the index (.ndx) file to use (default: index.ndx)"
    echo "  -m max_frames    Subsample each trajectory down to ~this many frames"
    echo "                   to speed up the distance calculation (default: 1000;"
    echo "                   use 0 to disable subsampling and use every frame)"
    echo "  -c center_group  Index group NUMBER to cluster + center on before"
    echo "                   measuring distances, so the assembly stays whole"
    echo "                   across PBC (default: 1, which is Protein in GROMACS'"
    echo "                   default numbering)"
}

while getopts "i:m:c:h" opt; do
    case $opt in
        i)
            index_file=$OPTARG
            ;;
        m)
            max_frames=$OPTARG
            ;;
        c)
            center_group=$OPTARG
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

# Parallelism = the CPUs this job was actually allocated. Prefer the SLURM
# allocation (cgroup-limited) over nproc, which on a shared node reports every
# core on the machine, not just the ones we may safely use.
jobs=${SLURM_CPUS_PER_TASK:-${SLURM_JOB_CPUS_PER_NODE:-$(nproc 2>/dev/null || echo 1)}}

# Each gmx distance job is essentially serial; pin OpenMP to 1 thread so that
# J parallel jobs don't oversubscribe the cores by each spawning many threads.
export OMP_NUM_THREADS=1

# Standard GROMACS groups that should NOT be treated as custom selections.
default_groups=(
    "System" "Protein" "Protein-H" "C-alpha" "Backbone" "MainChain"
    "MainChain+Cb" "MainChain+H" "SideChain" "SideChain-H" "Prot-Masses"
    "non-Protein" "Other" "Water" "SOL" "non-Water" "Ion" "Water_and_ions"
    "DNA" "RNA" "NA" "CL" "K" "MG" "CA" "ZN" "NA+" "CL-"
)

# Process a single trajectory: find its custom groups, build the pairwise COM
# selection, and write all distances to one .xvg. Run as a background job; each
# invocation is its own subshell, so variables and output files never collide.
process_traj() {
    traj="$1"
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
        return
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
        return
    fi

    # Work out the subsampling step: we want ~max_frames frames total, so we
    # read one frame every dt_new ps. This -dt is handed to the trjconv pass
    # below, which writes the trimmed (and centered) trajectory.
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

    # Preprocess with trjconv before measuring distances. Raw trajectories let
    # the assembly drift across a box face, so one group's atoms (or a whole
    # chain) can end up imaged on the far side of the box and the COM distance
    # spikes. -pbc cluster reassembles everything in "$center_group" into one
    # contiguous cluster, then -center puts it in the middle of a compact box.
    # We also apply the subsampling here (-dt) so the trimmed, centered .xtc is
    # what gmx distance reads.
    centered="${traj%.xtc}_centered.xtc"

    echo "Centering on group $center_group with trjconv (-pbc cluster -center)..."
    # trjconv prompts, in order: clustering group, centering group, output group.
    printf '%s\n%s\n%s\n' "$center_group" "$center_group" "0" \
        | gmx trjconv -f "$traj" -s "$tpr" -n "$index_file" "${dt_opt[@]}" \
            -pbc cluster -center -ur compact -o "$centered"

    if [ ! -f "$centered" ]; then
        echo " trjconv did not produce $centered; skipping $traj."
        echo " (Does group $center_group exist in $index_file?)"
        return
    fi

    echo "Running gmx distance..."
    # The trajectory is now whole and centered, so no -dt is needed here.
    # whole_mol_com: COM from intact molecules; -pbc/-rmpbc: minimum-image distances.
    gmx distance -f "$centered" -s "$tpr" -n "$index_file" \
        -seltype whole_mol_com -selrpos whole_mol_com -pbc -rmpbc \
        -select "$select_str" -oall results/distance_${out_base}_${index_file%.ndx}.xvg
}

# Find all trajectories in subdirectories (actual files, not directories).
trajectories=( */*.xtc )

if [ ${#trajectories[@]} -eq 0 ]; then
    echo "No *.xtc files found in subdirectories."
    exit 1
fi

mkdir -p results

echo "Processing ${#trajectories[@]} trajectories, up to $jobs in parallel..."

# Launch one background job per trajectory, capping concurrency at $jobs.
# `wait -n` blocks until any single job finishes before we start the next.
running=0
for traj in "${trajectories[@]}"; do
    process_traj "$traj" &
    running=$((running + 1))
    if [ "$running" -ge "$jobs" ]; then
        wait -n
        running=$((running - 1))
    fi
done

# Wait for the remaining jobs to finish before exiting.
wait
echo "----------------------------------------"
echo "All trajectories done. Results in results/"
