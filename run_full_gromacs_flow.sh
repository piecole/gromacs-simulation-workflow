#!/bin/bash -l
#SBATCH --partition=gpu
#SBATCH --time=48:00:00
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=16
#SBATCH --mem-per-gpu=32G
#SBATCH --job-name=gromacs
#SBATCH --output=gromacs-workflow-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL   # When to send the emails

# Help flag
if [ "$1" == "-h" ]; then
    echo "This script sets up and runs a full GROMACS simulation workflow."
    echo "Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s] [-g gro_file] [-t] [-d]"
    echo "Required arguments:"
    echo "  <input_file>    : Base name of the input PDB file (without .pdb extension)"
    echo "  <indicator>     : Unique identifier for this run"
    echo "Optional arguments:"
    echo "  -r             : Resume from production run"
    echo "  -m custom_mdp.mdp : Specify a custom production MDP file"
    echo "  -n custom_index.ndx : Specify a custom index file"
    echo "  -s             : Stop the process just before pdb2gmx"
    echo "  -g gro_file    : Continue the process from after pdb2gmx with the specified .gro file"
    echo "  -t             : Force generation of new TPR file using existing checkpoint"
    echo "  -d 'y n ...'   : Disulphide answer string. When given, pdb2gmx runs with -ss and"
    echo "                   these answers are fed to its interactive S-S prompts (one y/n per"
    echo "                   candidate cysteine pair, in the order pdb2gmx asks). Run pdb2gmx on"
    echo "                   the pdb manually first to confirm the order/count. When omitted, -ss"
    echo "                   is not used and all proximal cysteines bond automatically. Example:"
    echo "                   -d 'y n y' (use all 'n' to suppress every disulphide)."
    echo "Requirements:"
    echo "[input].pdb, ions.mdp, minim.mdp, nvt.mdp, npt.mdp, md.mdp (or other file), md_energy.mdp, index.ndx (or other file)"
    exit 0
fi

# Check for required arguments
if [ $# -lt 2 ]; then
    echo "Error: Missing required arguments"
    echo "Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s] [-g gro_file] [-t] [-d]"
    exit 1
fi

echo "------- git info -------"
echo "commit:"
git rev-parse HEAD
echo "branch:"
git branch --show-current
echo "status:"
git status
echo "------------------------"

# Parse arguments
input_file=$1
indicator=$2
shift 2

# Initialize variables
production_mdp="md.mdp"
index_file="index.ndx"
resume_flag=false
stop_before_pdb2gmx=false
gro_file=""
force_new_tpr=false
form_disulphide=""

# Parse optional arguments
while getopts "rm:n:sg:td:" opt; do
    case $opt in
        r)
            resume_flag=true
            ;;
        m)
            if [ -z "$OPTARG" ]; then
                echo "Error: -m requires a custom MDP file argument"
                exit 1
            fi
            production_mdp="$OPTARG"
            ;;
        n)
            if [ -z "$OPTARG" ]; then
                echo "Error: -n requires a custom index file argument"
                exit 1
            fi
            index_file="$OPTARG"
            ;;
        s)
            stop_before_pdb2gmx=true
            ;;
        g)
            if [ -z "$OPTARG" ]; then
                echo "Error: -g requires a gro file argument"
                exit 1
            fi
            gro_file="$OPTARG"
            ;;
        t)
            force_new_tpr=true
            ;;
        d)
            form_disulphide=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            echo "Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s] [-g gro_file] [-t] [-d]"
            exit 1
            ;;
    esac
done

# Define the run directory.
run_dir="run_${input_file}_${indicator}"

# Print all the arguments.
echo "Input file: $input_file"
echo "Indicator: $indicator"
echo "Resume flag: $resume_flag"
echo "Production MDP: $production_mdp"
echo "Index file: $index_file"
echo "Stop before pdb2gmx: $stop_before_pdb2gmx"
echo "Gro file: $gro_file"
echo "Force new TPR: $force_new_tpr"
echo "Form disulphides: $form_disulphide"
echo "Run directory: $run_dir"

# Make the run directory if it doesn't exist, then copy any mdp files that are
# not already present in it (so a pre-populated run directory is left intact).
mkdir -p "$run_dir"
missing_file=""
for mdp in ions.mdp minim.mdp nvt.mdp npt.mdp "$production_mdp" md_energy.mdp; do
    if [ ! -f "$run_dir/$mdp" ]; then
        cp "$mdp" "$run_dir" || missing_file+="$mdp "
    fi
done

if [ -n "$missing_file" ]; then
    echo "Error: Missing file(s): $missing_file"
    exit 1
fi

# Print info about GPU allocation.
echo "Allocated GPU(s): $CUDA_VISIBLE_DEVICES"
echo "GPU details:"
nvidia-smi --query-gpu=index,name,driver_version --format=csv,noheader

### In this section I'm tweaking between loading the slurm default gromacs
### or a more recent docker image using singularity.

module load gromacs/2021.5-gcc-11.4.0-cuda-11.8.0

# Switched to using docker image of more recent gromacs version, due to old version
# only seemingly being able to do one restraint at a time.
# mpirun used instead of srun to avoid compatibility issues with the mpi

# shopt -s expand_aliases
# alias gmx='mpirun singularity exec --nv ../gromacs_latest.sif gmx'

###

export OMP_NUM_THREADS=16

# Specify the force field folder as working directory.
export GMXLIB
GMXLIB=$(pwd)

# If resume mode is requested, verify that the run directory and production TPR exist.
if [ "$resume_flag" = true ]; then
    echo "Resume flag detected. Skipping setup steps."
    if [ ! -d "$run_dir" ]; then
         echo "Error: Run directory '$run_dir' does not exist. Cannot resume."
         exit 1
    fi
    cd "$run_dir" || exit 1
    # if [ ! -f "md_0_1.tpr" ]; then
    #     echo "Error: Production TPR file not found in '$run_dir'. Cannot resume."
    #     exit 1
    # fi
else
    echo "No resume flag detected. Running full pre-production setup."

    # Change to the run directory.
    cd "$run_dir" || exit 1

    # When not using a custom .gro (-g), prepare the cleaned input pdb that
    # pdb2gmx will read. With -g the input pdb is not needed at all, so we
    # neither require it nor overwrite struc_clean.pdb.
    if [ -z "$gro_file" ]; then
        if [ ! -f "../$input_file.pdb" ]; then
             echo "Error: Structure file '$input_file.pdb' not found."
             exit 1
        fi
        # Clean up: remove HOH and other small molecules from the input pdb.
        grep -Ev "HOH|GOL|SO4|PEG|EDO|ACT|DMS|MPD|TRS|MES|HEPES" "../$input_file.pdb" > struc_clean.pdb
    fi

    # Stop before pdb2gmx if -s flag is set
    if [ "$stop_before_pdb2gmx" = true ]; then
        echo "Stopping before pdb2gmx as requested with -s flag."
        exit 0
    fi

    cp ../charmm36-jul2022.ff . -r

    # Structure used for the box/solvation step. Defaults to the pdb2gmx output
    # but is overridden by the -g .gro file when provided.
    structure_gro="struc_processed.gro"

    # Skip pdb2gmx if a .gro file is provided with -g flag
    if [ -n "$gro_file" ]; then
        echo "Skipping pdb2gmx and using provided .gro file: $gro_file"
        if [ ! -f "$gro_file" ]; then
            echo "Error: Specified .gro file '$gro_file' not found in $run_dir."
            exit 1
        fi
        structure_gro="$gro_file"
    else
        # Build gmx termini string
        # First count the number of chains in the pdb file.
        num_chains=$(awk '/^ATOM/ {chains[substr($0,22,1)] = 1} END {print length(chains)}' struc_clean.pdb)
        # Then build the gmx termini string.
        termini_string=""
        for ((i=1; i<=num_chains; i++)); do
            termini_string+="1 0 "
        done
        echo "$num_chains chains found in pdb file. Running pdb2gmx with termini string: $termini_string."
        # Run pdb2gmx.

        if [ -n "$form_disulphide" ]; then
            echo "Forming disulphides with option string: $form_disulphide"
            disulphide_option="-ss"
            # -ss prompts (one y/n per candidate SG-SG pair) come before the
            # terminus prompts, so prepend the answers to the termini string.
            termini_string="${form_disulphide} ${termini_string}"
        else
            disulphide_option=""
        fi

        echo $termini_string | gmx pdb2gmx -f struc_clean.pdb \
                                       -o struc_processed.gro \
                                       -p topol.top \
                                       -water spce \
                                       -ff charmm36-jul2022 \
                                       -merge all \
                                       -ter -ignh $disulphide_option 2> pdb2gmx_error.log
        if [ $? -ne 0 ]; then
             echo "pdb2gmx failed. Check $run_dir/pdb2gmx_error.log for details."
             echo -e "\n---- Directory listing ----" >> "pdb2gmx_error.log"
             ls -halt >> "pdb2gmx_error.log"
             exit 1
        fi
    fi

    # Make solvation/ion addition idempotent. gmx solvate and gmx genion edit
    # the [ molecules ] section of topol.top in place, so re-running setup in a
    # directory that was already solvated would append duplicate SOL/ion lines
    # (topology atom count then no longer matches the coordinates). We snapshot
    # the pristine protein-only topology on the first run and restore it from
    # that snapshot on every subsequent run, before solvating.
    if [ ! -f topol.orig.top ]; then
        cp topol.top topol.orig.top
        echo "Saved pristine topology snapshot to topol.orig.top."
    else
        echo "Restoring pristine topology from topol.orig.top before solvation."
        cp topol.orig.top topol.top
    fi

    # Create a box.
    gmx editconf -f "$structure_gro" -o struc_newbox.gro -c -d 1.0 -bt cubic

    # Solvate.
    echo "Solvating..."
    gmx solvate -cp struc_newbox.gro -cs spc216.gro -o struc_solv.gro -p topol.top 2> solvate_error.log
    if [ $? -ne 0 ]; then
         echo "Solvation failed. Check $run_dir/solvate_error.log for details."
         exit 1
    fi

    # Add ions.
    echo "Adding ions..."
    gmx grompp -f ions.mdp -c struc_solv.gro -p topol.top -o ions.tpr 2> ions_grompp_error.log
    if [ $? -ne 0 ]; then
         echo "grompp for ions failed. Check $run_dir/ions_grompp_error.log for details."
         exit 1
    fi
    echo "SOL" | gmx genion -s ions.tpr \
                          -o struc_solv_ions.gro \
                          -p topol.top \
                          -pname NA -nname CL -neutral 2> genion_error.log
    if [ $? -ne 0 ]; then
         echo "genion failed. Check $run_dir/genion_error.log for details."
         exit 1
    fi

    # Energy minimization.
    echo "Energy minimization..."
    gmx grompp -f minim.mdp -c struc_solv_ions.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em -s em.tpr -ntmpi 1
    echo "Potential" | gmx energy -f em.edr -o potential.xvg
    if [ $? -ne 0 ]; then
         echo "Energy minimization failed."
         exit 1
    fi

    # NVT equilibration.
    echo "Temperature equilibration..."
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -deffnm nvt -s nvt.tpr -ntmpi 1 -ntomp 16
    echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg
    if [ $? -ne 0 ]; then
         echo "NVT equilibration failed."
         exit 1
    fi

    # NPT equilibration (with modified npt.mdp without positional restraints).
    echo "Pressure equilibration..."
    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt \
              -p topol.top -o npt.tpr
    gmx mdrun -deffnm npt -s npt.tpr -ntmpi 1 -ntomp 16
    echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg
    echo "Density" | gmx energy -f npt.edr -o density.xvg
    if [ $? -ne 0 ]; then
         echo "NPT equilibration failed."
         exit 1
    fi
fi

#---------------------------------------------------
# Production run and energy calculation as one function.
run_production() {
    # If final .gro file exists, skip the production run.
    if [ -f "md_0_1.gro" ] && [ "$force_new_tpr" = false ]; then
         echo "Final .gro file found. Skipping production run."
    elif [ -f "md_0_1.cpt" ] && [ "$force_new_tpr" = false ]; then # If a checkpoint file exists, resume the production run.
         echo "Checkpoint file found."
         if [ "$force_new_tpr" = true ]; then
             echo "Force flag detected. Generating new TPR file from checkpoint."
             # grep production mdp for gen-vel, if gen-vel is yes then
             # return an error. This must be false if you want to continue
             # from a checkpoint properly including the velocities.
             gen_vel=$(sed -n 's/^gen_vel[ \t]*=[ \t]*\([^ \t;]*\).*/\1/p' "$production_mdp")

             if [ "${gen_vel,,}" = "yes" ]; then
                 echo "Error: gen_vel = yes in $production_mdp. Please set gen_vel = no to continue from a checkpoint."
                 exit 1
             fi

             gmx grompp -f "$production_mdp" -c npt.gro -t npt.cpt -p topol.top -n "$index_file" -o md_0_1.tpr
         fi
         echo "Resuming production MD from checkpoint. See readout at $run_dir/md_0_1_error.log"
         # Execute the production run with pullx and pullf files if they exist.
         # For some reason omitting this only gives an error when a run is resumed.
         if [ -f "md_0_1_pullx.xvg" ] && [ -f "md_0_1_pullf.xvg" ]; then
              gmx mdrun -deffnm md_0_1 \
                        -s md_0_1.tpr \
                        -cpi md_0_1.cpt \
                        -px md_0_1_pullx.xvg \
                        -pf md_0_1_pullf.xvg \
                        -ntmpi 1 -ntomp 16 \
                        -nb gpu \
                        -pme gpu \
                        -v -stepout 10000 2> md_0_1_error.log
         else
              gmx mdrun -deffnm md_0_1 \
                        -s md_0_1.tpr \
                        -cpi md_0_1.cpt \
                        -ntmpi 1 -ntomp 16 \
                        -nb gpu \
                        -pme gpu \
                        -v -stepout 10000 2> md_0_1_error.log
         fi
    else # If no checkpoint file exists, start the production run from the beginning.
         echo "No checkpoint file found. Starting production MD from the beginning."
         # Check if index file exists before proceeding
         if [ ! -f "../$index_file" ]; then
             echo "Index file '$index_file' not found."
             echo "Please create an index file using 'gmx make_ndx' before proceeding with the production run. Indicate index file with the -n flag."
             echo "You can then resume the workflow from here by using the -r flag."
             echo "----------------------------------------"
             echo "Load:"
             echo "> module load gromacs/2021.5-gcc-11.4.0-cuda-11.8.0"
             echo "To identify chains:"
             echo "> grep -v SOL $run_dir/npt.gro | grep CA"
             echo "To create index file:"
             echo "> gmx make_ndx -f $run_dir/npt.gro -o index.ndx"
             echo "To resume the workflow:"
             echo "> sbatch run_full_gromacs_flow.sh $input_file $indicator -r"
             echo "----------------------------------------"
             exit 1
         else
             cp "../$index_file" .
             echo "Index file '$index_file' found. Proceeding with production run."
         fi
         # Create the production run TPR.
         gmx grompp -f "$production_mdp" -c npt.gro -t npt.cpt -p topol.top -n "$index_file" -o md_0_1.tpr
         # Execute production run
         echo "Executing production run. See readout at $run_dir/md_0_1_error.log"
         gmx mdrun -deffnm md_0_1 \
                   -s md_0_1.tpr \
                   -ntmpi 1 -ntomp 16 \
                   -nb gpu \
                   -pme gpu \
                   -v -stepout 10000 2> md_0_1_error.log
    fi
    if [ $? -ne 0 ]; then # If the production run fails, exit with an error message.
         echo "Production MD run failed. Check $run_dir/md_0_1_error.log for details."
         exit 1
    fi

    # Energy calculation steps.
    gmx grompp -f md_energy.mdp -c md_0_1.gro \
              -p topol.top -o md_0_1_energy.tpr -n "$index_file"
    echo "Executing energy calculation. See readout at $run_dir/md_0_1_energy.log"
    gmx mdrun -s md_0_1_energy.tpr \
              -rerun md_0_1.xtc \
              -deffnm md_0_1_energy \
              -ntmpi 1 -nb cpu 2> md_0_1_energy.log
    if [ $? -ne 0 ]; then
         echo "Energy calculation failed. Check $run_dir/md_0_1_energy.log for details."
         exit 1
    fi

    echo "Getting energy groups from md_energy.mdp."
    groups=$(grep "energygrps" md_energy.mdp | awk -F'=' '{print $2}' | xargs)
    group1=$(echo "$groups" | awk '{print $1}')
    group2=$(echo "$groups" | awk '{print $2}')
    # Return the energy groups. Handle errors if not found.
    if [ -z "$group1" ] || [ -z "$group2" ]; then
         echo "Error: Energy groups not found in md_energy.mdp. Define as energygrps = protein_1 protein_2 in the mdp file. With groups referenced in the index file."
         exit 1
    fi
    echo "Energy groups: $group1 and $group2"

    echo "Extracting interaction energy."
    echo -e "Coul-SR:${group1}-${group2}\nLJ-SR:${group1}-${group2}" | \
    gmx energy -f md_0_1_energy.edr -o interaction_energy.xvg

    echo "Calculating interaction surface area (buried surface area / BSA)."
    # SASA is computed for each binding partner in isolation and for the
    # complex. gmx sasa ignores atoms outside the -surface selection, so the
    # two single-partner runs report free-state SASA. The buried (interface)
    # surface area then follows from BSA = SASA_1 + SASA_2 - SASA_complex.
    gmx sasa -f md_0_1.xtc -s md_0_1_energy.tpr -n "$index_file" \
             -surface "group \"$group1\"" \
             -o "sasa_${group1}.xvg" 2> "sasa_${group1}.log"
    gmx sasa -f md_0_1.xtc -s md_0_1_energy.tpr -n "$index_file" \
             -surface "group \"$group2\"" \
             -o "sasa_${group2}.xvg" 2> "sasa_${group2}.log"
    gmx sasa -f md_0_1.xtc -s md_0_1_energy.tpr -n "$index_file" \
             -surface "group \"$group1\" or group \"$group2\"" \
             -o "sasa_complex.xvg" 2> "sasa_complex.log"

    if [ ! -f "sasa_${group1}.xvg" ] || [ ! -f "sasa_${group2}.xvg" ] || [ ! -f "sasa_complex.xvg" ]; then
         echo "Error: SASA calculation failed. Check sasa_*.log for details."
         exit 1
    fi

    # Combine the per-frame SASA values into the buried surface area.
    grep -v '^[#@]' "sasa_${group1}.xvg" | awk '{print $1, $2}' > sasa_1.tmp
    grep -v '^[#@]' "sasa_${group2}.xvg" | awk '{print $2}'      > sasa_2.tmp
    grep -v '^[#@]' "sasa_complex.xvg"   | awk '{print $2}'      > sasa_c.tmp

    {
        echo "@    title \"Interaction (buried) surface area\""
        echo "@    xaxis  label \"Time (ps)\""
        echo "@    yaxis  label \"Area (nm\\S2\\N)\""
        echo "@ s0 legend \"Buried surface area (BSA)\""
        echo "@ s1 legend \"Interface area (BSA/2)\""
        paste sasa_1.tmp sasa_2.tmp sasa_c.tmp | awk '{
            bsa = $2 + $3 - $4;
            printf "%s %.4f %.4f\n", $1, bsa, bsa/2.0
        }'
    } > interaction_surface_area.xvg

    # Report the mean buried surface area over the trajectory.
    mean_bsa=$(paste sasa_1.tmp sasa_2.tmp sasa_c.tmp | awk '{
        bsa = $2 + $3 - $4; sum += bsa; n++
    } END { if (n > 0) printf "%.4f", sum/n }')
    echo "Mean buried surface area: ${mean_bsa} nm^2 (mean interface area: $(awk "BEGIN{printf \"%.4f\", ${mean_bsa}/2.0}") nm^2)"

    rm -f sasa_1.tmp sasa_2.tmp sasa_c.tmp
    }

# Execute the production run (and energy calculation) block.
run_production

echo "Production run, energy, and surface area calculation complete."