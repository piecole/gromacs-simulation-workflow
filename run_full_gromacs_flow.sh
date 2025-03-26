#!/bin/bash -l
#SBATCH --partition=gpu
#SBATCH --time=48:00:00
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=16
#SBATCH --mem-per-gpu=32G
#SBATCH --job-name=gromacs
#SBATCH --output=gromacs-workflow-%j.out
#SBATCH --mail-user=pierre.coleman@kcl.ac.uk
#SBATCH --mail-type=BEGIN,END,FAIL   # When to send the emails

# Help flag
if [ "$1" == "-h" ]; then
    echo "This script sets up and runs a full GROMACS simulation workflow."
    echo "Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s] [-g gro_file]"
    echo "Requirements:"
    echo "[input].pdb, ions.mdp, minim.mdp, nvt.mdp, npt.mdp, md.mdp (or other file), md_energy.mdp, index.ndx (or other file)."
    echo "Optional arguments:"
    echo "-m custom_mdp.mdp: Specify a custom production MDP file."
    echo "-n index.ndx: Specify an index file."
    echo "-r: Resume from production run."
    echo "-s: Stop the process just before pdb2gmx."
    echo "-g gro_file: Continue the process from after pdb2gmx with the specified .gro file."
    exit 0
fi

# Parse arguments
input_file=$1
indicator=$2
shift 2

# Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx]
production_mdp="md.mdp"
index_file="index.ndx"   # Default .ndx file
resume_flag=false
stop_before_pdb2gmx=false
gro_file=""
while getopts "rm:n:sg:" opt; do
    case $opt in
        r)
            resume_flag=true
            ;;
        m)
            production_mdp="$OPTARG"
            ;;
        n)
            index_file="$OPTARG"
            ;;
        s)
            stop_before_pdb2gmx=true
            ;;
        g)
            gro_file="$OPTARG"
            ;;
        \?)
            echo "Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s] [-g gro_file]"
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
echo "Run directory: $run_dir"


# Copy mdp files and make rundir  if rundir doesn't exist.
if [ ! -d "$run_dir" ]; then
    mkdir -p "$run_dir"
    missing_file=""
    cp ions.mdp "$run_dir" || missing_file+="ions.mdp "
    cp minim.mdp "$run_dir" || missing_file+="minim.mdp "
    cp nvt.mdp "$run_dir" || missing_file+="nvt.mdp "
    cp npt.mdp "$run_dir" || missing_file+="npt.mdp "
    cp "$production_mdp" "$run_dir" || missing_file+="$production_mdp "
    cp md_energy.mdp "$run_dir" || missing_file+="md_energy.mdp "

    if [ -n "$missing_file" ]; then
        echo "Error: Missing file(s): $missing_file"
        rm -r "$run_dir"
        exit 1
    fi
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

# Check for required arguments.
if [ $# -lt 2 ]; then
    echo "Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx]"
    echo "Optional arguments:"
    echo "-m custom_mdp.mdp: Specify a custom production MDP file."
    echo "-n custom_index.ndx: Specify a custom index file."
    echo "-r: Resume from production run."
    echo "-s: Stop the process just before pdb2gmx."
    echo "-g gro_file: Continue the process from after pdb2gmx with the specified .gro file."
    exit 1
fi

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

    # Check for structure file
    if [ ! -f "$input_file.pdb" ]; then
         echo "Error: Structure file '$input_file.pdb' not found."
         exit 1
    fi
    # Clean up: remove HOH lines from the pdb file and copy it to the run directory.
    grep -Ev "HOH|GOL|SO4|PEG|EDO|ACT|DMS|MPD|TRS|MES|HEPES" "$input_file.pdb" > "$run_dir/struc_clean.pdb"

    # Change to the run directory.
    cd "$run_dir" || exit 1

    # Stop before pdb2gmx if -s flag is set
    if [ "$stop_before_pdb2gmx" = true ]; then
        echo "Stopping before pdb2gmx as requested with -s flag."
        exit 0
    fi

    cp ../charmm36-jul2022.ff .

    # Skip pdb2gmx if a .gro file is provided with -g flag
    if [ -n "$gro_file" ]; then
        echo "Skipping pdb2gmx and using provided .gro file: $gro_file"
        if [ ! -f "$gro_file" ]; then
            echo "Error: Specified .gro file '$gro_file' not found."
            exit 1
        fi
    else
        # Run pdb2gmx.
        echo "1 0 1 0 1 0" | gmx pdb2gmx -f struc_clean.pdb \
                                       -o struc_processed.gro \
                                       -p topol.top \
                                       -water spce \
                                       -ff charmm36-jul2022 \
                                       -ter -ignh -merge all 2> pdb2gmx_error.log
        if [ $? -ne 0 ]; then
             echo "pdb2gmx failed. Check $run_dir/pdb2gmx_error.log for details."
             echo -e "\n---- Directory listing ----" >> "pdb2gmx_error.log"
             ls -halt >> "pdb2gmx_error.log"
             exit 1
        fi
    fi

    # Create a box.
    gmx editconf -f struc_processed.gro -o struc_newbox.gro -c -d 1.0 -bt cubic

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
    echo "16" | gmx genion -s ions.tpr \
                          -o struc_solv_ions.gro \
                          -p topol.top \
                          -pname NA -nname CL -neutral

    # Energy minimization.
    echo "Energy minimization..."
    gmx grompp -f minim.mdp -c struc_solv_ions.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em -s em.tpr -ntmpi 1
    echo "11 0" | gmx energy -f em.edr -o potential.xvg
    if [ $? -ne 0 ]; then
         echo "Energy minimization failed."
         exit 1
    fi

    # NVT equilibration.
    echo "Temperature equilibration..."
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -deffnm nvt -s nvt.tpr
    echo "17 0" | gmx energy -f nvt.edr -o temperature.xvg
    if [ $? -ne 0 ]; then
         echo "NVT equilibration failed."
         exit 1
    fi

    # NPT equilibration (with modified npt.mdp without positional restraints).
    echo "Pressure equilibration..."
    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt \
              -p topol.top -o npt.tpr
    gmx mdrun -deffnm npt -s npt.tpr
    echo "18 0" | gmx energy -f npt.edr -o pressure.xvg
    echo "24 0" | gmx energy -f npt.edr -o density.xvg
    if [ $? -ne 0 ]; then
         echo "NPT equilibration failed."
         exit 1
    fi
fi

#---------------------------------------------------
# Production run and energy calculation as one function.
run_production() {
    if [ -f "md_0_1.cpt" ]; then # If a checkpoint file exists, resume the production run.
         echo "Checkpoint file found. Resuming production MD from checkpoint."
         # Execute the production run with pullx and pullf files if they exist.
         # For some reason omitting this only gives an error when a run is resumed.
         if [ -f "md_0_1_pullx.xvg" ] && [ -f "md_0_1_pullf.xvg" ]; then
              gmx mdrun -deffnm md_0_1 \
                        -s md_0_1.tpr \
                        -cpi md_0_1.cpt \
                        -px md_0_1_pullx.xvg \
                        -pf md_0_1_pullf.xvg \
                        -v -stepout 10000 2> md_0_1_error.log
         else
              gmx mdrun -deffnm md_0_1 \
                        -s md_0_1.tpr \
                        -cpi md_0_1.cpt \
                        -v -stepout 10000 2> md_0_1_error.log
         fi
    else # If no checkpoint file exists, start the production run from the beginning.
         echo "No checkpoint file found. Starting production MD from the beginning."
         # Check if index file exists before proceeding
         if [ ! -f "../$index_file" ]; then
             echo "Index file '$index_file' not found."
             echo "Please create an index file using 'gmx make_ndx' before proceeding with the production run. Indicate index file with the -n flag."
             echo "You can then resume the workflow from here by using the -r flag."
             exit 1
         else
             cp "../$index_file" .
             echo "Index file '$index_file' found. Proceeding with production run."
         fi
         # Create the production run TPR.
         gmx grompp -f "$production_mdp" -c npt.gro -t npt.cpt -p topol.top -n "$index_file" -o md_0_1.tpr
         # Execute production run
         gmx mdrun -deffnm md_0_1 \
                   -s md_0_1.tpr \
                   -v -stepout 10000 2> md_0_1_error.log
    fi
    if [ $? -ne 0 ]; then # If the production run fails, exit with an error message.
         echo "Production MD run failed. Check $run_dir/md_0_1_error.log for details."
         exit 1
    fi

    # Energy calculation steps.
    gmx grompp -f md_energy.mdp -c md_0_1.gro \
              -p topol.top -o md_0_1_energy.tpr -n "$index_file"
    gmx mdrun -s md_0_1_energy.tpr \
              -rerun md_0_1.xtc \
              -deffnm md_0_1_energy \
              -ntmpi 1 -nb cpu 2> md_0_1_energy.log
    if [ $? -ne 0 ]; then
         echo "Energy calculation failed. Check $run_dir/md_0_1_energy.log for details."
         exit 1
    fi

    echo "Extracting interaction energy."
    echo "22 23 0" | gmx energy -f md_0_1_energy.edr -o interaction_energy.xvg
}

# Execute the production run (and energy calculation) block.
run_production
