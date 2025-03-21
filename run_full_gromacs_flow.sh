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
    echo "Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s 'CYS A 58 CYS B 158']"
    echo "Requirements:"
    echo "[input].pdb, ions.mdp, minim.mdp, nvt.mdp, npt.mdp, md.mdp (or other file), md_energy.mdp, index.ndx (or other file)."
    echo "Optional arguments:"
    echo "-m custom_mdp.mdp: Specify a custom production MDP file."
    echo "-n custom_index.ndx: Specify a custom index file."
    echo "-s 'CYS A 58 CYS B 158': Add a disulfide bond between specified cysteines (can be used multiple times for multiple bonds)."
    echo "resume=TRUE: Skip to production run."
    exit 0
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
    echo "Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s 'CYS A 58 CYS B 158']"
    echo "Optional arguments:"
    echo "-m custom_mdp.mdp: Specify a custom production MDP file."
    echo "-n custom_index.ndx: Specify a custom index file."
    echo "-s 'CYS A 58 CYS B 158': Add a disulfide bond between specified cysteines (can be used multiple times for multiple bonds)."
    echo "resume=TRUE: Skip setting up the folder and equilibriating system."
    exit 1
fi

# Parse arguments
input_file=$1
indicator=$2
shift 2

# Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s 'CYS A 58 CYS B 158']
production_mdp="md.mdp"
index_file="index.ndx"   # Default .ndx file
resume_flag=false
disulfide_bonds=()
while getopts "rm:n:s:" opt; do
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
            disulfide_bonds+=("$OPTARG")
            ;;
        \?)
            echo "Usage: sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s 'CYS A 58 CYS B 158']"
            exit 1
            ;;
    esac
done

# Define the run directory.
run_dir="run_${input_file}_${indicator}"

# Specify the force field folder as working directory.
export GMXLIB
GMXLIB=$(pwd)

# Copy mdp files and make folder
set_up_folder() {
    if [ ! -d "$run_dir" ]; then
        mkdir -p "$run_dir"
    fi
    missing_file=""
    cp ions.mdp "$run_dir" || missing_file+="ions.mdp "
    cp minim.mdp "$run_dir" || missing_file+="minim.mdp "
    cp nvt.mdp "$run_dir" || missing_file+="nvt.mdp "
    cp npt.mdp "$run_dir" || missing_file+="npt.mdp "
    cp "$production_mdp" "$run_dir" || missing_file+="$production_mdp "
    cp md_energy.mdp "$run_dir" || missing_file+="md_energy.mdp "
    cp "$index_file" "$run_dir" || missing_file+="$index_file "

    if [ -n "$missing_file" ]; then
        echo "Error: Missing file(s): $missing_file"
        exit 1
    fi
}

# If resume mode is requested, verify that the run directory and production TPR exist.
if [ "$resume_flag" = true ]; then
    echo "Resume flag detected. Skipping setup steps."
    if [ ! -d "$run_dir" ]; then
         echo "Error: Run directory '$run_dir' does not exist. Cannot resume."
         exit 1
    fi
    cd "$run_dir" || exit 1
    if [ ! -f "md_0_1.tpr" ]; then
         echo "Error: Production TPR file not found in '$run_dir'. Cannot resume."
         exit 1
    fi
else
    echo "No resume flag detected. Running full pre-production setup."
    
    # Make the run directory if it doesnt exist, and copy mdp files and index file into it.
    set_up_folder

    # Check for structure file
    if [ ! -f "$input_file.pdb" ]; then
         echo "Error: Structure file '$input_file.pdb' not found."
         exit 1
    fi
    # Clean up: remove HOH lines from the pdb file and copy it to the run directory.
    grep -v HOH "$input_file.pdb" > "$run_dir/struc_clean.pdb"
    
    # Add disulfide bond records to the PDB file if specified
    if [ ${#disulfide_bonds[@]} -gt 0 ]; then
        echo "Adding ${#disulfide_bonds[@]} disulfide bond(s) to PDB file"
        # Create a temporary file with SSBOND records
        tmp_file=$(mktemp)
        bond_count=1
        
          # Find line number of first ATOM or HETATM entry
          first_atom_line=$(grep -n "^ATOM\|^HETATM" "$run_dir/struc_clean.pdb" | head -1 | cut -d: -f1)
          
          if [ -n "$first_atom_line" ]; then
               # Copy header records
               head -n $((first_atom_line - 1)) "$run_dir/struc_clean.pdb" > "$tmp_file"
               
               # Add each disulfide bond record
               for bond in "${disulfide_bonds[@]}"; do
               # Parse the bond specification
               read -r res1 chain1 num1 res2 chain2 num2 <<< "$bond"
               echo "SSBOND   $bond_count $res1 $chain1   $num1    $res2 $chain2  $num2" >> "$tmp_file"
               bond_count=$((bond_count + 1))
               done
               
               # Add the ATOM/HETATM lines and the rest of the file
               tail -n +$first_atom_line "$run_dir/struc_clean.pdb" >> "$tmp_file"
          else
               # No ATOM/HETATM lines found, just add SSBOND at the top
               for bond in "${disulfide_bonds[@]}"; do
               # Parse the bond specification
               read -r res1 chain1 num1 res2 chain2 num2 <<< "$bond"
               echo "SSBOND   $bond_count $res1 $chain1   $num1    $res2 $chain2  $num2" >> "$tmp_file"
               bond_count=$((bond_count + 1))
               done
               
               # Add the entire original file
               cat "$run_dir/struc_clean.pdb" >> "$tmp_file"
          fi
        
        # Replace the original file with the new one containing SSBOND records
        mv "$tmp_file" "$run_dir/struc_clean.pdb"
        echo "Added disulfide bond records to PDB file"
    fi
    
    # Change to the run directory.
    cd "$run_dir" || exit 1

    # Run pdb2gmx.
    echo "1 0 1 0 1 0" | gmx pdb2gmx -f struc_clean.pdb \
                                   -o struc_processed.gro \
                                   -p topol.top \
                                   -water spce \
                                   -ff charmm36-jul2022 \
                                   -ter -ignh 2> pdb2gmx_error.log
    if [ $? -ne 0 ]; then
         echo "pdb2gmx failed. Check $run_dir/pdb2gmx_error.log for details."
         echo -e "\n---- Directory listing ----" >> "pdb2gmx_error.log"
         ls -halt >> "pdb2gmx_error.log"
         exit 1
    fi

    # Create a box.
    gmx editconf -f struc_processed.gro -o struc_newbox.gro -c -d 1.0 -bt cubic

    # Solvate.
    gmx solvate -cp struc_newbox.gro -cs spc216.gro -o struc_solv.gro -p topol.top 2> solvate_error.log
    if [ $? -ne 0 ]; then
         echo "Solvation failed. Check $run_dir/solvate_error.log for details."
         exit 1
    fi

    # Add ions.
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
    gmx grompp -f minim.mdp -c struc_solv_ions.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em -s em.tpr -ntmpi 1
    echo "11 0" | gmx energy -f em.edr -o potential.xvg
    if [ $? -ne 0 ]; then
         echo "Energy minimization failed."
         exit 1
    fi

    # NVT equilibration.
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -deffnm nvt -s nvt.tpr
    echo "17 0" | gmx energy -f nvt.edr -o temperature.xvg
    if [ $? -ne 0 ]; then
         echo "NVT equilibration failed."
         exit 1
    fi

    # NPT equilibration (with modified npt.mdp without positional restraints).
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
    gmx grompp -f "$production_mdp" -c md_0_1.gro \
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
