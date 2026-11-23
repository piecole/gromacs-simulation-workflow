#!/bin/bash -l
#SBATCH --partition=gpu
#SBATCH --time=12:00:00
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=16
#SBATCH --mem-per-gpu=32G
#SBATCH --job-name=gromacs
#SBATCH --output=gromacs-test-%j.out

module load gromacs/2026.1

# Cluster provides the MPI build as gmx_mpi (not thread-MPI gmx).
GMX=gmx_mpi
if ! command -v "$GMX" >/dev/null 2>&1; then
  echo "Error: '$GMX' not found after module load. Check the module name or GMX setting."
  exit 1
fi

export OMP_NUM_THREADS=64

if [ -e "md_0_1.cpt" ]; then
  "$GMX" mdrun -cpi md_0_1.cpt -deffnm md_0_1
else
  "$GMX" mdrun -deffnm md_0_1 -v -stepout 10000
fi
