#!/bin/bash -l
#SBATCH --partition=gpu
#SBATCH --time=12:00:00
#SBATCH --gpus=1
#SBATCH --cpus-per-gpu=16
#SBATCH --mem-per-gpu=32G
#SBATCH --job-name=gromacs
#SBATCH --output=gromacs-test-%j.out

module load gromacs/2021.5-gcc-11.4.0-cuda-11.8.0
export OMP_NUM_THREADS=64

if [ -e "md_0_1.cpt" ]; then
  gmx mdrun -cpi md_0_1.cpt -deffnm md_0_1
else
  gmx mdrun -deffnm md_0_1 -v -stepout 10000
fi
