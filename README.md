# GROMACS Simulation Workflow

This repository contains scripts and configuration files for running molecular dynamics simulations using GROMACS.

## Files

- `run_full_gromacs_flow.sh`: Main script for running the complete GROMACS workflow
- `run_sim.sh`: Script for running individual simulations
- `gmx_distances.sh`: Script for calculating distances
- `get_distances.sh`: Script for extracting distance data
- Various `.mdp` files for different simulation stages:
  - `ions.mdp`: Ion addition parameters
  - `minim.mdp`: Energy minimization parameters
  - `nvt.mdp`: NVT equilibration parameters
  - `npt.mdp`: NPT equilibration parameters
  - `md.mdp`: Production MD parameters
  - `md_energy.mdp`: Energy calculation parameters
  - `md_restraint.mdp`: Restraint parameters

## Usage

To run a complete simulation workflow:

```bash
sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s 'CYS A 58 CYS B 158']
```

### Arguments:
- `input_file`: Name of the input PDB file (without .pdb extension)
- `indicator`: Identifier for the run
- `-r`: Resume flag to skip setup steps
- `-m`: Specify custom production MDP file
- `-n`: Specify custom index file
- `-s`: Add disulfide bond between specified cysteines

### Requirements:
- GROMACS 2021.5 or later
- SLURM job scheduler
- GPU access
- Required input files:
  - `[input].pdb`
  - `ions.mdp`
  - `minim.mdp`
  - `nvt.mdp`
  - `npt.mdp`
  - `md.mdp`
  - `md_energy.mdp`
  - `index.ndx`

## License

[Add your chosen license here] 