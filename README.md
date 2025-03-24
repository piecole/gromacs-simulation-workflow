# GROMACS Simulation Workflow

This repository contains a SLURM-based workflow script for running GROMACS molecular dynamics simulations. The script automates the process of setting up and running MD simulations, including system preparation, equilibration, and production runs.

## Requirements

- GROMACS 2021.5 or later
- SLURM workload manager
- Required MDP files:
  - `ions.mdp`
  - `minim.mdp`
  - `nvt.mdp`
  - `npt.mdp`
  - `md.mdp` (or custom production MDP file)
  - `md_energy.mdp`
  - `index.ndx` (or custom index file)

## Usage

```bash
sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s] [-g gro_file]
```

### Arguments
- `input_file`: Base name of your input PDB file (without .pdb extension)
- `indicator`: Identifier for this specific run
- `-r`: Resume flag to skip to production run
- `-m`: Specify a custom production MDP file
- `-n`: Specify a custom index file
- `-s`: Stop the process just before pdb2gmx
- `-g gro_file`: Continue the process from after pdb2gmx with the specified .gro file

## Important Notes

This script is not universally applicable and requires customization for different systems:

1. **Energy Groups**: The `md.mdp` file needs to be customized with appropriate energy groups for your specific system. The current version may not include all necessary energy groups for your analysis.

2. **Force Field**: The script currently uses the CHARMM36-Jul2022 force field. If you need to use a different force field, you'll need to modify the `pdb2gmx` command in the script.

3. **System-Specific Parameters**: Various parameters in the MDP files (temperature, pressure, timestep, etc.) may need to be adjusted based on your specific system and requirements.

4. **Index Groups**: The script uses a default `index.ndx` file. You may need to create custom index groups for your specific analysis needs.

## Customization Required

Before using this script for a new system, you should:

1. Review and modify the `md.mdp` file to include appropriate energy groups
2. Check and adjust the force field settings if needed
3. Verify that the index groups in `index.ndx` match your analysis requirements
4. Adjust simulation parameters in the MDP files as needed
5. Consider whether you need to modify the box size or solvation settings

## Example

```bash
# Basic usage
sbatch run_full_gromacs_flow.sh protein_A run1

# Using custom MDP and index files
sbatch run_full_gromacs_flow.sh protein_A run1 -m custom_md.mdp -n custom_index.ndx

# Resuming from production run
sbatch run_full_gromacs_flow.sh protein_A run1 -r

# Stopping before pdb2gmx
sbatch run_full_gromacs_flow.sh protein_A run1 -s

# Continuing with a custom GRO file
sbatch run_full_gromacs_flow.sh protein_A run1 -g path/to/custom.gro
```

## Workflow Control

The script provides several options to control the workflow:

1. **Complete Workflow**: By default, runs the entire process from PDB to production.
2. **Pre-Processing Only**: Use the `-s` flag to prepare the directory and stop before pdb2gmx, useful when you need to do mandual gmx2pdb such as specifying disulphide bonds.
3. **Post-Processing with Custom Structure**: Use the `-g` flag with a .gro file to skip pdb2gmx and continue with your pre-processed structure.
4. **Resume Production**: Use the `-r` flag to skip straight to the production run phase. Most useful if compute time runs out before simulation finishes