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
sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r(resume)] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s 'CYS A 58 CYS B 158']
```

### Arguments
- `input_file`: Base name of your input PDB file (without .pdb extension)
- `indicator`: Identifier for this specific run
- `-r`: Resume flag to skip to production run
- `-m`: Specify a custom production MDP file
- `-n`: Specify a custom index file
- `-s`: Add disulfide bond between specified cysteines (can be used multiple times)

## Important Notes

This script is not universally applicable and requires customization for different systems:

1. **Energy Groups**: The `md.mdp` file needs to be customized with appropriate energy groups for your specific system. The current version may not include all necessary energy groups for your analysis.

2. **Force Field**: The script currently uses the CHARMM36-Jul2022 force field. If you need to use a different force field, you'll need to modify the `pdb2gmx` command in the script.

3. **System-Specific Parameters**: Various parameters in the MDP files (temperature, pressure, timestep, etc.) may need to be adjusted based on your specific system and requirements.

4. **Index Groups**: The script uses a default `index.ndx` file. You may need to create custom index groups for your specific analysis needs.

5. **Disulfide Bonds**: While the script can handle disulfide bonds, you need to specify them explicitly using the `-s` option.

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

# With disulfide bonds
sbatch run_full_gromacs_flow.sh protein_A run1 -s 'CYS A 58 CYS B 158'

# Using custom MDP and index files
sbatch run_full_gromacs_flow.sh protein_A run1 -m custom_md.mdp -n custom_index.ndx

# Resuming from production run
sbatch run_full_gromacs_flow.sh protein_A run1 -r
``` 