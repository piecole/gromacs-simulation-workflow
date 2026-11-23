# GROMACS Simulation Workflow

This repository contains a SLURM-based workflow script for running GROMACS molecular dynamics simulations. The script automates the process of setting up and running MD simulations, including system preparation, equilibration, and production runs.

## Requirements

- GROMACS 2026.1 (`module load gromacs/2026.1`)
- SLURM workload manager
- GPU with CUDA support

## Usage

```bash
sbatch run_full_gromacs_flow.sh <input_file> <indicator> [-r] [-m custom_mdp.mdp] [-n custom_index.ndx] [-s] [-g gro_file] [-t] [-d 'y n ...'] [-c '0 0 ...']
```

### Arguments
- `input_file`: Base name of your input PDB file (without .pdb extension)
- `indicator`: Identifier for this specific run
- `-r`: Resume flag to skip to production run
- `-m`: Specify a custom production MDP file
- `-n`: Specify a index file
- `-s`: Stop the process just before pdb2gmx
- `-g gro_file`: Continue the process from after pdb2gmx with the specified .gro file
- `-t`: Force generation of new TPR file using existing checkpoint. Useful when you want to modify simulation parameters (e.g., extend simulation time) while keeping the same trajectory.
- `-d 'y n ...'`: Disulphide answer string. When given, pdb2gmx runs with `-ss` and these y/n answers are fed to its interactive S-S prompts (one per candidate cysteine pair, in the order pdb2gmx asks). Run pdb2gmx manually first to confirm the order/count. When omitted, `-ss` is not used and all proximal cysteines bond automatically (e.g. `-d 'y n y'`, or all `n` to suppress every disulphide).
- `-c '0 0 ...'`: Termini answer string. When given, pdb2gmx runs with `-ter` and these answers are fed to its terminus prompts, N-terminus then C-terminus for each chain in chain order. For this force field the menu is `0` = NH3+ / COO-, `1` = NH2 / COOH, and the C-terminal amide cap CT2 is `2`. A chain starting on GLY or PRO is offered different N-terminal options, so run pdb2gmx manually once to confirm the numbers. When omitted, `-ter` is not used and every chain gets the charged default (NH3+ and COO-).

## Workflow Steps

1. **System Setup**:
   - Creates a run directory with unique identifier
   - Copies all required MDP files
   - Cleans input PDB file by removing water and other small molecules

2. **Structure Processing**:
   - Runs pdb2gmx with CHARMM36-Jul2022 force field
   - Creates simulation box with 1.2 nm padding
   - Solvates system with TIP3P water model
   - Adds ions to neutralize the system and bring the salt to 0.15 M NaCl

3. **Equilibration**:
   - Energy minimization
   - NVT equilibration (temperature)
   - NPT equilibration (pressure)

4. **Production**:
   - Production MD run with GPU acceleration
   - Energy calculations
   - Interaction energy analysis between specified groups
   - Buried surface area (BSA) between the two groups via `gmx sasa` (BSA = SASA_group1 + SASA_group2 - SASA_complex)

## Important Notes

This script is not universally applicable and requires customization for different systems:

1. **Energy Groups**: The `md_energy.mdp` file needs to be customized with appropriate energy groups for your specific system. The script extracts interaction energy between two groups specified in this file.

2. **Force Field**: The script currently uses the CHARMM36-Jul2022 force field. If you need to use a different force field, you'll need to modify the `pdb2gmx` command in the script.

3. **System-Specific Parameters**: Various parameters in the MDP files (temperature, pressure, timestep, etc.) may need to be adjusted based on your specific system and requirements.

4. **Index Groups**: An `.ndx` index file is required for the production run. The script will check for the presence of the index file before starting the production run and stop if it's not found, allowing you to create the necessary index groups for your specific analysis needs. Use `-n` to specify index file in parent folder.

## Customization Required

Before using this script for a new system, you should:

1. Review and modify the `md_energy.mdp` file to include appropriate energy groups
2. Check and adjust the force field settings if needed
3. Adjust simulation parameters in the MDP files as needed
4. Consider whether you need to modify the box size or solvation settings

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
sbatch run_full_gromacs_flow.sh protein_A run1 -r -g path/to/custom.gro

# Forming specific disulphides (answers fed to pdb2gmx -ss prompts)
sbatch run_full_gromacs_flow.sh protein_A run1 -d 'y n y'

# Choosing termini (NH3+/COO- on chains 1 and 3, NH3+/CT2 on chains 2 and 4)
sbatch run_full_gromacs_flow.sh protein_A run1 -c '0 0 0 2 0 0 0 2'

# Extending simulation with modified parameters
sbatch run_full_gromacs_flow.sh protein_A run1 -r -t
sbatch run_full_gromacs_flow.sh protein_A run1 -r -t -m extended_md.mdp
```

## Workflow Control

The script provides several options to control the workflow:

1. **Complete Workflow**: By default, runs the entire process from PDB to production.
2. **Pre-Processing Only**: Use the `-s` flag to prepare the directory and stop before pdb2gmx, useful when you need to do manual gmx2pdb such as specifying disulphide bonds. Also may need to make a correct index file here.
3. **Post-Processing with Custom Structure**: Use the `-g` flag with a .gro file to skip pdb2gmx and continue with your pre-processed structure.
4. **Resume Production**: Use the `-r` flag to skip straight to the production run phase. Most useful if compute time runs out before simulation finishes.

## Output Files

The script generates several important output files:
- `interaction_energy.xvg`: Contains Coulomb and Lennard-Jones interaction energies between specified groups
- `interaction_surface_area.xvg`: Per-frame buried surface area (BSA) and interface area (BSA/2); intermediate `sasa_*.xvg` files are also written
- `md_0_1.xtc`: Trajectory file
- `md_0_1.gro`: Final structure
- Various energy and analysis files (temperature.xvg, pressure.xvg, density.xvg)
- Checkpoint files (.cpt) for resuming interrupted simulations
