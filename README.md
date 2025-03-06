# Ensembe-Driven Molecular Dynamics (EDMD)
This module can be used to analyse the Phi and Psi dihedral (or torsion) angle distribution in a protein structural ensemble (e.g., from Chemical-Shift-Rosetta), define potential energy functions (PEFs), and replace the original dihedral energy terms in GROMACS for molecular dynamics (MD) simulations.

1. Set up the system for the MD simulation from the very best structure in the ensemble. 

2. Set your configuration in the "EDMD_config.json" file.

3. Run the main.py, where you can add the JSON file by the -c, or --config flag. This is a pipeline to run "save_dihedrals.py", "fit_dihedrals.py" and "create_tables.py". Optionally it can also call "visualize_dihedrals.py" and  "visualize_pef.py" if "VISUALIZE": True in the JSON.

save_dihedrals.py: The dihedral angles in your ensemble will be measured and saved to a pickle.

fit_dihedrals.py: The probability density functions (PDF) will be defined for each backbone dihedral angle, according to the dihedral angle distributions using kernel density estimation. Finally, the PEFs will be created.

create_tables.py: You need to have a ".gro" file and a ".top" file about your solvated system. By running this script you will get a ".new.top" file, which you should use as a topology file for your GROMACS MD simulation.

visualize_dihedrals.py: Optionally, you can prepare figures about the dihedral angle distribution for every residue.

visualize_pef.py: You can look at the angle distributions and the PEFs in case of each residue.

# EDMD_config.json file
ROSETTA_RESULTS_FOLDER: (string) Path of the directory containing the ExtractedPDBs folder with the individual PDB files of the ensemble and a "name.scores.txt" containing model names and Rosetta-scores.

GMX_FOLDER: (string) Path of the folder, where you want to run the MD simulation with the modified force field and where you have your TOP and GRO files.

RESI_IDX_SHIFT: (int) Shift the residue numbering (if it was e.g. trimmed).

VISUALIZE: (bool) Set True, if you want to run the visualize_dihedrals.py and visualize_pef.py scripts as well.

SCORE_SCALE: (float) Set to scale the Rosetta-score for weighting during the PEF definition.

TEMPERATURE: (float) Temperature of your simulation in Kelvin. Needed for the Boltzman-inversion during the PEF definition.

GRO_FILENAME: (string) Name of your GRO file.

TOP_FILENAME: (string) Name of you processed TOP file (created e.g. by gmx grompp -pp flag in gromacs).
