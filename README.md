# Ensembe-Driven Molecular Dynamics (EDMD)
This module can be used to analyse the Phi and Psi dihedral (or torsion) angle distribution 
in a protein structural ensemble (e.g., from Chemical-Shift-Rosetta), define potential energy functions (PEFs), 
and replace the original dihedral energy terms in GROMACS for molecular dynamics (MD) simulations.

1. Set up the system for the MD simulation from the very best structure in the ensemble. 

2. Set your configuration in the "EDMD_config.json" file.

3. Run the main.py, where you can add the JSON file by the -c, or --config flag. This is a pipeline to run "save_dihedrals.py", "fit_dihedrals.py" and "create_tables.py". Optionally it can also call "visualize_dihedrals.py" and  "visualize_pef.py" if "VISUALIZE": True in the JSON.

## How to install?
Start by updating the pip version:
```bash
python3 -m pip install --upgrade pip
```

Install the EDMD_package:
```bash
python3 -m pip install EDMD_package
```

## EDMD_config.json file
- `ROSETTA_RESULTS_FOLDER: str` Path of the directory containing the ExtractedPDBs folder with the individual PDB files of the ensemble and a "name.scores.txt" containing model names and Rosetta-scores.

- `GMX_FOLDER: str` Path of the folder, where you want to run the MD simulation with the modified force field and where you have your TOP and GRO files.

- `RESI_IDX_SHIFT: int` Shift the residue numbering (if it was e.g. trimmed).

- `VISUALIZE: bool` Set True, if you want to run the visualize_dihedrals.py and visualize_pef.py scripts as well.

- `SCORE_SCALE: float or int` Set to scale the Rosetta-score for weighting during the PEF definition.

- `TEMPERATURE: float or int` Temperature of your simulation in Kelvin. Needed for the Boltzmann-inversion during the PEF definition.

- `GRO_FILENAME: str` Name of your GRO file, which is ready to run the simulation, so is solvated, etc.

- `TOP_FILENAME: str` Name of you processed TOP file (created e.g. by `gmx grompp -pp` flag in gromacs).

## How to use?
If you have set the EDMD_config.json , you can simply call:
```bash
python3 -m EDMD_package -c {path_to_JSON}
```

You can also call induvidual scripts:
```bash
pathon3 -m EDMD_package -c {path_to_JSON} -fn {name_of_script}
```

## Individual scripts
`save_dihedrals.py` The dihedral angles in your ensemble will be measured and saved to a pickle.

`fit_dihedrals.py` The probability density functions (PDF) will be defined for each backbone dihedral angle, 
according to the dihedral angle distributions using kernel density estimation. Finally, the PEFs will be created.

`create_tables.py` You need to have a ".gro" file and a ".top" file about your solvated system. 
By running this script you will get a ".new.top" file, which you should use as a topology file for your GROMACS MD simulation.

`visualize_dihedrals.py` Optionally, you can prepare figures about the dihedral angle distribution for every residue.

`visualize_pef.py` You can look at the angle distributions and the PEFs in case of each residue.

`format_scores_csrosetta` If you have already rescored the models (e.g. by running score_jd2 of Rosetta3) and 
you have a score.cs file, you can call this script to format the model names and scores so "save_dihedrals.py" will be able to read them 
and "fit_dihedrals.py" can use the Rosetta-score for weighting.


