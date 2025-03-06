# NEW_NAME
This modul can be used to analyse the Phi and Psi dihedral (or torsion) angle distribution in a protein structural ensemble (e.g. from Chemical-Shift-Rosetta), and to define potential energy functions (PEFs), to exchange them for the original torsional energy terms in GROMACS for a molecular dynamics (MD) simulation.


1. Set up the system for the simulation from the very best structure in the ensemble. 

2. Set your configuration in the "config.json" file.

3. Run the main.py, where you can add the JSON file by the -c, or --config flag. This is a pipeline to run "save_dihedrals.py", "fit_dihedrals.py" and "create_tables.py". Optionally it can also call "visualize_dihedrals.py" and  "visualize_pef.py" if "VISUALIZE": True in the JSON.

save_dihedrals.py: The dihedral angles in your ensemble will be measured and saved to a pickle.

fit_dihedrals.py: The probability density functions (PDF) will be defined for each backbone dihedral angle, according to the dihedral angle distributions using kernel density estimation. Finally, the PEFs will be created.

create_tables.py: You need to have a ".gro" file and a ".top" file about your solvated system. By running the "create_tables.py" you will get a ".new.top" file, which you should use as a topology file for your GROMACS MD simulation.

visualize_dihedrals.py: Optionally, you can prepare figures about the dihedral angle distribution for every residue.

visualize_pef.py: You can look at the angle distributions and the PEFs in case of each residue.

