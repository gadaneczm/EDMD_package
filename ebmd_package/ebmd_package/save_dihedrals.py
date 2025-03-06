"""
This script measures and saves the Phi and Psi dihedral angles of an ensemble of separate pdb files.
"""

import json
import numpy as np
import pickle
from pathlib import Path
from math import atan2, pi
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Chain import Chain


def load_config(config_path: Path):
    """Load configuration from the provided JSON file."""

    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def get_filenames(rosetta_folder: Path):
    """
    Collect the filenames of the pdb files.

    :return: A list of pdb filenames and a list of the belonging Rosetta scores.
    """

    # Read in the Rosetta scores and filenames for the models
    with open(Path(rosetta_folder / "ExtractedPDBs" / "name.scores.txt"), "r") as f:
        scores_data = f.read()

    scores_data = list(filter(lambda x: len(x) != 0 and
                              not any(x.startswith(word) for word in ["SEQUENCE", "description"]),
                              scores_data.split("\n")))

    # Sort the filenames according to the scores, to print the best and worst ones
    scores_data = list(map(lambda line: list(filter(lambda x: len(x) != 0, line.split(" "))), scores_data))
    scores_data = list(map(lambda x: (x[0], float(x[1])), scores_data))
    scores_data.sort(key=lambda x: x[1])

    print(f"First 3 file elements: {scores_data[:3]}")
    print(f"Last 3 file elements: {scores_data[-3:]}")

    file_names = list(map(lambda x: f"{x[0]}.pdb", scores_data))
    scores = list(map(lambda x: x[1], scores_data))

    return file_names, scores


def progress_bar(percentage: float, length: int) -> str:

    n_of_hashtags = int(percentage * length)

    out = "["
    out += n_of_hashtags * "#"
    out += (length - n_of_hashtags) * " "
    out += "]"
    out += f" {percentage:.2%}"
    return out


def get_dihedral(r1: np.ndarray, r2: np.ndarray, r3: np.ndarray, r4: np.ndarray) -> float:
    """
    Define Phi or Psi dihedral angles for the right residues from the coordinates of 4 backbone atoms.

    :param r1: Coordinates of C(i-1) for Phi, N(i) for Psi
    :param r2: Coordinates of N(i) for Phi, CA(i) for Psi
    :param r3: Coordinates of CA(i) for Phi, C(i) for Psi
    :param r4: Coordinates of C(i) for Phi, N(i+1) for Psi
    :return: The dihedral angle in degrees.
    """

    u1: np.ndarray = r2 - r1
    u2: np.ndarray = r3 - r2
    u3: np.ndarray = r4 - r3

    u12: np.ndarray = np.cross(u1, u2)
    u23: np.ndarray = np.cross(u2, u3)

    atan2_arg1 = np.linalg.norm(u2) * np.dot(u1, u23)
    atan2_arg2 = np.dot(u12, u23)

    return atan2(atan2_arg1, atan2_arg2)


def main(config_path: Path):

    print("\nsave_dihedrals.py is running:")

    # Load config data
    config = load_config(config_path)

    # Access global variables
    rosetta_results_folder = Path(config.get("ROSETTA_RESULTS_FOLDER"))
    resi_idx_shift: int = config.get("RESI_IDX_SHIFT")

    # Read filenames and Rosetta scores
    file_names, scores = get_filenames(rosetta_results_folder)

    scores = np.array(scores)
    angles_dict = dict()
    print("Measuring dihedrals...")

    # Iterate through the pdb files in the ROSETTA_FOLDER
    for file_idx, file_name in enumerate(file_names):

        # Read in the model (1st model, chain A) from a pdb file
        kras: Chain = PDBParser(QUIET=True).get_structure("kras",
                                                          str(rosetta_results_folder / "ExtractedPDBs" / file_name))[0]["A"]

        # Iterate through the residues and measuring the dihedral angles
        for resi_idx in range(1, len(kras) - 1):

            current_resi = kras.child_list[resi_idx]

            prev_c = kras.child_list[resi_idx - 1]["C"]

            curr_n = current_resi["N"]
            curr_ca = current_resi["CA"]
            curr_c = current_resi["C"]

            next_n = kras.child_list[resi_idx + 1]["N"]

            phi = get_dihedral(prev_c.coord, curr_n.coord, curr_ca.coord, curr_c.coord) * 180 / pi
            psi = get_dihedral(curr_n.coord, curr_ca.coord, curr_c.coord, next_n.coord) * 180 / pi

            # Collect the dihedral data to a dictionary
            pdb_resi_id = current_resi.full_id[3][1] + resi_idx_shift
            phi_key = f"{pdb_resi_id}-{current_resi.resname} PHI"
            psi_key = f"{pdb_resi_id}-{current_resi.resname} PSI"

            if phi_key in angles_dict:
                angles_dict[phi_key].append(phi)
            else:
                angles_dict[phi_key] = [phi, ]

            if psi_key in angles_dict:
                angles_dict[psi_key].append(psi)
            else:
                angles_dict[psi_key] = [psi, ]

        print("\r", end="")
        print(progress_bar((file_idx+1) / len(file_names), 30), end=", ")
        print(f"Filename {file_name} done...", end="")
    print()

    # Write the results to a pickle for later usage
    angles_dict = {key: np.array(value) for key, value in angles_dict.items()}

    with open(rosetta_results_folder / "angles_csr.pickle", "wb") as f:
        pickle.dump((angles_dict, scores), f)
