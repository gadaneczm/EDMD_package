"""
Creating tables with tabulated potentials for every Phi & Psi dihedral angle in each residue.
This should be run after preparing the system for the simulation (having a GRO and TOP file)
"""

import json
import os
import shutil
import pickle
import numpy as np
from typing import List, Dict, Tuple, Iterable
from pathlib import Path
import sys


def load_config(config_path: Path):
    """Load configuration from the provided JSON file."""

    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def create_id_key(ids: Iterable[int]):
    return frozenset({frozenset(ids[:3]), frozenset(ids[1:])})


def parse_gro(gro_path: Path, keys: List[str]) -> Dict[str, Tuple[int, int, int, int]]:
    """
    Parse the GRO file to collect dihedral IDs and atom indexes for each Phi/Psi in the sequence.

    :param gro_path: Path of the GRO file.
    :param keys: List of amino acids (e.g. 1-MET).
    :return: Dictionary containing the atom indexes for Phi and Psi dihedral angles in each residue.
    """

    with open(gro_path, "r") as f:
        gro_data = f.read()

    gro_data = gro_data.split("\n")[2:-2]

    gro_data = {
        (int(line[:5]), line[10:15].strip()): int(line[15:20])
        for line in gro_data
    }

    out = dict()
    for resi_name in keys:

        resi_idx = int(resi_name.split("-")[0])

        phi_ids = (
            gro_data[(resi_idx - 1, "C")],
            gro_data[(resi_idx, "N")],
            gro_data[(resi_idx, "CA")],
            gro_data[(resi_idx, "C")],
        )

        psi_ids = (
            gro_data[(resi_idx, "N")],
            gro_data[(resi_idx, "CA")],
            gro_data[(resi_idx, "C")],
            gro_data[(resi_idx + 1, "N")],
        )

        out[resi_name + " PHI"] = phi_ids
        out[resi_name + " PSI"] = psi_ids

    return out


def parse_top(top_path: Path) -> (List[str], Dict[frozenset[frozenset[int, int, int], frozenset[int, int, int]], int]):
    """
    Read the TOP file and collect the row indexes for the appropriate dihedral angles.

    :param top_path: Path of the TOP file.
    :return: List containing the lines of the TOP file, and a dictionary with the atom index pairs
             and row indexes for each dihedral angle.
    """

    with open(top_path, "r") as f:
        top_data = f.read()

    top_data = top_data.split("\n")

    # Get the index pairs for the boundaries of the dihedral sections.
    dih_header_idxs = [None, ]
    for idx, line in enumerate(top_data):

        if line.startswith("[") and type(dih_header_idxs[-1]) is int:
            dih_header_idxs[-1] = (dih_header_idxs[-1], idx)

        if line.startswith("[ dihedrals ]"):
            dih_header_idxs.append(idx)

    if type(dih_header_idxs[-1]) is int:
        dih_header_idxs[-1] = (dih_header_idxs[-1], len(top_data))

    dih_header_idxs = dih_header_idxs[1:]  # remove dummy None

    # Create the atom id quartets and assign the top file row indices to them.
    dih_ids_to_row_idxs = dict()
    for section_start, section_end in dih_header_idxs:

        section = top_data[section_start + 1:section_end]

        for idx, line in enumerate(section):

            line = line.split(";")[0]
            line = list(filter(lambda x: len(x) != 0, line.split(" ")))

            if len(line) == 0:
                continue

            key = list(map(int, line[:4]))
            # key = tuple(key) if key[0] < key[-1] else tuple(key[::-1])
            key = create_id_key(key)

            dih_ids_to_row_idxs[key] = section_start + 1 + idx

    return top_data, dih_ids_to_row_idxs


def get_xvg(ref_points: np.ndarray, pef_data: np.ndarray, dpef_data: np.ndarray) -> str:
    """
    Create a table of angles, and belonging PEF and dPEF values, defined by fit_dihedrals.py.

    :param ref_points: List of angles.
    :param pef_data: List of PEF values.
    :param dpef_data: List of dPEF values.
    :return: Table containing the angles, PEFs, dPEFs.
    """

    out = ""
    for x, y, dy in zip(ref_points, pef_data, dpef_data, strict=True):
        out += f"{x:.0f}\t{y:.10f}\t{-dy:.10f}\n"

    out += f"180\t{pef_data[0]:.10f}\t{-dpef_data[0]:.10f}\n"

    return out


def progress_bar(percentage: float, length: int) -> str:

    n_of_hashtags = int(percentage * length)

    out = "["
    out += n_of_hashtags * "#"
    out += (length - n_of_hashtags) * " "
    out += "]"
    out += f" {percentage:.2%}"
    return out


def main(config_path: Path):

    print("\ncreate_tables.py is running:")

    # Load config data
    config = load_config(config_path)

    # Access global variables
    rosetta_results_folder = Path(config.get("ROSETTA_RESULTS_FOLDER"))
    gmx_folder = Path(config.get("GMX_FOLDER"))
    score_scale: float = config.get("SCORE_SCALE")
    top_filename: str = config.get("TOP_FILENAME")
    gro_filename: str = config.get("GRO_FILENAME")
    temperature: float = config.get("TEMPERATURE")

    # Import PEFs from pickle, written by fit_dihedrals.py
    x_values: np.ndarray
    pef_dpef_data: Dict[str, Tuple[np.ndarray, np.ndarray]]

    pickle_name = f"pef_dpef_data_scoreScale{score_scale:.0f}_{temperature}K"
    with open(rosetta_results_folder / f"{pickle_name}.pickle", "rb") as f:
        x_values, pef_dpef_data = pickle.load(f)

    # Collect the amino acid sequence
    keys = list({key[:-4] for key in pef_dpef_data.keys()})
    keys.sort(key=lambda x: int(x.split("-")[0]))

    # Verify if the GRO and TOP files exist
    if not os.path.exists(gmx_folder / gro_filename):
        sys.exit(f"There is no GRO file in the GMX_FOLDER")
    if not os.path.exists(gmx_folder / top_filename):
        sys.exit(f"There is no TOP file in the GMX_FOLDER")

    # Parse the GRO file to collect dihedral IDs and atom indexes for each Phi/Psi in the sequence
    resi_to_ids = parse_gro(gmx_folder / gro_filename, keys)

    # Read the TOP file and collect the atom index pairs and row indexes for each dihedral angles
    top_data, dih_ids_to_row_idxs = parse_top(gmx_folder / top_filename)

    # Create the necessary dictionaries
    if not os.path.exists(gmx_folder / "for_gmx"):
        os.mkdir(gmx_folder / "for_gmx")

    out_folder = f"for_gmx/tables_of_potentials"
    if not os.path.exists(gmx_folder / out_folder):
        os.mkdir(gmx_folder / out_folder)

    print("Writing tables...")

    # Iterate through the Phi/Psi angles
    for table_idx, angle_name in enumerate(resi_to_ids):

        atom_ids_tuple = resi_to_ids[angle_name]

        atom_ids_frset = create_id_key(atom_ids_tuple)
        top_row_idx = dih_ids_to_row_idxs[atom_ids_frset]

        force_constant = 1

        # Modify the lines of the TOP file describing the Phi/Psi dihedral potentials
        new_row = " ".join(map(lambda x: f"{x:5d}", atom_ids_tuple))  # atom ids
        new_row += "     8"  # tabulated dihedral function type
        new_row += f" {table_idx:5d}"  # table index
        new_row += f" {force_constant}"
        new_row += f" ;{top_data[top_row_idx]}"

        top_data[top_row_idx] = new_row

        pef_dpef_table = get_xvg(x_values, *pef_dpef_data[angle_name])

        with open(gmx_folder / f"{out_folder}/table_d{table_idx}.xvg", "w") as f:
            f.write(pef_dpef_table)

        print("\r", end="")
        print(progress_bar((table_idx+1) / len(resi_to_ids), 30), end=", ")
        print(f"Filename {angle_name} done...", end="")
    print()

    # Create a NEW.TOP file
    top_data = "\n".join(top_data)

    new_top_filename = f"{top_filename[:-9]}.new.top"
    with open(gmx_folder / f"for_gmx/{new_top_filename}", "w") as f:
        f.write(top_data)

    # Copy the necessary files to a "for_gmx" folder
    shutil.copy(gmx_folder / gro_filename, gmx_folder / "for_gmx" / gro_filename)
    shutil.copy(gmx_folder / "../../../gromacs_processing/edmd/EDMD_config.py",
                gmx_folder / f"{out_folder}/EDMD_config.py")
