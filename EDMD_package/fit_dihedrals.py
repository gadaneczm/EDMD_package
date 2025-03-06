"""
This script defines a probability density function (PDF) according to the
dihedral angle distributions of every residue using kernel density estimation (KDE).
After that the script defines the potential energy functions (PEFs) weighted by Rosetta scores.
"""

import json
from pathlib import Path
import pickle
from typing import Dict
import numpy as np

DELTA_ANGLE = 1.0
KERNEL_WIDTH_SCALE = 1


def load_config(config_path: Path):
    """Load configuration from the provided JSON file."""

    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def kde_kernel(distances: np.ndarray, kernel_width: float):
    """
    Kernel function definition.

    :param distances: The dih. angle distances (already corrected for periodicity).
    :param kernel_width: The width of the cosine kernel.
    :return: The PDF of a single kernel.
    """

    # Definition of the kernel function
    alpha = 1 + 1 / np.tan(kernel_width * np.pi / 360) ** 2

    kernel_values = np.cos(distances * np.pi / 360) ** alpha

    return kernel_values


def get_pdf(data: np.ndarray, scores: np.ndarray, score_scale: float) -> np.ndarray:
    """
    Definition of the Probability Density Function (PDF).

    :param data: Phi or Psi dihedral angles from the structure ensemble
    :param scores: Rosetta scores of the structure ensemble for weighting
    :param score_scale: Scaling factor for weighting
    :return: Probability function for Phi or Psi dihedral angles
    """

    ref_points = np.arange(-180, 180, DELTA_ANGLE)

    dmx = np.abs(ref_points[:, np.newaxis] - data[np.newaxis, :])
    mask = dmx > 180
    dmx[mask] = 360 - dmx[mask]

    # The contribution of the datapoints are weighted according to the Rosetta-scores
    weights = np.exp(-scores / score_scale)
    weights /= np.sum(weights)

    kernel_width = KERNEL_WIDTH_SCALE * 360 / len(data) ** (1 / 3)
    kernel_values = kde_kernel(dmx, kernel_width)
    kernel_values = np.sum(kernel_values * weights, axis=1)

    return kernel_values


def get_pef(data: np.ndarray, temperature: float) -> np.ndarray:
    """
    Definition of Potential Energy Function (PEF).

    :param data: PDF for the Phi or Psi dihedral angles
    :param temperature: Temperature (in Kelvin)
    :return: PEF
    """

    out = -1.9872036e-3 * temperature * np.log(data)
    out_mean = np.mean(out)
    out = out - out_mean

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

    print("\nfit_dihedrals.py is running:")

    # Load config data
    config = load_config(config_path)

    # Access global variables
    rosetta_results_folder: Path = Path(config.get("ROSETTA_RESULTS_FOLDER"))
    score_scale: float = config.get("score_scale")
    temperature: float = config.get("temperature")

    # Read in the dihedral data from a pickle
    with open(rosetta_results_folder / f"angles_csr.pickle", "rb") as f:
        data: Dict[str, np.ndarray]
        scores: np.ndarray
        data, scores = pickle.load(f)

    keys = list({key[:-4] for key in data.keys()})
    keys.sort(key=lambda x: int(x.split("-")[0]))

    # Define the x values of the histogram from -180 to 180 degrees in increments of one degree
    x_values = np.arange(-180, 180, DELTA_ANGLE)

    fit_data = (x_values, dict())

    print("Fitting PEF...")

    # Iterate through the residues to define PEFs
    resi_name: str
    for counter, resi_name in enumerate(keys):

        phi_key = resi_name + " PHI"
        psi_key = resi_name + " PSI"

        data_phi = data[phi_key]
        data_psi = data[psi_key]

        # Define PDFs using KDE
        pdf_phi = get_pdf(data_phi, scores, score_scale)
        pdf_psi = get_pdf(data_psi, scores, score_scale)

        # Define PEFs from PDFs
        pef_phi = get_pef(pdf_phi, temperature)
        pef_psi = get_pef(pdf_psi, temperature)

        # Calculate dPEFs
        dpef_phi = pef_phi[1:] - pef_phi[:-1]
        dpef_phi = np.append(dpef_phi, (dpef_phi[0] + dpef_phi[-1]) / 2) / DELTA_ANGLE

        dpef_psi = pef_psi[1:] - pef_psi[:-1]
        dpef_psi = np.append(dpef_psi, (dpef_psi[0] + dpef_psi[-1]) / 2) / DELTA_ANGLE

        fit_data[1][phi_key] = (pef_phi, dpef_phi)
        fit_data[1][psi_key] = (pef_psi, dpef_psi)

        print("\r", end="")
        print(progress_bar((counter+1) / len(keys), 30), end=", ")
        print(f"{resi_name} is done...", end="")
    print()

    # Save the PEF and dPEF data in a pickle
    pickle_name = f"pef_dpef_data_scoreScale{score_scale:.0f}_{temperature}K"
    with open(rosetta_results_folder / f"{pickle_name}.pickle", "wb") as f:
        pickle.dump(fit_data, f)
