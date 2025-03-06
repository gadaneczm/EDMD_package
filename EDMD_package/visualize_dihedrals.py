"""
This script prepares figures about the dihedral angle distribution for every residue in the structure ensemble.
"""

import json
from pathlib import Path
import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict

HIST_COLOR = "tab:blue"


def load_config(config_path: Path):
    """Load configuration from the provided JSON file."""

    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def tick_shifter(tick):
    """
    Change angle values from radian to degree and shift them between -180 and 180.

    :param tick: Angles in radian.
    :return: Angles in degree.
    """

    return f"{(180 / np.pi) * (tick if tick < np.pi else tick - 2 * np.pi):.0f}"


def progress_bar(percentage: float, length: int) -> str:

    n_of_hashtags = int(percentage * length)

    out = "["
    out += n_of_hashtags * "#"
    out += (length - n_of_hashtags) * " "
    out += "]"
    out += f" {percentage:.2%}"
    return out


def main(config_path: Path):

    print("\nvisualize_dihedrals.py is running:")

    # Load config data
    config = load_config(config_path)

    # Access global variables
    rosetta_results_folder = Path(config.get("ROSETTA_RESULTS_FOLDER"))

    # Read in the Rosetta dihedral angles from the pickle, written in save_dihedrals.py
    with open(rosetta_results_folder / f"angles_csr.pickle", "rb") as f:
        data: Dict[str, np.ndarray]
        data, _ = pickle.load(f)

    # Define the x values of the histogram
    hist_width = 2
    hist_x = np.arange(-180, 180 + hist_width, hist_width)
    plot_hist_x = (np.pi / 180) * (hist_x[:-1] + hist_x[1:]) / 2

    # Create the empty plots
    fig, ax = plt.subplots(1, 2, subplot_kw=dict(projection="polar"))
    plt.subplots_adjust(wspace=0.5)

    # Define keys to iterate the data
    keys = list({key[:-4] for key in data.keys()})
    keys.sort(key=lambda x: int(x.split("-")[0]))

    # Create a folder for the figures, if it doesn't exist yet
    if not os.path.exists(rosetta_results_folder / f"angle_figures_csr"):
        os.mkdir(rosetta_results_folder / f"angle_figures_csr")

    print("Saving figures...")

    # Iterate through the residues to create figures
    for counter, resi_name in enumerate(keys):

        phi_values = data[resi_name + " PHI"]
        psi_values = data[resi_name + " PSI"]

        # Define the y values of teh histograms for both dihedral angle types
        hist_phi, _ = np.histogram(phi_values, bins=hist_x)
        hist_psi, _ = np.histogram(psi_values, bins=hist_x)

        # Cleare the previous data from the figure
        [axis.cla() for axis in ax]

        # Set the titles
        ax[0].set_title(resi_name + " PHI")
        ax[1].set_title(resi_name + " PSI")

        # Plot the histograms
        ax[0].bar(plot_hist_x, hist_phi / np.max(hist_phi),
                  width=hist_width * np.pi / 180,
                  bottom=1,
                  color=HIST_COLOR)

        ax[1].bar(plot_hist_x, hist_psi / np.max(hist_psi),
                  width=hist_width * np.pi / 180,
                  bottom=1,
                  color=HIST_COLOR)

        # Set the x-labels
        x0_ticks = ax[0].get_xticks()
        ax[0].set_xticks(x0_ticks, labels=list(map(tick_shifter, x0_ticks)))

        x1_ticks = ax[1].get_xticks()
        ax[1].set_xticks(x1_ticks, labels=list(map(tick_shifter, x1_ticks)))

        # Set the y-labels
        ax[0].set_yticks(ax[0].get_yticks(), labels=["" for _ in ax[0].get_yticks()])
        ax[1].set_yticks(ax[1].get_yticks(), labels=["" for _ in ax[1].get_yticks()])

        # Save the figure
        fig.savefig(rosetta_results_folder / f"angle_figures_csr/{resi_name}.png", dpi=300)

        print("\r", end="")
        print(progress_bar((counter+1) / len(keys), 30), end=", ")
        print(f"Figure: {resi_name} done...", end="")
    print()
