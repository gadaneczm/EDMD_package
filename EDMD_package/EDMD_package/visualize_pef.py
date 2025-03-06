"""
This script prepares figures about the dihedral angle distribution and
the belonging potential energy functions (PEFs)
for every residue in the structure ensemble.
"""

import json
from pathlib import Path
import os
import pickle
import numpy as np
from typing import Dict, Tuple
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

HIST_COLOR_CSR = "tab:blue"
TAG = ""
LEGEND = True


def load_config(config_path: Path):
    """Load configuration from the provided JSON file."""

    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def progress_bar(percentage: float, length: int) -> str:

    n_of_hashtags = int(percentage * length)

    out = "["
    out += n_of_hashtags * "#"
    out += (length - n_of_hashtags) * " "
    out += "]"
    out += f" {percentage:.2%}"
    return out


def main(config_path: Path):

    print("\nvizualize_pef.py is running:")

    # Load config data
    config = load_config(config_path)

    # Access global variables
    rosetta_results_folder = Path(config.get("ROSETTA_RESULTS_FOLDER"))
    score_scale: float = config.get("SCORE_SCALE")
    temperature: float = config.get("TEMPERATURE")

    # Read in the Rosetta dihedral angles from the pickle, written in save_dihedrals.py
    with open(rosetta_results_folder / f"angles_csr.pickle", "rb") as f:
        angle_data, _ = pickle.load(f)

    # Read in PEF data from a pickle, written in fit_dihedrals.py
    pickle_name = f"pef_dpef_data_scoreScale{score_scale:.0f}_{temperature}K"
    with open(rosetta_results_folder / f"{pickle_name}.pickle", "rb") as f:
        x_values: np.ndarray
        pef_dpef_data: Dict[str, Tuple[np.ndarray, np.ndarray]]
        x_values, pef_dpef_data = pickle.load(f)

    # Optional: add legend-tag for the folder name
    if LEGEND:
        legend_tag = "_legend"
    else:
        legend_tag = ""

    # Create a folder for the figures, if it doesn't exist yet
    out_folder = f"pef_figures_scoreScale{score_scale:.0f}_{temperature}K{legend_tag}{TAG}"
    if not os.path.exists(rosetta_results_folder / out_folder):
        os.mkdir(rosetta_results_folder / out_folder)

    # Set colors for later plotting
    pef_color = "red"
    dpef_color = "orange"
    so_color = "purple"

    # Define legends
    pef_patch = mpatches.Patch(color=pef_color, label="PEF [kcal / mol]")
    single_patches = mpatches.Patch(color=so_color, label="single occurrence")
    dpef_patch = mpatches.Patch(color=dpef_color, label="derivative of PEF [kcal / (mol * degree)]")
    hist_patch = mpatches.Patch(color=HIST_COLOR_CSR, label="CS-Rosetta dihedral angle distribution")

    # Define the x values of the histogram
    hist_width = 2
    hist_x = np.arange(-180, 180 + hist_width, hist_width)
    plot_hist_x = (hist_x[:-1] + hist_x[1:]) / 2

    # Define keys to iterate the data
    keys = list({key[:-4] for key in pef_dpef_data.keys()})
    keys.sort(key=lambda x: int(x.split("-")[0]))

    print("Saving pef figures...")

    # Iterate through the residues to create figures
    for counter, resi_name in enumerate(keys):

        phi_key = resi_name + " PHI"
        psi_key = resi_name + " PSI"

        phi_pef, dphi_pef = pef_dpef_data[phi_key]
        psi_pef, dpsi_pef = pef_dpef_data[psi_key]

        # Create the empty plots
        fig, ax = plt.subplots(2, 2)

        # Define figures with right-side axis on the top
        axr0 = ax[0, 0].twinx()
        axr1 = ax[0, 1].twinx()

        # Plot for PEFs
        axr0.plot(x_values, phi_pef, c=pef_color)
        axr1.plot(x_values, psi_pef, c=pef_color)

        # Scatter plot for single occurrences
        ax[0, 0].scatter(angle_data[phi_key], np.ones_like(angle_data[phi_key]), alpha=0.3, c=so_color, marker="|")
        ax[0, 1].scatter(angle_data[psi_key], np.ones_like(angle_data[psi_key]), alpha=0.3, c=so_color, marker="|")

        # Define the histograms fo CS-Rosetta dihedral angle distribution
        hist_phi_y, _ = np.histogram(angle_data[phi_key], bins=hist_x)
        hist_phi_y = hist_phi_y / np.max(hist_phi_y)
        hist_psi_y, _ = np.histogram(angle_data[psi_key], bins=hist_x)
        hist_psi_y = hist_psi_y / np.max(hist_psi_y)

        # Plot the histograms on a bar diagram
        ax[0, 0].bar(plot_hist_x, hist_phi_y, width=hist_width, color=HIST_COLOR_CSR, alpha=1)
        ax[0, 1].bar(plot_hist_x, hist_psi_y, width=hist_width, color=HIST_COLOR_CSR, alpha=1)

        # Set the title of the whole figure
        sep_pos = resi_name.find("-")
        resi_name_back = f"{resi_name[-3:]}-{resi_name[:sep_pos]}"
        ax[0, 0].set_title(resi_name_back + " $\mathrm{{\phi}}$", fontsize=15, y=1.05)
        ax[0, 1].set_title(resi_name_back + " $\mathrm{{\psi}}$", fontsize=15, y=1.05)

        # Plot the dPEF on the bottom figures
        ax[1, 0].plot(x_values, dphi_pef, c=dpef_color)
        ax[1, 1].plot(x_values, dpsi_pef, c=dpef_color)

        # Set the x-/y-labels of all 4 figures
        ax[0, 0].tick_params(axis="x", labelsize=11)
        ax[0, 0].tick_params(axis="y", labelsize=11, labelcolor=HIST_COLOR_CSR)
        ax[0, 1].tick_params(axis="x", labelsize=11)
        ax[0, 1].tick_params(axis="y", labelsize=11, labelcolor=HIST_COLOR_CSR)
        ax[1, 0].tick_params(axis="x", labelsize=11)
        ax[1, 0].tick_params(axis="y", labelsize=11, labelcolor=dpef_color)
        ax[1, 1].tick_params(axis="x", labelsize=11)
        ax[1, 1].tick_params(axis="y", labelsize=11, labelcolor=dpef_color)
        axr0.tick_params(axis="both", labelsize=11, labelcolor=pef_color)
        axr1.tick_params(axis="both", labelsize=11, labelcolor=pef_color)

        if LEGEND:
            # Optional: set legends describing the data on the figures
            plt.legend(handles=[single_patches, pef_patch, dpef_patch, hist_patch],
                       loc="lower center",
                       bbox_transform=fig.transFigure,
                       ncol=2,
                       bbox_to_anchor=(0.5, 0.01),
                       fontsize=10)

            # Set the text on the bottom describing all the x-axis
            fig.text(0.5, 0.15, "dihedral angle [degree] ", ha="center", fontsize=12)

            plt.subplots_adjust(wspace=0.5, hspace=0.25, bottom=0.25)

        else:
            # Set the text on the bottom describing all the x-axis
            fig.text(0.5, 0.01, "dihedral angle [degree] ", ha="center", fontsize=12)

            plt.subplots_adjust(wspace=0.5, hspace=0.2)

        # Save the figure
        fig.savefig(rosetta_results_folder / f"{out_folder}/{resi_name}.png", dpi=300)

        plt.close(fig)

        print("\r", end="")
        print(progress_bar((counter + 1) / len(keys), 30), end=", ")
        print(f"Filename {resi_name} done...", end="")
    print()
