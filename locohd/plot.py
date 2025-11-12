import pickle
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

TITLE = "LoCoHD of FMN's Phosphate anchor atom\n in FBP NOE-RASREC-Rosetta refined trajectories"
JSON = "LoCoHD_FMN_dt1ns_v3.pickle"
SOURCES = [
    f"/rhome/PROTMOD/gadaneczm/EDMD_processing/results_gmx_runs/fbp-fmn/csrr1_fbp-fmn_opc_edmd_310K_1/{JSON}",
    f"/rhome/PROTMOD/gadaneczm/EDMD_processing/results_gmx_runs/fbp-fmn/csrr1_fbp-fmn_opc_edmd_310K_2/{JSON}",
    f"/rhome/PROTMOD/gadaneczm/EDMD_processing/results_gmx_runs/fbp-fmn/csrr1_fbp-fmn_opc_fforig_310K_1/{JSON}",
    f"/rhome/PROTMOD/gadaneczm/EDMD_processing/results_gmx_runs/fbp-fmn/csrr1_fbp-fmn_opc_fforig_310K_2/{JSON}"
    ]
LABELS = ["EDMD-1", "EDMD-2", "MD-1", "MD-2"]
COLORS = ["deepskyblue", "royalblue", "limegreen", "green"]
ATOMS = ["M51", "FMN123", "FMN123"]
DATA_IDX = 0

def main():

    fig, ax = plt.subplots()

    for idx in [0,1]:

        with open(Path(SOURCES[idx]), "rb") as f:
            ref_anchors, scores, times = pickle.load(f)

        times = [int(time / 1000) for time in times]

        ax.plot(times, [score[DATA_IDX] for score in scores], color=COLORS[idx], label=LABELS[idx])

    file_name = f"LoCoHD_csrr-FBP_EDMD_{ATOMS[DATA_IDX]}-{ref_anchors[DATA_IDX][0]}_v2"

    ax.legend(loc="upper left", fontsize=18)

    ax.set_ylim(0, 0.12)

    ax.tick_params(axis='both', labelsize=16)

    plt.tight_layout()

    plt.show()

    fig.savefig(f"/rhome/PROTMOD/gadaneczm/EDMD_benchmark/{file_name}.png", dpi=300)

    print()


if __name__ == "__main__":
    main()