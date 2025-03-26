"""
Analysing the dihedral angle distribution of the structure ensemble
to define Potential Energy Functions (PEFs) and dPEFs.
"""

import json
import argparse
import sys
import os
from pathlib import Path

from .save_dihedrals import main as save_dihedrals_main
from .visualize_dihedrals import main as visualize_dihedrals_main
from .fit_dihedrals import main as fit_dihedrals_main
from .visualize_pef import main as visualize_pef_main
from .create_tables import main as create_tables_main
from .format_scores_csrosetta import main as format_scores_main


def load_config(config_path: Path):
    """Load configuration from the provided JSON file."""

    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def write_config():
    """Write a blank EDMD_config.json file."""

    new_json = {
        "ExtractedPDBs_FOLDER": "<path to the ExtractedPDBs folder, containing the CS-Rosetta structures>",
        "GMX_FOLDER": "<path to the folder containing the GRO and the processed TOP files for gromacs>",
        "RESI_IDX_SHIFT": 0,
        "VISUALIZE": True,
        "SCORE_SCALE": 10,
        "TEMPERATURE": 310,
        "GRO_FILENAME": "<GRO filename>",
        "PROCESSED_TOP_FILENAME": "<processed TOP filename>"
    }

    current_dir = os.getcwd()

    with open(Path(f"{current_dir}/EDMD_config.json"), "w") as f:
        json.dump(new_json, f)

    print("EDMD_config.json was written.")


def main():

    parser = argparse.ArgumentParser(
        description="CLI for analysing dihedral angles and generating PEFs."
    )
    parser.add_argument(
        "-c", "--config", type=Path, default=Path("EDMD_config.json"),
        help="A path pointing to a config JSON file."
    )
    parser.add_argument(
        "-fn", "--function_name", type=str, default=None,
        help="The name of the individual script you want to call (e.g., save_dihedrals)."
    )
    parser.add_argument(
        "-w", "--write_config", action="store_true",
        help="Write a blank config JSON file."
    )

    args = parser.parse_args()

    if args.write_config:
        write_config()
        sys.exit(0)

    if not os.path.exists(args.config):
        print("No EDMD_config.json was found.")
        print("Use -w to write a blank config file or -c to specify an existing one.")
        sys.exit(1)

    # Mapping function names to their respective functions
    function_map = {
        "save_dihedrals": save_dihedrals_main,
        "fit_dihedrals": fit_dihedrals_main,
        "create_tables": create_tables_main,
        "visualize_dihedrals": visualize_dihedrals_main,
        "visualize_pef": visualize_pef_main,
        "format_scores": format_scores_main
    }

    if args.function_name:
        # Check if function name is valid
        if args.function_name in function_map:
            function_map[args.function_name](args.config)
            sys.exit(0)
        else:
            print(f"Invalid function name '{args.function_name}'!\n")
            print("Available functions:")
            for name in function_map.keys():
                print(f"  - {name}")
            sys.exit(1)

    # Default: Run full workflow
    save_dihedrals_main(args.config)
    fit_dihedrals_main(args.config)
    create_tables_main(args.config)

    # Load config and check if visualization is enabled
    config = load_config(args.config)
    if config.get("VISUALIZE", False):
        visualize_dihedrals_main(args.config)
        visualize_pef_main(args.config)

    print("\nAll done!")


if __name__ == "__main__":
    main()
