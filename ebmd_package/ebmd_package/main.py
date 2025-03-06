"""
Analysing the dihedral angle distribution of the structure ensemble
to define Potential Energy Functions (PEFs) and dPEFs.
"""

import json
import argparse
import sys
from pathlib import Path
from typing import Optional

from .save_dihedrals import main as save_dihedrals_main
from .visualize_dihedrals import main as visualize_dihedrals_main
from .create_tables import main as create_tables_main
from .fit_dihedrals import main as fit_dihedrals_main
from .visualize_pef import main as visualize_pef_main


def load_config(config_path: Path):
    """Load configuration from the provided JSON file."""

    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=Path, default=Path("ebMD_config.json"),
                        help="A path pointing to an ebMD CONFIGURATION json file.")
    parser.add_argument("-fn", "--function_name", default=None, type=Optional[str],
                        help="The name of the individual script you want to call, e.g. save_dihedrals.")

    args = parser.parse_args()

    if args.function_name == "save_dihedrals":
        save_dihedrals_main(args.config)
        sys.exit(0)
    elif args.function_name == "fit_dihedrals":
        fit_dihedrals_main(args.config)
        sys.exit(0)
    elif args.function_name == "create_tables":
        create_tables_main(args.config)
        sys.exit(0)
    elif args.function_name == "visualize_dihedrals":
        visualize_dihedrals_main(args.config)
        sys.exit(0)
    elif args.function_name == "visualize_pef":
        visualize_pef_main(args.config)
        sys.exit(0)

    elif args.function_name is not None:
        print("Invalid function name!")
        sys.exit(1)

    # Calling scripts to analyse the structure ensemble and create PEFs
    save_dihedrals_main(args.config)
    fit_dihedrals_main(args.config)
    create_tables_main(args.config)

    # Loading global variable
    config = load_config(args.config)
    visualize: bool = config.get("VISUALIZE")

    # Optional visualization
    if visualize:
        visualize_dihedrals_main(args.config)
        visualize_pef_main(args.config)

    print("\nAll done!")


if __name__ == "__main__":
    main()
