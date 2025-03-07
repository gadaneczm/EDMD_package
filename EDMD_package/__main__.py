"""
Analysing the dihedral angle distribution of the structure ensemble
to define Potential Energy Functions (PEFs) and dPEFs.
"""

import json
import argparse
import sys
from pathlib import Path

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

    parser = argparse.ArgumentParser(
        description="CLI for analysing dihedral angles and generating PEFs."
    )
    parser.add_argument(
        "-c", "--config", type=Path, default=Path("EDMD_config.json"),
        help="A path pointing to an ebMD CONFIGURATION json file."
    )
    parser.add_argument(
        "-fn", "--function_name", type=str, default=None,
        help="The name of the individual script you want to call (e.g., save_dihedrals)."
    )

    args = parser.parse_args()

    # Mapping function names to their respective functions
    function_map = {
        "save_dihedrals": save_dihedrals_main,
        "fit_dihedrals": fit_dihedrals_main,
        "create_tables": create_tables_main,
        "visualize_dihedrals": visualize_dihedrals_main,
        "visualize_pef": visualize_pef_main,
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
