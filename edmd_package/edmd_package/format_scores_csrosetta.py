"""
This script is to format the score.sc file: it will only keep the "description" (name of teh model)
and the "total_score" (Rosetta-score of the given model) columns to write a "name.scores.txt" file.
"""

import json
import argparse
from pathlib import Path
import os
import sys

def load_config(config_path: Path):
    """Load configuration from the provided JSON file."""

    with open(config_path, "r") as file:
        config = json.load(file)
    return config


def main():

    print("\nformat_scores_csrosetta.py is running:")

    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=Path, default=Path("EDMD_config.json"),
                        help="A path pointing to an ebMD CONFIGURATION json file.")

    args = parser.parse_args()

    # Load config data
    config = load_config(args.config)

    # Access global variables
    rosetta_results_folder = Path(config.get("ROSETTA_RESULTS_FOLDER"))

    # Check if score.sc exists with the scores of the models
    if not os.path.exists(rosetta_results_folder / "score.sc"):
        sys.exit(f"There is no score.sc file in {rosetta_results_folder}. "
                 f"Run e.g. <path of Rosetta3 folder>/main/source/bin/score_jd2.linuxgccrelease "
                 f"-in:file:silent decoys.out")

    # Check if the folder with the pdb files exists
    if not os.path.exists(rosetta_results_folder / "ExtractedPDBs"):
        sys.exit(f"There is no ExtractedPDBs directory in {rosetta_results_folder}")

    out_filename = "name.scores.txt"

    if os.path.exists(rosetta_results_folder / f"ExtractedPDBs/{out_filename}"):
        print(f"{rosetta_results_folder}/ExtractedPDBs/{out_filename} already exists")
        return

    with open(rosetta_results_folder / "score.sc", "r") as f:
        score_data = f.read()

    # Format the score.sc file
    score_data = list(map(lambda line:
                          list(filter(lambda x:
                                      len(x) != 0 and line[:6] == "SCORE:",
                                      line.replace("\t", " ").split(" "))),
                          score_data.split("\n")))

    score_data = list(filter(lambda line: len(line) != 0, score_data))

    score_data.sort(key=lambda line: line[1], reverse=True)

    new_score_data = ""

    for line in score_data:

        new_score_data += f"{line[-1]} {line[1]}\n"

    with open(rosetta_results_folder / f"ExtractedPDBs/{out_filename}", "w") as f:
        f.write(new_score_data)


if __name__ == "__main__":
    main()
