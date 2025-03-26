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


def main(config_path: Path):

    print("\nformat_scores_csrosetta.py is running:")

    # Load config data
    config = load_config(config_path)

    # Access global variables
    extractedpdbs_path = Path(config.get("ExtractedPDBs_FOLDER"))

    # Check if score.sc exists with the scores of the models
    score_path = Path(extractedpdbs_path.parent / "score.sc")
    if not os.path.exists(score_path):
        sys.exit(f"There is no score.sc file in {score_path}.\n"
                 f"Run e.g. <path of Rosetta3 folder>/main/source/bin/score_jd2.linuxgccrelease "
                 f"-in:file:silent decoys.out")

    # Check if the folder with the pdb files exists
    if not os.path.exists(extractedpdbs_path):
        sys.exit(f"There is no {extractedpdbs_path}")

    out_filename = "name.scores.txt"

    if os.path.exists(extractedpdbs_path / f"{out_filename}"):
        print(f"{extractedpdbs_path}/{out_filename} already exists")
        return

    with open(score_path, "r") as f:
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

        new_line = f"{line[-1]} {line[1]}\n"

        if new_line not in new_score_data:
            new_score_data += new_line

    with open(extractedpdbs_path / f"{out_filename}", "w") as f:
        f.write(new_score_data)


if __name__ == "__main__":
    main()
