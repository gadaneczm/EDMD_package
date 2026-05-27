from pathlib import Path
import re

FOLDER = Path("/rhome/PROTMOD/gadaneczm/Poky/Lists/")
PATTERN = re.compile(r"^kras-wt.*VW_uv\.list$")
SEQ_LEN = 169

def main():
    for file_path in FOLDER.iterdir():
        if file_path.is_file() and PATTERN.match(file_path.name):
            print(f"Matched: {file_path.name}")

            with open(file_path, "r") as f:
                data = f.read()

                data = list(map(lambda line: list(filter(lambda x: len(x) != 0, line.split(" "))), data.split("\n")))

                last_resi = data[-2][0]
                last_resi_nb = int(last_resi[1:-3])
                correction = SEQ_LEN - last_resi_nb

                resi_counter = 0

                new_text = ""

                for line in data:

                    if len(line) == 0:
                        new_text += "\n"
                        continue

                    if line[0] == "Assignment":
                        new_text += f"\t{line[0]}\t{line[1]}\t{line[2]}\t{line[3]}\n"
                        continue

                    resi_counter += 1

                    resi = line[0]

                    resi_type = resi[0]
                    resi_nb = int(resi[1:-3])
                    atoms = resi[-3:]

                    new_resi_nb = resi_nb + correction

                    while new_resi_nb > resi_counter:
                        new_text += f"\tX{resi_counter}N-H\t0\t0\t0\n"
                        resi_counter += 1

                    new_text += f"\t{resi_type}{new_resi_nb}{atoms}\t{line[1]}\t{line[2]}\t{line[3]}\n"

            with open(FOLDER / f"{file_path.stem}.new.list", "w") as f:
                f.write(new_text)


if __name__ == "__main__":
    main()