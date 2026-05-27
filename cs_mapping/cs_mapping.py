from pathlib import Path
from typing import List
import math
import statistics
import matplotlib.pyplot as plt

MUTANT = "g12d"
FOLDER = Path("/rhome/PROTMOD/gadaneczm/kras_foldamers/kras-g12d_ph-foldamers_2022/")
PATTERN = f"kras-{MUTANT}_ph-*20220117.new.list"


def read_data(file_path: Path) -> (List[float], List[float], List[float]):

    with open(file_path, "r") as f:
        data = f.read()

    data = list(filter(lambda line: len(line) != 0 and "Data" not in line, data.split("\n")))
    data = list(map(lambda line: list(filter(lambda x: len(x) != 0, line.split("\t"))), data))

    resi_list, cs_n_list, cs_h_list, height_list = [], [], [], []

    last_resi_height = int(data[-1][3])  # for data height normalization

    for line in data:
        resi_list.append(line[0])
        cs_n_list.append(float(line[1]))
        cs_h_list.append(float(line[2]))
        height_list.append(int(line[3]) / last_resi_height)  # normalized height

    return resi_list, cs_n_list, cs_h_list, height_list


def get_csp(ref_n: float, ref_h: float, ph_n: float, ph_h: float) -> float:

    if ref_n == 0 or ref_h == 0 or ph_n == 0 or ph_h == 0:  # if any cs is missing
        return 0

    delta_n = ph_n - ref_n
    delta_h = ph_h - ref_h

    csp = math.sqrt(0.5 * (delta_h**2 + 0.14 * delta_n**2))

    return csp


def get_height_reduction(ref_height: int, ph_height: int) -> float:

    if ref_height == 0 and ph_height == 0:  # if there was no assignment originally
        return 0

    elif ref_height != 0 and ph_height == 0:  # if the peak is disappeared
        return 100

    height_red = abs(1 - (ph_height / ref_height)) * 100  # peak height reduction in %, relative to the reference

    return height_red


def get_figures(csp_list: List[float], height_red_list: List[float], file_path: Path):

    fig, ax = plt.subplots(2, 1)

    csp_std = statistics.stdev(csp_list)

    for line in [csp_std * 2, csp_std * 3]:
        ax[0].axhline(y=line,
                      color="black",
                      linestyle="--",
                      linewidth=1)

    ax[0].bar(range(1, len(csp_list) + 1), csp_list)
    ax[0].set_ylabel("CSP [ppm]")

    for line in [30, 50]:
        ax[1].axhline(y=line,
                      color="black",
                      linestyle="--",
                      linewidth=1)

    ax[1].bar(range(1, len(csp_list) + 1), height_red_list, color="tab:orange")
    ax[1].set_ylabel("Rel. intensity\nreduction [%]")
    ax[1].set_xlabel("Residue index")



    name = file_path.name.split(".")[0]
    new_path = file_path.parent / "results" / name

    fig.suptitle(f"{name}")

    plt.tight_layout()

    #plt.show()

    fig.savefig(new_path, dpi=300)

    plt.close()


def main():

    ref_file = FOLDER / f"kras-{MUTANT}_ref_uv_20220117.new.list"

    ref_resi_list, ref_n_list, ref_h_list, ref_height_list = read_data(ref_file)

    results_txt = ""
    csp_csv = ""
    height_csv = ""

    for file_path in FOLDER.glob(PATTERN):
        print(f"Matched: {file_path.name}")

        ph_resi_list, ph_n_list, ph_h_list, ph_height_list = read_data(file_path)

        csp_list, height_red_list = [], []

        results_txt += f"\n{file_path.name}\n"
        csp_csv += f"\n{file_path.name}"
        height_csv += f"\n{file_path.name}"

        for (ref_n, ref_h, ref_height, ph_n, ph_h, ph_height) in zip(ref_n_list, ref_h_list, ref_height_list,
                                                             ph_n_list, ph_h_list, ph_height_list):

            csp = get_csp(ref_n, ref_h, ph_n, ph_h)

            csp_list.append(csp)

            csp_csv += f",{csp}"

            height = get_height_reduction(ref_height, ph_height)

            height_red_list.append(height)

            height_csv += f",{height}"

        csp_std = statistics.stdev(csp_list)

        mid_csp, big_csp = "", ""
        mid_height_red, big_height_red = "", ""

        for resi, csp, height_red in zip(ref_resi_list, csp_list, height_red_list):

            if csp > csp_std * 3:
                big_csp += f"{resi[:-3]} "
            elif csp > csp_std * 2:
                mid_csp += f"{resi[:-3]} "

            if height_red > 50:
                big_height_red += f"{resi[:-3]} "
            elif height_red > 30:
                mid_height_red += f"{resi[:-3]} "

        results_txt += f"\tfast exchange - shift\n"
        results_txt += f"\t\tbig shift (>3 std): {big_csp}\n"
        results_txt += f"\t\tmid shift (>2 std): {mid_csp}\n"

        results_txt += f"\tmid exchange - broadening\n"
        results_txt += f"\t\tbig broadening (>50%): {big_height_red}\n"
        results_txt += f"\t\tmid broadening (>30%): {mid_height_red}\n"

        get_figures(csp_list, height_red_list, file_path)

    with open(FOLDER / "results" / f"kras-{MUTANT}_results.txt" , "w") as f:
        f.write(results_txt)

    with open(FOLDER / "results" / f"kras-{MUTANT}_csp_20220117.csv" , "w") as f:
        f.write(csp_csv)

    with open(FOLDER / "results" / f"kras-{MUTANT}_height_20220117.csv" , "w") as f:
        f.write(height_csv)


if __name__ == "__main__":
    main()
