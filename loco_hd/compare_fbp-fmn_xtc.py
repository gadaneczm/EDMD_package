import pickle
from typing import List, Dict
from pathlib import Path
from time import time
import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe
from MDAnalysis.coordinates.base import Timestep
from loco_hd import *  # Use the "locohd" interpreter!

XRAY_REF_PATH = Path("/rhome/PROTMOD/gadaneczm/EDMD_benchmark/fbp-fmn/1flm_ref.pdb")
SOURCE_DIR = Path("/rhome/PROTMOD/gadaneczm/EDMD_processing/results_gmx_runs/fbp-fmn/csrr1_fbp-fmn_opc_fforig_310K_2")
TRAJECTORY_PATH = SOURCE_DIR / "results_data/csrr1_fbp-fmn_opcMD_mcc.xtc"
STRUCTURE_PATH = SOURCE_DIR / "results_data/csrr1_fbp-fmn_opcMD_mcc_dump0.pdb"
PRIMITIVE_TYPING_SCHEME_PATH = Path("/rhome/PROTMOD/gadaneczm/PycharmProjects/loco_hd/all_atom_fmn.config.json")
OUT_NPY = "LoCoHD_FMN_dt1ns_v3.pickle"

class MDPrimitiveAssigner(PrimitiveAssigner):
    def assign_from_universe(self, frame: Timestep, universe: Universe) -> List[PrimitiveAtomTemplate]:

        out = list()

        resi: mda.core.groups.Residue
        for resi in universe.residues:

            resi_name = resi.resname
            resi_id = ("prot", 0, "A", (" ", resi.ix, " "))

            for tse in self.scheme:

                if not tse.match_resi(resi_name):
                    continue

                atom_names = list()
                atom_coords = list()

                atom: mda.core.groups.Atom
                for atom in resi.atoms:

                    if not tse.match_atom(atom.name):
                        continue

                    atom_names.append(atom.name)
                    atom_coords.append(frame.positions[atom.ix])

                if tse.atom_counter == "any":
                    pass
                elif tse.atom_counter == len(atom_coords):
                    pass
                else:
                    continue

                centroid = np.mean(atom_coords, axis=0)
                pras = PrimitiveAtomSource(resi_id, resi_name, atom_names)
                prat = PrimitiveAtomTemplate(tse.primitive_type, centroid, pras)
                out.append(prat)

        return out


def prat_to_pra(prat: PrimitiveAtomTemplate) -> PrimitiveAtom:

    resi_id = prat.atom_source.source_residue
    resname = prat.atom_source.source_residue_name
    source = f"{resi_id[2]}/{resi_id[3][1]}-{resname}"

    return PrimitiveAtom(
        prat.primitive_type,
        source,  # this is the tag field!
        prat.coordinates
    )


def select_anchors(pra_template_list: List[PrimitiveAtomTemplate]) -> Dict[str, int]:

    out = dict()
    for prat_idx, prat in enumerate(pra_template_list):

        if prat.primitive_type == "FMN_P":
            out["P"] = prat_idx
        elif (prat.atom_source.source_residue_name == "FMN"
              and
              prat.atom_source.source_atom[0] == "C2"):
            out["C2"] = prat_idx
        elif (prat.atom_source.source_residue_name == "MET"
              and
              prat.atom_source.source_residue[3][1] == 50  # indexing starts from 0
              and
              prat.atom_source.source_atom[0] == "SD"):
            out["SD"] = prat_idx

    return out


def main():
    print(f"Running in {SOURCE_DIR}")

    # Read the structure files and the trajectory
    ref = mda.Universe(str(XRAY_REF_PATH), str(XRAY_REF_PATH))
    universe = mda.Universe(str(STRUCTURE_PATH), str(TRAJECTORY_PATH))

    # Define primitive typing scheme
    primitive_assigner = PrimitiveAssigner(PRIMITIVE_TYPING_SCHEME_PATH)

    # Set pairing rule
    tag_pairing_rule = TagPairingRule({"accept_same": True})

    # Define the LoCoHD instance
    lchd = LoCoHD(primitive_assigner.all_primitive_types, tag_pairing_rule=tag_pairing_rule)

    # Define the primitive typing scheme
    primitive_assigner = MDPrimitiveAssigner(PRIMITIVE_TYPING_SCHEME_PATH)

    # Get primitive atoms for the reference
    pra_templates_ref = primitive_assigner.assign_from_universe(ref.trajectory[0], ref)
    ref_anchors = list(select_anchors(pra_templates_ref).items())
    pra_ref = list(map(prat_to_pra, pra_templates_ref))

    scores, times = [], []

    delta_frame = 250  # dt = 1 ns

    full_trajectory_length = len(universe.trajectory[::delta_frame])
    start_real_time = time()

    frame: Timestep
    for frame_idx, frame in enumerate(universe.trajectory[::delta_frame]):

        # Get the primitive atoms for the current frame
        pra_templates_frame = primitive_assigner.assign_from_universe(frame, universe)
        frame_anchors = select_anchors(pra_templates_frame)
        pra_frame = list(map(prat_to_pra, pra_templates_frame))

        anchor_pairs = [
            (idx1, frame_anchors[key]) for key, idx1 in ref_anchors
        ]

        # anchor_pairs = [
        #     (idx_a, idx_b)
        #     for idx_a, prat_a in enumerate(pra_templates_ref)
        #     if prat_a.primitive_type == "FMN_P"
        #     for idx_b, prat_b in enumerate(pra_templates_frame)
        #     if prat_b.primitive_type == "FMN_P"
        # ]

        lchd_scores = lchd.from_primitives(
            pra_ref,
            pra_frame,
            anchor_pairs,
            10.  # upper distance cutoff at 10 angströms
        )

        scores.append(lchd_scores)
        times.append(frame.time)

        # Print out time statistics
        current_real_time = time()
        eta = full_trajectory_length - frame_idx  # number of remaining frames
        eta *= (current_real_time - start_real_time) / (frame_idx + 1)  # current time / frame rate
        print(f"\r{(frame_idx + 1) / full_trajectory_length:.1%} at time {frame.time:.0f} ps. ETA: {eta:.1f} s.", end="")

    with open(SOURCE_DIR / OUT_NPY, "wb") as f:
        pickle.dump((ref_anchors, scores, times), f)


if __name__ == "__main__":
    main()
