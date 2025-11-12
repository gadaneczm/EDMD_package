from pathlib import Path
from Bio.PDB.PDBParser import PDBParser
from loco_hd import *


def prat_to_pra(prat: PrimitiveAtomTemplate) -> PrimitiveAtom:

    resi_id = prat.atom_source.source_residue
    resname = prat.atom_source.source_residue_name
    source = f"{resi_id[2]}/{resi_id[3][1]}-{resname}"

    return PrimitiveAtom(
        prat.primitive_type,
        source,  # this is the tag field!
        prat.coordinates
    )


def main():

    ref = PDBParser(QUIET=True).get_structure("s1",
                                              "/rhome/PROTMOD/gadaneczm/EDMD_benchmark/fbp-fmn/1flm_ref.pdb")
    traj = PDBParser(QUIET=True).get_structure("s2",
                                               "/rhome/PROTMOD/gadaneczm/EDMD_processing/results_gmx_runs/"
                                               "fbp-fmn/csra1_fbp-fmn_opc_edmd_310K_1/results_data/"
                                               "csra1_fbp-fmn_opcMD_mcc_dt20ns.pdb")

    primitive_assigner = PrimitiveAssigner(Path("/rhome/PROTMOD/gadaneczm/PycharmProjects/loco_hd/all_atom_fmn_p.config.json"))

    tag_pairing_rule = TagPairingRule({"accept_same": False})

    lchd = LoCoHD(primitive_assigner.all_primitive_types, tag_pairing_rule=tag_pairing_rule)

    pra_templates1 = primitive_assigner.assign_primitive_structure(ref)

    pra1 = list(map(prat_to_pra, pra_templates1))

    scores = []

    for idx, frame in enumerate(traj):

        pra_templates2 = primitive_assigner.assign_primitive_structure(frame)

        anchor_pairs = [
            (idx_a, idx_b)
            for idx_a, prat_a in enumerate(pra_templates1)
            if prat_a.primitive_type == "FMN_P"
            for idx_b, prat_b in enumerate(pra_templates2)
            if prat_b.primitive_type == "FMN_P"
        ]

        pra2 = list(map(prat_to_pra, pra_templates2))

        lchd_scores = lchd.from_primitives(
            pra1,
            pra2,
            anchor_pairs,
            10.  # upper distance cutoff at 10 angströms
        )

        scores.append(lchd_scores[0])

    print(scores)


if __name__ == "__main__":
    main()
