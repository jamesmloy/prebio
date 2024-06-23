from argparse import ArgumentParser
from pathlib import Path
from prebio import *
import numpy as np
import pickle as pkl
import os
from json import dumps

log_level = os.environ.get("LOG_LEVEL", "INFO")


def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--in-folder",
        type=Path,
        help="Folder containing the cif files to run predictions on.",
        required=True,
    )

    parser.add_argument(
        "--out-folder",
        type=Path,
        help="The folder where the pickle files will be saved",
        required=True,
    )

    parser.add_argument(
        "--channel-first",
        action="store_true",
        help="Channel first data ordering (default channel last, ignored for graphs)",
    )

    parser.add_argument(
        "--graph",
        action="store_true",
        help="Generate graph based microenvironment",
    )

    return parser.parse_args()


def make_collection_plus(
    ac: AtomicCollection, ss: dict, additional_columns: list[str]
) -> dict:
    target_atom_idx = ac.get_atoms_at(ss["chain_id"], ss["label"], ss["res_seq_num"])
    atoms = ac.get_atoms()

    return {
        "label": ss["label"],
        "atomic_collection": ac.as_json(additional_columns),
        "target": target_atom_idx,
        "target_alpha_carbon": next(
            (i for i in target_atom_idx if atoms[i].atom_name() == "CA"),
            None,
        ),
        "snapshot": ss,
    }


def save_voxelized(
    in_folder: Path, out_folder: Path, the_file: Path, channel_first: bool
):    
    try:
        snapshots = snapshots_from_cif(the_file)
        vp = VoxelPipeline(
            snapshots=snapshots,
            cif_folder=in_folder,
            channel_first=channel_first,
        )
        boxes = [b for b in vp()]
        to_pickle = {"snapshots": snapshots, "boxes": np.asarray(boxes)}

        suffix = "_ch_first" if channel_first else "_ch_last"
        pkl_file = out_folder / (the_file.name.split(".", 1)[0] + suffix + ".pkl")
        with pkl_file.open("wb") as f:
            pkl.dump(to_pickle, f)
    except Exception as e:
        print(
            f"Failed to generate boxes for {the_file.name} because an exception occurred\n{e}"
        )
    else:
        print(f"Saved pickle file for {the_file.name} to {pkl_file.name}.")


def save_graph(in_folder: Path, out_folder: Path, the_file: Path):
    try:
        snapshots = snapshots_from_cif(the_file)
        ac_maker = AtomicCollectionMaker(cif_folder=in_folder, filter_radius=10.0)
        ac_list = [
            dumps(
                make_collection_plus(
                    ac_maker.make_collection(ss),
                    ss,
                    ac_maker.additional_params,
                )
            )
            for ss in snapshots
        ]
        jsonl_file = out_folder / (the_file.name.split(".", 1)[0] + ".jsonl")
        jsonl_file.write_text("\n".join(ac_list))
    except Exception as e:
        print(
            f"Failed to generate boxes for {the_file.name} because an exception occurred\n{e}"
        )
    else:
        print(f"Saved jsonl file for {the_file.name} to {jsonl_file.name}.")


def main(in_folder: Path, out_folder: Path, channel_first: bool, graph: bool):
    out_folder.mkdir(exist_ok=True, parents=True)

    for the_file in in_folder.iterdir():
        if the_file.suffix == ".cif":
            print(f"Processing file {the_file.name}")
            if graph:
                save_graph(
                    in_folder=in_folder, out_folder=out_folder, the_file=the_file
                )
            else:
                save_voxelized(
                    in_folder=in_folder,
                    out_folder=out_folder,
                    the_file=the_file,
                    channel_first=channel_first,
                )


if __name__ == "__main__":
    args = vars(get_args())
    main(**args)
