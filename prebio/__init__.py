from prebio.prebioimpl import (
    VoxelPipeline as PipelineImpl,
    from_protein_chain_residue,
    AtomicCollection as AtomicCollectionCpp,
)
import json
import numpy as np
from pathlib import Path
import gemmi
from pydantic import BaseModel


AMINO_ACIDS = set(
    [
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
    ]
)


def snapshots_from_cif(cif: Path):
    cif_block = gemmi.cif.read(str(cif))[0]

    model = gemmi.make_structure_from_block(cif_block)[0]

    snapshots = []
    for chain in model:
        for residue in chain.first_conformer():
            if (
                residue.name.upper() in AMINO_ACIDS
                and residue.find_atom("CA", "*") is not None
            ):
                snapshots.append(
                    {
                        "type": "FILE_CHAIN_RESIDUE",
                        "filename": cif.name,
                        "chain_id": chain.name,
                        "res_seq_num": residue.seqid.num,
                        "label": residue.name.upper(),
                    }
                )

    return snapshots


def snapshots_from_cif_folder(cif_folder: Path):
    snapshots = []
    for the_file in cif_folder.iterdir():
        if the_file.suffix == ".cif":
            snapshots.extend(snapshots_from_cif(the_file))
    return snapshots


class VoxelPipeline(PipelineImpl):
    def __init__(self, snapshots: list[dict], cif_folder: Path, channel_first: bool):
        str_snapshots = [json.dumps(ss) for ss in snapshots]
        super().__init__(str_snapshots, str(cif_folder), channel_first)

    def __call__(self):
        while self.has_more_data():
            res_tuple = self.next()
            yield np.asarray(res_tuple[0])


class AtomicCollection:
    def __init__(self, ac_cpp: AtomicCollectionCpp):
        self.ac = ac_cpp

    def get_atoms(self):
        return self.ac.get_atoms()

    def get_atoms_at(self, chain_id, label, seq_num):
        return self.ac.get_atoms_at(chain_id, label, seq_num)

    def as_json(self, addl_params=[]):
        ac_json = {}
        ac_json["x"] = [a.x() for a in self.get_atoms()]
        ac_json["y"] = [a.y() for a in self.get_atoms()]
        ac_json["z"] = [a.z() for a in self.get_atoms()]
        ac_json["element"] = [a.element() for a in self.get_atoms()]
        ac_json["residue_name"] = [a.residue_name() for a in self.get_atoms()]
        ac_json["chain_id"] = [a.chain_id() for a in self.get_atoms()]
        ac_json["res_seq_num"] = [a.residue_sequence_number() for a in self.get_atoms()]
        ac_json["atom_name"] = [a.atom_name() for a in self.get_atoms()]

        col_prefix = "_atom_site."
        prefixed_params = []
        for ap in addl_params:
            if ap[: len(col_prefix)] != col_prefix:
                prefixed_params.append(col_prefix + ap)
            else:
                prefixed_params.append(ap)

        for key in prefixed_params:
            ac_json[key] = [a[key] for a in self.get_atoms()]

        return ac_json


class AtomicCollectionMaker(BaseModel):
    cif_folder: Path
    filter_radius: float
    include_target: bool = True
    include_hetatm: bool = True
    include_water: bool = False
    additional_params: list[str] = ["_atom_site.fw2_charge", "_atom_site.FreeSASA_value"]

    def make_collection(
        self, snapshot: dict, doc: gemmi.cif.Document | None = None
    ) -> AtomicCollection:
        if doc is None:
            filename = snapshot.get("filename")
            if not filename:
                raise RuntimeError(
                    "Could not make collection for snapshot "
                    f"{snapshot} because not document was provided "
                    "and the snapshot does not have the filename"
                )
            doc = gemmi.cif.read_file(str(self.cif_folder / filename))

        ac = from_protein_chain_residue(
            snapshot,
            doc,
            self.filter_radius,
            self.additional_params,
            self.include_target,
            self.include_hetatm,
            self.include_water,
        )

        return AtomicCollection(ac_cpp=ac)
