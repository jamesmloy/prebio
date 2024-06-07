from collections import namedtuple


Residue = namedtuple("Residue", ["asym", "seq", "comp", "label_seq"])

Ligand = namedtuple("Ligand", ("comp_id", "asym_id", "seq_id"))

Atom = namedtuple("Atom", ("atom_id", "x", "y", "z", "type_symbol", "comp_id"))

AtomSiteAtom = namedtuple(
    "CifAtom",
    (
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_is_or_equiv",
        "pdbx_formal_charge",
        "auth_seq_id",
        "auth_comp_id",
        "auth_asym_id",
        "auth_atom_id",
        "pdbx_PDB_model_num",
    ),
)

ChemComp = namedtuple("ChemComp", ("id", "type", "name", "formula", "mass"))
# -- Deprecated
LigandData = namedtuple("LigandData", ("file", "ligand", "H_insert_idx"))

