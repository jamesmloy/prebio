"""danny305 7/7/14"""

import gemmi
import json
import requests

from subprocess import run
from collections import Counter, namedtuple
from pathlib import Path
from typing import List
from tempfile import TemporaryDirectory

import chargefw2_python

from ApChem import ChimeraXException, ChargeFW2Exception, FreeSASAException

from ApChem.cif import AtomSiteAtom, Residue, ChemComp
from ApChem.utils import rename_cif_by_split
from ApChem.cif.doc_split import run_as_seprate_model_docs


cation_formal_charges_radii = {
    "LI": ("1.00", "1.82"),
    "NA": ("1.00", "2.27"),
    "K": ("1.00", "2.75"),
    "MG": ("2.00", "1.73"),
    "CA": ("2.00", "2.31"),
    "ZN": ("2.00", "1.39"),
    "FE": ("3.00", "1.56"),
    "FE2": ("2.00", "1.56"),
    "CU": ("2.00", "1.40"),
    "CU1": ("1.00", "1.40"),
    "CU3": ("3.00", "1.40"),
    "CR": ("3.00", "1.66"),
    "RU": ("3.00", "1.78"),
    "RH": ("1.00", "1.73"),
    "MO": ("0.00", "1.90"),
    "AG": ("1.00", "1.72"),
    "AU": ("1.00", "1.66"),
    "AU3": ("3.00", "1.66"),
    "CO": ("2.00", "1.52"),
    "NI": ("2.00", "1.49"),
    "MN": ("2.00", "1.61"),
    "MN3": ("3.00", "1.61"),
    "IR": ("4.00", "1.80"),
    "IR3": ("3.00", "1.80"),
    "PB": ("2.00", "2.02"),
    "HG": ("2.00", "1.55"),
}

ATOM_SITE_TAGS = [
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
    "B_iso_or_equiv",
    "pdbx_formal_charge",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_PDB_model_num",
]


def find_mmcif_tag_idx(mmcif_category, tag):
    tags = [tag.split(".")[1] for tag in mmcif_category.tags]

    try:
        return tags.index(tag)

    except ValueError:
        return None


def fetch_cif_from_api(pdb: str, cif_dir=Path("./cif"), pdb_redo=True, redo_fail=False):
    assert len(pdb) == 4

    try:
        assert pdb_redo is True

        print("Fetching CIF file from PDB-REDO...")

        respond = requests.get(f"https://pdb-redo.eu/db/{pdb}/{pdb}_final.cif")
        if respond.ok:
            cif_file = cif_dir / f"{pdb}_final.cif"

            with cif_file.open("wt") as f:
                f.write(respond.text)

            reorder_pdb_redo_columns(cif_file)

        else:
            msg = f"Unable to retrieve cif file from pdb-redo: {pdb}."
            if redo_fail:
                raise RuntimeError(msg)

            print(msg + " Attempting RCSB...")

            raise FileNotFoundError(msg)

    except (AssertionError, FileNotFoundError):
        print("Fetching CIF file from RCSB...")

        cif_file = cif_dir / f"{pdb}.cif"

        respond = requests.get(f"https://files.rcsb.org/download/{pdb}.cif")
        if respond.ok:
            with cif_file.open("wt") as f:
                f.write(respond.text)

        else:
            raise FileNotFoundError(f"Unable to retrieve cif file from RCSB: {pdb}")

    return cif_file


def reorder_pdb_redo_columns(redo_cif):
    doc = gemmi.cif.read_file(str(redo_cif))

    block = doc.sole_block()

    atom_site = block.find_mmcif_category("_atom_site.")

    assert (
        atom_site.width() == 21
    ), f"atom_site does not have 21 columns: {atom_site.width()}"

    group_idx = find_mmcif_tag_idx(atom_site, "group_PDB")
    id_idx = find_mmcif_tag_idx(atom_site, "id")
    type_idx = find_mmcif_tag_idx(atom_site, "type_symbol")
    lab_atom_id_idx = find_mmcif_tag_idx(atom_site, "label_atom_id")
    lab_alt_id_idx = find_mmcif_tag_idx(atom_site, "label_alt_id")
    lab_comp_id_idx = find_mmcif_tag_idx(atom_site, "label_comp_id")
    lab_asym_id_idx = find_mmcif_tag_idx(atom_site, "label_asym_id")
    lab_entity_id_idx = find_mmcif_tag_idx(atom_site, "label_entity_id")
    lab_seq_id_idx = find_mmcif_tag_idx(atom_site, "label_seq_id")
    ins_idx = find_mmcif_tag_idx(atom_site, "pdbx_PDB_ins_code")
    x_idx = find_mmcif_tag_idx(atom_site, "Cartn_x")
    y_idx = find_mmcif_tag_idx(atom_site, "Cartn_y")
    z_idx = find_mmcif_tag_idx(atom_site, "Cartn_z")
    occupancy_idx = find_mmcif_tag_idx(atom_site, "occupancy")
    b_factor_idx = find_mmcif_tag_idx(atom_site, "B_iso_or_equiv")
    formal_charge_idx = find_mmcif_tag_idx(atom_site, "pdbx_formal_charge")
    auth_seq_id_idx = find_mmcif_tag_idx(atom_site, "auth_seq_id")
    auth_comp_id_idx = find_mmcif_tag_idx(atom_site, "auth_comp_id")
    auth_asym_id_idx = find_mmcif_tag_idx(atom_site, "auth_asym_id")
    auth_atom_id_idx = find_mmcif_tag_idx(atom_site, "auth_atom_id")
    model_num_idx = find_mmcif_tag_idx(atom_site, "pdbx_PDB_model_num")

    reordered_atom_site = [
        AtomSiteAtom(
            atom[group_idx],
            atom[id_idx],
            atom[type_idx],
            atom[lab_atom_id_idx],
            atom[lab_alt_id_idx],
            atom[lab_comp_id_idx],
            atom[lab_asym_id_idx],
            atom[lab_entity_id_idx],
            atom[lab_seq_id_idx],
            atom[ins_idx],
            atom[x_idx],
            atom[y_idx],
            atom[z_idx],
            atom[occupancy_idx],
            atom[b_factor_idx],
            atom[formal_charge_idx],
            atom[auth_seq_id_idx],
            atom[auth_comp_id_idx],
            atom[auth_asym_id_idx],
            atom[auth_atom_id_idx],
            atom[model_num_idx],
        )
        for atom in atom_site
    ]

    formatted_atom_site = block.init_mmcif_loop("_atom_site.", ATOM_SITE_TAGS)

    for atom in reordered_atom_site:
        formatted_atom_site.add_row(atom)

    block.set_pair("_apchem.pdb_redo_correct_column_order", "True")

    doc.write_file(str(redo_cif), gemmi.cif.Style.Pdbx)

    return redo_cif


"""
chimera_addh.py
"""


def run_chimeraX(inp_cif: Path, out_cif: Path = None):
    """Add hydrogens to cif file using ChimeraX.

    :param inp_cif: cif file you want to add hydrogens to
    :type inp_cif: Path
    :raises ValueError: if in files do not have cif extension
    :raises ValueError: if any of the input files do not exist
    """

    inp_cif = Path(inp_cif)

    addh_cmd_fmt = "open {}; addh; save {} format mmcif; exit"

    if inp_cif.suffix != ".cif":
        raise ValueError(f"File {inp_cif} is not a cif file.")

    if not inp_cif.exists():
        raise ValueError(f"File {inp_cif} does not exist")

    if out_cif is None:
        out_cif = inp_cif.with_suffix(".chimX.cif")

    # This line fails in the singularity image
    # Path('/non/existent/directory/.local/share/ChimeraX').mkdir(mode=0o777, parents=True, exist_ok=True)

    # chimerax doesnt like saving cif files if the protein has
    # chains with ids that are non-alphanumeric.  this just removes
    # all the non alnum characters.
    in_doc = gemmi.cif.read_file(str(inp_cif))
    block = in_doc.sole_block()
    atom_site_table = block.find_mmcif_category("_atom_site.")
    auth_asym_id_idx = find_mmcif_tag_idx(atom_site_table, "auth_asym_id")
    lab_asym_id_idx = find_mmcif_tag_idx(atom_site_table, "label_asym_id")

    for atom in atom_site_table:
        auth = "".join([c for c in atom[auth_asym_id_idx] if c.isalnum()])
        label = "".join([c for c in atom[lab_asym_id_idx] if c.isalnum()])
        atom[auth_asym_id_idx] = auth
        atom[lab_asym_id_idx] = label

    with TemporaryDirectory() as tmp_dir_str:
        tmp_cif = "/".join([tmp_dir_str, inp_cif.name])
        in_doc.write_file(tmp_cif)

        try:
            proc = run(
                [
                    "chimerax",
                    "--nogui",
                    "--cmd",
                    addh_cmd_fmt.format(tmp_cif, out_cif),
                ],
                capture_output=True,
            )

        except (FileNotFoundError, NotADirectoryError):
            print("Chimerax failed")
            print(f"Standard output:\n{proc.stdout}\n\n")
            print(f"Standard error:{proc.stderr}\n\n")
            pass

        err = proc.stderr.decode()
        if err:
            print(f"Chimerax encountered error:\n{err}")
            raise ChimeraXException(str(inp_cif) + "\n" + err)

        return out_cif


"""
extract_chimX_H_atoms.py
"""


def is_polymer_atom(atom_row, group_idx=0):
    """Using group_PDB index = 0"""
    if atom_row[group_idx].strip().upper() == "ATOM":
        return True

    return False


def is_ligand_atom(atom_row, group_idx=0, comp_id_idx=5):
    """Using group_PDB index = 0; using label_comp_id index = 5"""
    if atom_row[group_idx].strip().upper() == "HETATM":
        if atom_row[comp_id_idx].strip().upper() != "HOH":
            return True

    return False


def is_HOH_atom(atom_row, group_idx=0, comp_id_idx=5):
    """Using group_PDB index = 0; using label_comp_id index = 5"""
    if atom_row[group_idx].strip().upper() == "HETATM":
        if atom_row[comp_id_idx].strip().upper() == "HOH":
            return True

    return False


def reorder_atoms_site(block, use_auth_comp_id=False):
    atom_site = block.find_mmcif_category("_atom_site.")

    group_idx = find_mmcif_tag_idx(atom_site, "group_PDB")
    model_idx = find_mmcif_tag_idx(atom_site, "pdbx_PDB_model_num")
    asym_idx = find_mmcif_tag_idx(atom_site, "auth_asym_id")
    seq_idx = find_mmcif_tag_idx(atom_site, "auth_seq_id")

    if use_auth_comp_id:
        comp_idx = find_mmcif_tag_idx(atom_site, "auth_comp_id")
    else:
        comp_idx = find_mmcif_tag_idx(atom_site, "label_comp_id")

    all_models = sorted(
        {atom[model_idx].strip() for atom in atom_site}, key=lambda x: int(x)
    )
    all_chains = sorted({atom[asym_idx].strip() for atom in atom_site})

    reordered_atoms = []
    for model in all_models:
        model_polymers = []
        model_ligands = []
        model_waters = []
        for chain in all_chains:
            chain_polymer = []
            chain_ligands = []
            chain_waters = []
            for atom in atom_site:
                if atom[model_idx].strip() == model and atom[asym_idx].strip() == chain:
                    if is_polymer_atom(atom, group_idx):
                        chain_polymer.append(list(atom))

                    elif is_ligand_atom(atom, group_idx, comp_idx):
                        chain_ligands.append(list(atom))

                    elif is_HOH_atom(atom, group_idx, comp_idx):
                        chain_waters.append(list(atom))

                    else:
                        raise RuntimeError(f"Atom out of place: {atom}")

            model_polymers.extend(sorted(chain_polymer, key=lambda x: int(x[seq_idx])))
            model_ligands.extend(sorted(chain_ligands, key=lambda x: int(x[seq_idx])))
            model_waters.extend(sorted(chain_waters, key=lambda x: int(x[seq_idx])))

        for atom in model_polymers:
            reordered_atoms.append(atom)

        for atom in model_ligands:
            reordered_atoms.append(atom)

        for atom in model_waters:
            reordered_atoms.append(atom)

    tags = [tag.split(".")[1] for tag in atom_site.tags]

    new_atom_site = block.init_mmcif_loop(atom_site.get_prefix(), tags)

    for atom in reordered_atoms:
        new_atom_site.add_row(atom)


def reset_atom_site_id_column(block):
    atom_site = block.find_mmcif_category("_atom_site.")

    for idx, row in enumerate(atom_site):
        row[find_mmcif_tag_idx(atom_site, "id")] = str(idx + 1)


def reorder_atoms_and_reset_id(block) -> None:
    reorder_atoms_site(block)

    reset_atom_site_id_column(block)


def fix_chimX_atom_site(cif_block, chimX_block):
    chimX_atom_site = chimX_block.find_mmcif_category("_atom_site.")

    group_idx = find_mmcif_tag_idx(chimX_atom_site, "group_PDB")
    id_idx = find_mmcif_tag_idx(chimX_atom_site, "id")
    type_idx = find_mmcif_tag_idx(chimX_atom_site, "type_symbol")
    lab_atom_id_idx = find_mmcif_tag_idx(chimX_atom_site, "label_atom_id")
    lab_alt_id_idx = find_mmcif_tag_idx(chimX_atom_site, "label_alt_id")
    lab_comp_id_idx = find_mmcif_tag_idx(chimX_atom_site, "label_comp_id")
    lab_asym_id_idx = find_mmcif_tag_idx(chimX_atom_site, "label_asym_id")
    lab_entity_id_idx = find_mmcif_tag_idx(chimX_atom_site, "label_entity_id")
    lab_seq_id_idx = find_mmcif_tag_idx(chimX_atom_site, "label_seq_id")
    ins_idx = find_mmcif_tag_idx(chimX_atom_site, "pdbx_PDB_ins_code")
    x_idx = find_mmcif_tag_idx(chimX_atom_site, "Cartn_x")
    y_idx = find_mmcif_tag_idx(chimX_atom_site, "Cartn_y")
    z_idx = find_mmcif_tag_idx(chimX_atom_site, "Cartn_z")
    occupancy_idx = find_mmcif_tag_idx(chimX_atom_site, "occupancy")
    b_factor_idx = find_mmcif_tag_idx(chimX_atom_site, "B_iso_or_equiv")
    # formal_charge_idx = None     # Added
    auth_seq_id_idx = find_mmcif_tag_idx(chimX_atom_site, "auth_seq_id")
    auth_comp_id_idx = find_mmcif_tag_idx(chimX_atom_site, "label_comp_id")  # Added
    auth_asym_id_idx = find_mmcif_tag_idx(chimX_atom_site, "auth_asym_id")
    auth_atom_id_idx = find_mmcif_tag_idx(chimX_atom_site, "label_atom_id")  # Added
    model_num_idx = find_mmcif_tag_idx(chimX_atom_site, "pdbx_PDB_model_num")

    cif_atom_site = cif_block.init_mmcif_loop(
        chimX_atom_site.get_prefix(), ATOM_SITE_TAGS
    )

    for atom in chimX_atom_site:
        cif_atom_site.add_row(
            AtomSiteAtom(
                atom[group_idx],
                atom[id_idx],
                atom[type_idx],
                atom[lab_atom_id_idx],
                atom[lab_alt_id_idx],
                atom[lab_comp_id_idx],
                atom[lab_asym_id_idx],
                atom[lab_entity_id_idx],
                atom[lab_seq_id_idx],
                atom[ins_idx],
                atom[x_idx],
                atom[y_idx],
                atom[z_idx],
                atom[occupancy_idx],
                atom[b_factor_idx],
                "?",
                atom[auth_seq_id_idx],
                atom[auth_comp_id_idx],
                atom[auth_asym_id_idx],
                atom[auth_atom_id_idx],
                atom[model_num_idx],
            )
        )


def extract_chimX_H_atoms(inp_cif, chimX_cif):
    inp_cif, chimX_cif = Path(inp_cif), Path(chimX_cif)

    cif_doc = gemmi.cif.read_file(str(inp_cif))
    chimX_doc = gemmi.cif.read_file(str(chimX_cif))

    cif_block = cif_doc.sole_block()
    chimX_block = chimX_doc.sole_block()

    reorder_atoms_site(chimX_block)

    fix_chimX_atom_site(cif_block, chimX_block)

    reset_atom_site_id_column(cif_block)

    out_cif = inp_cif.parent / (inp_cif.stem + "_H.cif")
    cif_doc.write_file(str(out_cif), gemmi.cif.Style.Pdbx)

    print(f"Extracted hydrogens... Wrote {out_cif}")

    return out_cif


"""
convert_Se_atoms.py
"""


def is_Se_atom(row, type_idx=2):
    if row[type_idx].strip() == "Se":
        return True

    return False


def is_Se_chem_comp(row, label_comp_idx, auth_comp_idx):
    Se_chem_comps = {"SEC", "MSE"}

    if row[label_comp_idx] in Se_chem_comps or row[auth_comp_idx] in Se_chem_comps:
        return True

    return False


def replace_Se_with_sulfur(block):
    atom_site = block.find_mmcif_category("_atom_site.")

    type_idx = find_mmcif_tag_idx(atom_site, "type_symbol")
    label_atom_idx = find_mmcif_tag_idx(atom_site, "label_atom_id")
    auth_atom_idx = find_mmcif_tag_idx(atom_site, "auth_atom_id")

    for atom in atom_site:
        if is_Se_atom(atom, type_idx):
            atom[type_idx] = "S"
            atom[label_atom_idx] = "S"
            atom[auth_atom_idx] = "S"


def replace_Se_residues(block):
    Se_chem_comp_map = {
        "MSE": "MET",
        "SEC": "CYS",
        "MET": "MET",
        "CYS": "CYS",
    }

    atom_site = block.find_mmcif_category("_atom_site.")

    group_idx = find_mmcif_tag_idx(atom_site, "group_PDB")
    label_comp_idx = find_mmcif_tag_idx(atom_site, "label_comp_id")
    auth_comp_idx = find_mmcif_tag_idx(atom_site, "auth_comp_id")

    for atom in atom_site:
        if is_Se_chem_comp(atom, label_comp_idx, auth_comp_idx):
            atom[group_idx] = "ATOM"
            atom[label_comp_idx] = Se_chem_comp_map.get(atom[label_comp_idx], "UNK")
            atom[auth_comp_idx] = Se_chem_comp_map.get(atom[auth_comp_idx], "UNK")


def correct_for_Se(inp_cif):
    inp_cif = Path(inp_cif)

    doc = gemmi.cif.read_file(str(inp_cif))

    block = doc.sole_block()

    replace_Se_with_sulfur(block)

    replace_Se_residues(block)

    reorder_atoms_site(block)

    reset_atom_site_id_column(block)

    out_cif = inp_cif.parent / (inp_cif.stem + "_corr_Se.cif")
    doc.write_file(str(out_cif), gemmi.cif.Style.Pdbx)

    print(f"Correcting for Selenium... Wrote {out_cif}")

    return out_cif


"""
correct_seq_id.py
"""


def correct_atom_site_seq_id(inp_cif: Path) -> Path:
    doc = gemmi.cif.read_file(str(inp_cif))
    block = doc.sole_block()

    model = gemmi.read_structure(str(inp_cif))[0]

    auth_residue_data = [
        Residue(chain.name, res.seqid.num, res.name, str(res.label_seq))
        for chain in model
        for res in chain.get_polymer()
    ]

    # check if we need to fix auth_seq_id column
    prev_asym_id = auth_residue_data[0].asym
    prev_seq_id = auth_residue_data[0].seq
    prev_label_seq_id = auth_residue_data[0].label_seq

    shift = 0
    found_diff = False

    new_seq_nums = {}
    for asym_id, seq_id, comp, label_seq_id in auth_residue_data:
        if asym_id != prev_asym_id:
            shift = 0

        if seq_id == prev_seq_id and label_seq_id != prev_label_seq_id:
            found_diff = True
            shift += 1

        new_seq_nums[asym_id, seq_id, comp, label_seq_id] = str(seq_id + shift)

        prev_asym_id = asym_id
        prev_seq_id = seq_id
        prev_label_seq_id = label_seq_id

        # print(f'{asym_id:<3s}{comp:<4s}{label_seq_id:<5s} old: {seq_id:5d} new: {seq_id + shift:5d}')

    updated_seq_ids = list()
    if found_diff:
        block.set_pair("_apchem.correct_auth_seq_id", "True")

        print(f"auth seq id need to be fixed for ", inp_cif.name)
        atom_site = block.find_mmcif_category("_atom_site.")

        group_pdb_idx = find_mmcif_tag_idx(atom_site, "group_PDB")
        auth_asym_idx = find_mmcif_tag_idx(atom_site, "auth_asym_id")
        auth_comp_idx = find_mmcif_tag_idx(atom_site, "auth_comp_id")
        auth_seq_idx = find_mmcif_tag_idx(atom_site, "auth_seq_id")
        label_seq_idx = find_mmcif_tag_idx(atom_site, "label_seq_id")

        for atom in atom_site:
            # this only corrects seq_id for polymers
            if atom[group_pdb_idx].strip() == "ATOM":
                if new_seq_nums.get(
                    (
                        atom[auth_asym_idx],
                        int(atom[auth_seq_idx]),
                        atom[auth_comp_idx],
                        atom[label_seq_idx],
                    ),
                    None,
                ):
                    updated_seq_id = new_seq_nums[
                        (
                            atom[auth_asym_idx],
                            int(atom[auth_seq_idx]),
                            atom[auth_comp_idx],
                            atom[label_seq_idx],
                        )
                    ]

                    updated_seq_ids.append(
                        (
                            atom[auth_asym_idx],
                            atom[auth_seq_idx],
                            atom[auth_comp_idx],
                            atom[label_seq_idx],
                            updated_seq_id,
                        )
                    )

                    atom[auth_seq_idx] = updated_seq_id

    else:
        block.set_pair("_apchem.correct_seq_id", "False")

        for res in auth_residue_data:
            updated_seq_ids.append(
                (res.asym, str(res.seq), res.comp, res.label_seq, str(res.seq))
            )

    loop = block.init_mmcif_loop(
        "_apchem_seq_id_corrections.",
        [
            "auth_asym_id",
            "orig_auth_seq_id",
            "auth_comp_id",
            "label_seq_id",
            "fixed_auth_seq_id",
        ],
    )

    updated_seq_ids = sorted(set(updated_seq_ids), key=lambda x: (x[0], int(x[-1])))
    for row in updated_seq_ids:
        loop.add_row(list(row))

    out_cif = inp_cif.parent / (inp_cif.stem + "_corr_seqID.cif")
    doc.write_file(str(out_cif), gemmi.cif.Style.Pdbx)

    print(f"Correct auth_seq_id... wrote: {out_cif}")

    return out_cif


"""
rm_atoms.py
"""


def elements_from_params(params):
    if params.exists():
        with params.open("rt") as f:
            params_data = json.load(f)
    else:
        raise FileNotFoundError(f"Parameter file passed does not exist: {params}")

    atom_data = params_data["atom"]["data"]

    params_elements = {param["key"][0].upper() for param in atom_data}

    return params_elements


def fetch_chem_comp_category(block):
    chem_comp = block.find_mmcif_category("_chem_comp.")

    if chem_comp:
        ChemComp = namedtuple("ChemComp", ("id", "type", "name", "formula", "mass"))

        return {
            ChemComp(comp[0], comp[1], comp[3], comp[5], comp[6]) for comp in chem_comp
        }

    else:
        raise RuntimeError(f"No chem_comp category in cif file {block.name}.cif")


def fetch_chem_comp_ids(block):
    chem_comps = fetch_chem_comp_category(block)

    return {comp.id for comp in chem_comps}


def fetch_non_polymer_chem_comp(block, exclude_water=True):
    chem_comps = fetch_chem_comp_category(block)

    non_polymers = {
        comp.id for comp in chem_comps if comp.type.strip() == "non-polymer"
    }

    if exclude_water:
        non_polymers.discard("HOH")

    return non_polymers


def fetch_exp_method(block):
    exp_method = block.find_values("_exptl.method")

    if exp_method:
        return exp_method[0]

    else:
        return None


def append_rm_atom_site_category(block, rm_atoms, tags, prefix="_rm_atom_site."):
    rm_atoms = [list(atom) for atom in rm_atoms]  # Prevents FPE

    table = block.find(prefix, tags)
    if table:
        loop = table.loop
    else:
        loop = block.init_loop(prefix, tags)

    for atom in rm_atoms:
        loop.add_row(atom)


def filter_atoms_to_rm_by_element(
    atom_site, allowed_elements=["C", "H", "N", "O", "S", "P", "F", "CL", "BR", "I"]
):
    """
    atom_site row list indice
    index 2 = _atom_site.type_symbol
    """

    type_idx = find_mmcif_tag_idx(atom_site, "type_symbol") or 2

    filtered_rows = [
        # TODO I should make the atom a list of str here rather than later in case of FPE
        (idx, atom)
        for idx, atom in enumerate(atom_site)
        if not atom[type_idx] in allowed_elements
    ]

    rm_idx, rm_rows = [], []
    if filtered_rows:
        rm_idx, rm_rows = zip(*filtered_rows)

    return rm_idx, rm_rows


def filter_atoms_to_rm_by_chem_comp_id(atom_site, chem_comp_ids, label_comp_id=False):
    """
    atom_site row list indice
    index 5  = _atom_site.label_comp_id
    index 17 = _atom_site.auth_comp_id
    """

    if isinstance(chem_comp_ids, (str, tuple, list, set)):
        if isinstance(chem_comp_ids, str):
            chem_comp_ids = [chem_comp_ids]
        if not isinstance(chem_comp_ids, set):
            chem_comp_ids = set(chem_comp_ids)
    else:
        raise TypeError(
            f'rm_residues parameter must be either a str of 1 comp_id" \
            "or a list/tuple/set of comp_id. rm_residues type: {type(chem_comp_ids)}'
        )

    label_idx = 5
    auth_idx = 17

    if label_comp_id:
        comp_idx = label_idx
    else:
        comp_idx = auth_idx

    filtered_rows = [
        (idx, atom)
        for idx, atom in enumerate(atom_site)
        if atom[comp_idx] in chem_comp_ids
    ]

    rm_idx, rm_rows = [], []
    if filtered_rows:
        rm_idx, rm_rows = zip(*filtered_rows)

    return rm_idx, rm_rows


def rm_atoms_from_atom_site(atom_site, rm_idx):
    for idx in rm_idx[::-1]:
        atom_site.remove_row(idx)


def rm_element_atoms(block, elements_allowed, prefix="_rm_atom_site."):
    atom_site = block.find_mmcif_category("_atom_site.")

    type_symbol = 2
    elem_idx = type_symbol

    rm_elem_idx, rm_elem_atoms = filter_atoms_to_rm_by_element(
        atom_site, elements_allowed
    )

    rm_elem_counts = Counter([row[elem_idx] for row in rm_elem_atoms])

    tags = [tag.split(".")[1] for tag in atom_site.tags]

    append_rm_atom_site_category(block, rm_elem_atoms, tags, prefix)

    rm_atoms_from_atom_site(
        block.find_mmcif_category("_atom_site."), rm_elem_idx
    )  # Prevents FPE
    # rm_atoms_from_atom_site(atom_site, rm_elem_idx) # May throw FPE

    return rm_elem_counts


def rm_chem_comp_atoms(block, chem_comp_ids, label=False):
    atom_site = block.find_mmcif_category("_atom_site.")

    label_idx = 5
    auth_idx = 17

    if label:
        comp_idx = label_idx
    else:
        comp_idx = auth_idx

    rm_water_idx, rm_water_atoms = filter_atoms_to_rm_by_chem_comp_id(
        atom_site, chem_comp_ids
    )

    water_count = Counter([row[comp_idx] for row in rm_water_atoms])

    append_rm_atom_site_category(block, rm_water_atoms, atom_site.tags)

    rm_atoms_from_atom_site(atom_site, rm_water_idx)

    return water_count


def rm_water_atoms(block, label=False):
    residue = {"HOH"}

    return rm_chem_comp_atoms(block, residue, label)


def keep_element_atoms(elements_allowed, add_elements):
    assert isinstance(add_elements, (tuple, set, list))

    elements_allowed = set(elements_allowed)

    for elem in add_elements:
        elements_allowed.add(elem)

    return elements_allowed


def keep_cation_atoms(elements_allowed):
    # TODO extend to additional elements
    cations = {"LI", "NA", "K", "MG", "CA", "CU", "FE", "NI", "CO", "MN", "ZN"}

    return keep_element_atoms(elements_allowed, cations)


def keep_halogen_atoms(elements_allowed):
    halogens = {"F", "CL", "BR", "I"}

    return keep_element_atoms(elements_allowed, halogens)


# Deprecated--there are other cation elements
def rm_cation_atoms(block, elements_allowed):
    elements_allowed = set(elements_allowed)

    cations = {"LI", "NA", "K", "MG", "CA", "CU", "FE", "NI", "CO", "MN", "ZN"}

    return rm_element_atoms(block, elements_allowed - cations)


def rm_halogen_atoms(block, elements_allowed):
    elements_allowed = set(elements_allowed)

    halogens = {"F", "CL", "BR", "I"}

    return rm_element_atoms(block, elements_allowed - halogens)


def add_back_chem_comp_to_atom_site(
    block, chem_comps, category="_rm_atom_site.", label=False
):
    other_atom_site = block.find_mmcif_category(category)

    atom_site = block.find_mmcif_category("_atom_site.")

    label_idx = find_mmcif_tag_idx(atom_site, "label_comp_id") or 5
    auth_idx = find_mmcif_tag_idx(atom_site, "auth_comp_id") or 17

    if label:
        comp_idx = label_idx
    else:
        comp_idx = auth_idx

    loop = atom_site.loop
    for row in other_atom_site:
        if row[comp_idx] in chem_comps:
            loop.add_row(list(row))


def add_back_water_atoms(block, category="_rm_atom_site.", label=False):
    add_back_chem_comp_to_atom_site(block, ["HOH"], category=category, label=label)


def rm_atoms_not_in_params(inp_cif, params, prefix="_rm_atom_site."):
    inp_cif = Path(inp_cif)
    params_dir = Path("/usr/local/share/parameters")

    params_file = params_dir / f"{params}.json"

    param_elements = elements_from_params(params_file)
    print(f"Elements found in params file: {param_elements}")

    doc = gemmi.cif.read_file(str(inp_cif))
    block = doc.sole_block()

    rm_elem_count = {}

    rm_elem_count = rm_element_atoms(block, param_elements, prefix)

    # Write out a file if atoms were removed
    rm_list = []
    if rm_elem_count:
        atoms_removed = Counter()
        if rm_elem_count:
            atoms_removed.update(rm_elem_count)

            rm_list.extend([elem.capitalize() for elem in rm_elem_count.keys()])

        else:
            print("No elements were removed.")

        print(f"Total Number of rows removed from atom_site: {atoms_removed}")

    else:
        print(f"Did not remove any atoms/elements.")

    rm_elem_str = "_rm_" + "_".join(rm_list)

    pdb_name, pdb_ext = inp_cif.stem, inp_cif.suffix

    out_cif = inp_cif.parent / (pdb_name + rm_elem_str + pdb_ext)

    doc.write_file(str(out_cif), gemmi.cif.Style.Pdbx)

    print(f"writing out: {out_cif}")

    return out_cif


"""
run_chargefw2.py
"""


def run_chargefw2_impl(inp_cif, method="eem", parameter_set="EEM_10_Cheminf_b3lyp_aim"):
    inp_cif = Path(inp_cif)

    out_cif = inp_cif.with_suffix(".fw2.cif")

    print(f"Running ChargeFW2....")
    try:
        protein = chargefw2_python.Molecules(str(inp_cif))

        charges = chargefw2_python.calculate_charges(protein, method, parameter_set)

        chargefw2_python.write_cif(protein, charges, str(inp_cif), str(inp_cif.parent))

        print(
            f"Calculated {method} charges with params {parameter_set}... Writing {out_cif}"
        )

    except Exception:
        pass

    return out_cif


def run_chargefw2(inp_cif, method="eem", parameter_set="EEM_10_Cheminf_b3lyp_aim"):
    gemmi_doc = gemmi.cif.read_file(str(inp_cif))
    run_as_seprate_model_docs(
        func=run_chargefw2_impl,
        gemmi_doc=gemmi_doc,
        method=method,
        parameter_set=parameter_set,
    )
    out_cif = Path(inp_cif).with_suffix(".fw2.cif")
    gemmi_doc.write_file(str(out_cif))
    return out_cif


"""
add_cations_back.py
"""


# TODO this needs to refactored so it uses the type_symbol rather than auth atom id
def add_back_cation_atoms(
    inp_cif, prefix="_cation_atom_site.", label=False, cleanup=False
):
    inp_cif = Path(inp_cif)

    doc = gemmi.cif.read_file(str(inp_cif))

    block = doc.sole_block()

    chem_comps = list(cation_formal_charges_radii.keys())

    other_atom_site = block.find_mmcif_category(prefix)

    atom_site = block.find_mmcif_category("_atom_site.")

    label_idx = find_mmcif_tag_idx(atom_site, "label_comp_id") or 5
    auth_idx = find_mmcif_tag_idx(atom_site, "auth_comp_id") or 17
    formal_charge_idx = find_mmcif_tag_idx(atom_site, "pdbx_formal_charge") or 15

    if label:
        comp_idx = label_idx
    else:
        comp_idx = auth_idx

    loop = atom_site.loop
    for atom in other_atom_site:
        if atom[comp_idx] in chem_comps:
            charge, radii = cation_formal_charges_radii[atom[comp_idx]]
            atom[formal_charge_idx] = str(int(float(charge)))

        else:
            charge, radii = "?", "?"

        atom = list(atom) + [charge, radii]

        loop.add_row(atom)

    if cleanup:
        other_atom_site.erase()

    stem = inp_cif.stem.split(".", 1)[0]

    out_cif = inp_cif.parent / (stem + "_cations" + ".fw2.cif").replace("__", "_")

    doc.write_file(str(out_cif), gemmi.cif.Style.Pdbx)

    print(f"Added back cations from {prefix}. Writing out: {out_cif}")

    return out_cif


def run_freesasa_impl(inp_cif):
    inp_cif = Path(inp_cif)

    if inp_cif.suffix != ".cif":
        raise ValueError(f"Input to freesasa cannot be a cif file: {inp_cif}")

    out_cif = inp_cif.with_suffix(".sasa.cif")

    # Old Implementation
    # Delete once the write_cif bindings are implemented
    freesasa_flags = [
        "--no-warnings",
        "--shrake-rupley",
        "--hydrogen",
        " --hetatm",
        "--cif",
        "--format=cif",
    ]

    run_cmd = "freesasa " + " ".join(freesasa_flags + [str(inp_cif)])
    print(f"Running FreeSASA: {run_cmd}")
    with out_cif.open("wt") as f:
        run(run_cmd.split(), stdout=f)

    if out_cif.stat().st_size != 0:
        print(f"Writing out: {out_cif}")

    return out_cif


def run_freesasa(inp_cif):
    inp_cif = Path(inp_cif)

    if inp_cif.suffix != ".cif":
        raise ValueError(f"Input to freesasa cannot be a cif file: {inp_cif}")

    out_cif = inp_cif.with_suffix(".sasa.cif")

    gemmi_doc = gemmi.cif.read_file(str(inp_cif))
    run_as_seprate_model_docs(func=run_freesasa_impl, gemmi_doc=gemmi_doc)
    gemmi_doc.write_file(str(out_cif))

    return out_cif


def process_pdb(pdb: str, cif_dir: Path, pdb_redo=True):
    if len(pdb) != 4 and not pdb.lower().endswith("cif"):
        raise ValueError("input pdb must be an RCSB pdb-id or a cif filename")

    cif_dir.mkdir(0o744, parents=True, exist_ok=True)

    if pdb.lower().endswith("cif"):
        cif_file = cif_dir / pdb

        if cif_file.exists():
            print(f"Using cached file ({cif_file}) for pipeline...")
            return cif_file

        else:
            raise FileNotFoundError(cif_file)

    else:
        return fetch_cif_from_api(pdb.lower(), cif_dir, pdb_redo)


def run_chimX_fw2_freesasa_pipeline(
    init_cif,
    fw2_method="eem",
    fw2_params="EEM_10_Cheminf_b3lyp_aim",
    prefix="_rm_atom_site.",
    remove_intermediate=True,
):
    rm_cations = False

    files_to_remove = []

    try:
        chimX_cif = run_chimeraX(init_cif)
        files_to_remove.append(chimX_cif)

        if not chimX_cif.exists() or chimX_cif.stat().st_size == 0:
            print(f"ChimeraX failed for {init_cif}... Removing cations and trying again...")

            rm_cation_cif = rm_atoms_not_in_params(init_cif, fw2_params, prefix)

            rm_cations = True

            chimX_cif = run_chimeraX(rm_cation_cif)
            files_to_remove.append(chimX_cif)

            if not chimX_cif.exists() or chimX_cif.stat().st_size == 0:
                print(f"ChimeraX Failed: {init_cif.stem}")
                raise ChimeraXException(init_cif.name)

        if rm_cations:
            H_cif = extract_chimX_H_atoms(rm_cation_cif, chimX_cif)
            files_to_remove.append(H_cif)

            Se_cif = correct_for_Se(H_cif)
            files_to_remove.append(Se_cif)

            charge_cif = run_chargefw2(Se_cif, fw2_method, fw2_params)
            files_to_remove.append(charge_cif)

        else:
            H_cif = extract_chimX_H_atoms(init_cif, chimX_cif)
            files_to_remove.append(H_cif)

            Se_cif = correct_for_Se(H_cif)
            files_to_remove.append(Se_cif)

            seqid_cif = correct_atom_site_seq_id(Se_cif)
            files_to_remove.append(seqid_cif)

            rm_cation_cif = rm_atoms_not_in_params(seqid_cif, fw2_params, prefix)
            files_to_remove.append(rm_cation_cif)

            charge_cif = run_chargefw2(rm_cation_cif, fw2_method, fw2_params)
            files_to_remove.append(charge_cif)

        if not charge_cif.exists() or charge_cif.stat().st_size == 0:
            print(f"ChargeFW2 Failed: {init_cif.stem}")
            raise ChargeFW2Exception(init_cif.name)

        cation_cif = add_back_cation_atoms(charge_cif, prefix)
        files_to_remove.append(cation_cif)

        sasa_cif = run_freesasa(cation_cif)
        files_to_remove.append(sasa_cif)

        if not sasa_cif.exists() or sasa_cif.stat().st_size == 0:
            print(f"Freesasa Failed: {init_cif.stem}")
            sasa_cif.unlink(missing_ok=True)
            raise FreeSASAException(init_cif.name)

        final_cif = rename_cif_by_split(sasa_cif, delimiter="_")

        print(f"Final cif: {final_cif.resolve()}")
    finally:
        if remove_intermediate:
            for f in files_to_remove:
                if f.exists():
                    print(f"Deleting intermediate file {f}")
                    f.unlink()

    return final_cif
