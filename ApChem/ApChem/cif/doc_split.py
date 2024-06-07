import gemmi
from collections import defaultdict
from typing import Any
from tempfile import TemporaryDirectory
from pathlib import Path


STD_ATOM_SITE_TAGS = [
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


def find_mmcif_tag_idx(mmcif_category: gemmi.cif.Table, tag: str) -> int:
    tags = [tag.split(".")[-1] for tag in mmcif_category.tags]

    try:
        return tags.index(tag)
    except ValueError:
        return -1


def init_atom_site_table(doc: gemmi.cif.Document, tags: list[str] = STD_ATOM_SITE_TAGS):
    return doc.sole_block().init_mmcif_loop("_atom_site.", tags)


def get_atom_site_table(doc: gemmi.cif.Document) -> gemmi.cif.Table:
    return doc.sole_block().find_mmcif_category("_atom_site.")


def get_atom_site_tags(doc: gemmi.cif.Document) -> list[str]:
    return [
        t.split(".")[-1]
        for t in doc.sole_block().find_mmcif_category("_atom_site.").tags
    ]


def group_atoms_by_model(doc: gemmi.cif.Document) -> dict[int, list[list[str]]]:
    atom_site_table = get_atom_site_table(doc)
    model_no_idx = find_mmcif_tag_idx(atom_site_table, "pdbx_PDB_model_num")
    model_to_atom = defaultdict(list)

    for atom_row in atom_site_table:
        model_no = atom_row[model_no_idx]
        model_to_atom[model_no].append([t for t in atom_row])

    return model_to_atom


def split_doc_by_models(doc: gemmi.cif.Document) -> dict[int, gemmi.cif.Document]:
    per_model_docs: dict[int, gemmi.cif.Document] = dict()

    atoms_by_model = group_atoms_by_model(doc)
    og_tags = get_atom_site_tags(doc)
    model_no_idx = find_mmcif_tag_idx(get_atom_site_table(doc), "pdbx_PDB_model_num")

    for model_no, atom_list in atoms_by_model.items():
        doc_for_model = gemmi.cif.read_file(doc.source)

        atom_site_loop = init_atom_site_table(doc_for_model, og_tags)

        for atom in atom_list:
            # we have to set the model number to 1 or else
            # chargefw2 freaks out

            atom[model_no_idx] = "1"
            atom_site_loop.add_row(atom)
        per_model_docs[model_no] = doc_for_model

    return per_model_docs


def combine_atom_site_tables(
    docs: dict[int, gemmi.cif.Document],
    doc_to_replace: gemmi.cif.Document | None = None,
) -> gemmi.cif.Document:

    if len(docs):
        # if a template doc was not provided, just use the first one
        if not doc_to_replace:
            doc_to_replace = gemmi.cif.read_file(docs[0].source)

        for doc in docs.values():
            new_atom_site_table = init_atom_site_table(
                doc_to_replace, get_atom_site_tags(doc)
            )
            break

        model_no_idx = find_mmcif_tag_idx(new_atom_site_table, "pdbx_PDB_model_num")

        for model_no, doc in docs.items():
            atoms_to_add = get_atom_site_table(doc)
            for atom in atoms_to_add:
                atom[model_no_idx] = model_no
                new_atom_site_table.add_row(atom)

        return doc_to_replace

    else:
        raise ValueError("Must pass a non-empty list of cif documents")


def run_as_seprate_model_docs(
    func: Any,
    gemmi_doc: gemmi.cif.Document,
    operate_on_file=True,
    logger=print,
    **other_func_kwargs
):
    docs_by_model = split_doc_by_models(gemmi_doc)

    # if there's only 1 model, dont bother splitting
    if len(docs_by_model) == 1:
        logger("Doc only contains 1 model-- skipping split")

        if operate_on_file:
            # if we're operating on a file, then we
            # need to put the contents of the temp file
            # in the original doc
            out_fname = func(Path(gemmi_doc.source), **other_func_kwargs)
            doc_map = dict()
            for k in docs_by_model.keys():
                doc_map[k] = gemmi.cif.read_file(str(out_fname))
            out_fname.unlink(missing_ok=True)
            combine_atom_site_tables(doc_map, gemmi_doc)
        else:
            # we're operating on the doc, then
            # just call the function
            func(gemmi_doc, **other_func_kwargs)

    else:
        doc_map = dict()
        logger(
            f"Doc contains {len(docs_by_model)}-- "
            f"splitting into {len(docs_by_model)} separate docs"
        )
        with TemporaryDirectory() as temp_dir_str:
            temp_dir = Path(temp_dir_str)
            for model_no, doc in docs_by_model.items():

                if operate_on_file:
                    in_fname = temp_dir / f"model-{str(model_no)}.cif"
                    doc.write_file(str(in_fname))
                    out_fname = func(in_fname, **other_func_kwargs)
                    doc_map[model_no] = gemmi.cif.read_file(str(out_fname))
                else:
                    func(doc, **other_func_kwargs)
                    doc_map[model_no] = doc

        combine_atom_site_tables(doc_map, gemmi_doc)
