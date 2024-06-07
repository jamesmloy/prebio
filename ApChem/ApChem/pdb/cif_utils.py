"""
Usage: python generate_modCIF.py --help

dependencies:   
    pdb2pqr v2.1 (only this version)
    freesasa v2.2.0a1 (only this version) (pip)
    gemmi (pip)
"""

import gemmi
import freesasa
import subprocess
import numpy as np
import requests
import pathlib
import argparse
import os


def gen_pdb_from_api(pdb, folder, redo=False, pdb_format=False):
    print("Fetching PDB from api...")

    def fetch_from_pdb_redo(pdb, folder, pdb_format=False):
        print("Fetching from PDB-REDO server...")

        pdb_filename = f"{pdb}_final_tot.pdb"

        if pdb_format:
            pdb_filename = f"{pdb}_final_tot.pdb"
        else:
            pdb_filename = f"{pdb}_final.cif"

        pdb_file = folder / pdb_filename

        respond = requests.get(f"https://pdb-redo.eu/db/{pdb}/{pdb_filename}")
        if respond.ok:
            with pdb_file.open("w") as f:
                f.write(respond.text)
            return pdb_file
        else:
            return False

    def fetch_from_rcsb(pdb, folder, pdb_format=False):
        print("Fetching from RCSB...")

        pdb_filename = ""
        if pdb_format:
            pdb_filename = f"{pdb}.pdb"
        else:
            pdb_filename = f"{pdb}.cif"

        pdb_file = folder / pdb_filename

        respond = requests.get(f"https://files.rcsb.org/download/{pdb_filename}")
        if respond.ok:
            with pdb_file.open("w") as f:
                f.write(respond.text)
            return pdb_file
        else:
            raise FileNotFoundError(
                f"Unable to retrieve pdb file from rcsb: {pdb_filename}"
            )

    assert len(pdb.stem) == 4
    folder = pathlib.Path(folder)

    pdb_file = ""
    if redo:
        pdb_file = fetch_from_pdb_redo(pdb, folder, pdb_format)
        if not pdb_file:
            pdb_file = fetch_from_rcsb(pdb, folder, pdb_format)
    else:
        pdb_file = fetch_from_rcsb(pdb, folder, pdb_format)

    return pdb_file


def setup_pdb_files(pdb_inp, folder="./data", redo=False, pdb_format=False):

    pdb_inp = pathlib.Path(pdb_inp)
    folder = pathlib.Path(folder)
    folder.mkdir(0o744, parents=True, exist_ok=True)

    if pdb_inp.suffix == "":
        pdb_inp = gen_pdb_from_api(pdb_inp, folder, redo=redo, pdb_format=pdb_format)
    elif pdb_inp.suffix == ".pdb" or pdb_inp.suffix == ".cif":
        new_pdb_inp = folder / pdb_inp.name
        new_pdb_inp.write_text(pdb_inp.read_text())
        pdb_inp = new_pdb_inp
    else:
        raise ValueError("Input protein must be a pdb (*.pdb) or cif (*.cif) file.")

    print(f"Input File: {pdb_inp}")

    pdb_files = {
        "input": pdb_inp.resolve(),
        "code": pdb_inp.stem[:4],
        "inp.cif": pdb_inp.with_suffix(".cif").resolve(),
        "pqr.pdb": pdb_inp.with_suffix(".pqr.pdb").resolve(),
        "pqr.cif": pdb_inp.with_suffix(".pqr.cif").resolve(),
        "pqr.log": pdb_inp.with_suffix(".pqr.log").resolve(),
        "sasa.pdb": pdb_inp.with_suffix(".sasa.pdb").resolve(),
        "sasa.cif": pdb_inp.with_suffix(".sasa.cif").resolve(),
        "output": (pdb_inp.parent / (pdb_inp.stem + "_mod.cif")).resolve(),
    }

    pqr_inp = pdb_files["pqr.pdb"]
    sasa_inp = pdb_files["sasa.pdb"]

    run_pdb2pqr(pdb_inp, pqr_inp)
    if (not pqr_inp.exists()) or (pqr_inp.stat().st_size == 0):
        raise RuntimeError(
            f"PQR file ({pqr_inp}) was not correctly generated. Cannot generate SASA file"
        )
    else:
        print(f"Finished generating pqr  file: {pqr_inp}")

    run_freesasa(pqr_inp, sasa_inp)
    if (not sasa_inp.exists()) or (sasa_inp.stat().st_size == 0):
        raise RuntimeError(f"SASA file ({sasa_inp}) was not correctly generated.")
    else:
        print(f"Finished generating sasa file: {sasa_inp}")

    return pdb_files


def run_pdb2pqr(inp_file, out_pqr):

    # Meant to run in docker
    pdb2pqr = pathlib.Path("/apchem/dependencies/pdb2pqr-2.1/pdb2pqr.py")

    pdb2pqr_flags = [
        "--ff=PARSE",
        "--chain",
        "--include-header",
    ]

    inp_file = pathlib.Path(inp_file).resolve()
    out_pqr = pathlib.Path(out_pqr).resolve()

    if inp_file.suffix == ".cif" and not "final" in inp_file.stem:
        inp_file = inp_file.with_suffix("")

    run_cmd = "python2  " + " ".join(
        [str(pdb2pqr)] + pdb2pqr_flags + [str(inp_file), str(out_pqr)]
    )
    print(f"Generating pqr: {run_cmd}")

    cwd = pathlib.Path.cwd()
    os.chdir(pathlib.Path(__file__).parent)

    if pdb2pqr.exists():
        p = subprocess.Popen(
            run_cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        p.communicate()
        os.chdir(cwd)
    else:
        raise FileNotFoundError(pdb2pqr)


def run_freesasa(in_pdb, out_sasa):
    in_pdb = pathlib.Path(in_pdb)
    out_sasa = pathlib.Path(out_sasa)

    if in_pdb.suffix == ".cif":
        raise ValueError(f"Input to freesasa cannot be a cif file: {in_pdb}")

    # Old Implementation prior to adding write_pdb python bindings
    # Delete once the write_pdb bindings are available on pypi
    freesasa_flags = ['--no-warnings',
                       ' --shrake-rupley',
                       ' --depth=atom',
                       '--hydrogen',
                       ' --hetatm',
                       '--format=pdb ']

    run_cmd = 'freesasa ' + ' '.join(freesasa_flags + [str(in_pdb)]) 
    print(f'Generating sasa: {run_cmd}')
    with out_sasa.open('wb') as sasa_out:
        p = subprocess.Popen(run_cmd.split(), stdout=sasa_out, stderr=subprocess.PIPE)
        p.communicate()

    # New Implementation but for some reason freesasa==2.2.0a1 cannot 
    # find Result.write_pdb() function for some reason.
    # options = {
    #     "hetatm": True,
    #     "hydrogen": True,
    #     "skip-unknown": False,
    # }

    # params = freesasa.Parameters()
    # params.setAlgorithm(freesasa.ShrakeRupley)

    # structure = freesasa.Structure(str(in_pdb), options=options)
    # result = freesasa.calc(structure, params)

    # result.write_pdb(str(out_sasa))

    # print(f"Generating sasa: {out_sasa}")


def write_pdb2cif(inp_pdb, out_cif):
    pdb_struct = gemmi.read_pdb(str(inp_pdb))
    pdb_struct.setup_entites()
    pdb_struct.assign_label_seq_id()
    cif_doc = pdb_struct.make_mmcif_document()
    cif_doc.write_file(str(out_cif))


def open_cif2doc(inp_cif):
    return gemmi.cif.read_file(str(inp_cif))


def update_category_from_rcsb_cif(doc, pdb_id, categories):

    if len(pdb_id) != 4:
        raise ValueError(
            f"PDB id ({pdb_id}) provided is not 4 chars but {len(pdb_id)} chars."
        )

    if isinstance(categories, str):
        categories = [categories]
    if not isinstance(categories, (list, tuple)):
        raise TypeError("categories parameter must be set to list/tuple.")

    respond = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif")

    if respond.ok:
        rcsb_block = gemmi.cif.read_string(respond.text).sole_block()
    else:
        raise RuntimeError(f"Failed ot download cif file from RCSB: {pdb_id}")

    block = doc.sole_block()

    for category in categories:
        cat_data = rcsb_block.get_mmcif_category(category)
        block.set_mmcif_category(category, cat_data)

    print(f"Finished updating categories: {categories} from the RCSB")

    return doc


def write_ciffile(doc, filename):
    doc.check_for_missing_values()
    doc.check_for_duplicates()
    doc.write_file(str(filename))


def merge_block_cifDoc(origDoc, newDoc, loopName="_atom_site."):
    block = origDoc.sole_block()
    newBlock = newDoc.sole_block()

    # Column tags for the atom_site loop
    newTags = list(newBlock.get_mmcif_category(loopName).keys())

    replace_val_lst = [
        str(list(k)[0]) for tag in newTags for k in newBlock.find(loopName, [tag])
    ]
    replace_val_lst = np.reshape(
        replace_val_lst, (len(newTags), int(len(replace_val_lst) / len(newTags)))
    )

    loop = block.init_loop(loopName, newTags)
    loop.set_all_values(replace_val_lst)

    return origDoc


def rename_column_cifDoc(doc, loopName, oldCol, newCol):
    block = doc.sole_block()
    tags = list(block.get_mmcif_category(loopName).keys())
    new_tags = [k.replace(oldCol, newCol) for k in tags]
    val_lst = [str(list(k)[0]) for tag in tags for k in block.find(loopName, [tag])]
    val_lst = np.reshape(val_lst, (len(tags), int(len(val_lst) / len(tags))))
    loop = block.init_loop(loopName, new_tags)
    loop.set_all_values(val_lst)

    return doc


def extract_column_cifDoc(doc, loopName, colName):
    block = doc.sole_block()
    values = [str(list(k)[0]) for k in block.find(loopName, [colName])]

    return values


def add_column_cifDoc(doc, newFieldName, newFieldValues, loopName="_atom_site."):
    block = doc.sole_block()
    tags = list(block.get_mmcif_category(loopName).keys()) + [newFieldName]
    val_lst = [
        str(list(k)[0]) for tag in tags for k in block.find(loopName, [tag])
    ] + newFieldValues
    val_lst = np.reshape(val_lst, (len(tags), int(len(val_lst) / len(tags))))
    loop = block.init_loop(loopName, tags)
    loop.set_all_values(val_lst)

    return doc


def clean_files(pdb_files):

    file_types = ["input", "pqr.pdb", "pqr.cif", "pqr.log", "sasa.cif", "sasa.pdb"]

    for f_type, f_path in pdb_files.items():
        if f_type in file_types and f_path.exists():
            f_path.unlink()


def generate_modCIF_from_pdb(
    protein,
    data_dir="./data",
    redo=False,
    pdb_format=True,
    keep_temp_files=False,
    ignore_cache=False,
):

    protein = pathlib.Path(protein)
    data_dir = pathlib.Path(data_dir)

    cached_file = pathlib.Path("")
    if not ignore_cache:
        if protein.suffix == "" and len(protein.stem) == 4:
            cached_rcsb_file = data_dir / (protein.stem + "_mod.cif")
            cached_redo_file = data_dir / (protein.stem + "_final_tot_mod.cif")
        elif protein.suffix in [".pdb", ".cif"] and (
            not protein.name.endswith("_mod.cif")
        ):
            cached_rcsb_file = data_dir / (protein.stem[:4] + "_mod.cif")
            cached_redo_file = data_dir / (protein.stem[:4] + "_final_tot_mod.cif")
        elif protein.name.endswith("_mod.cif"):
            cached_rcsb_file = data_dir / protein.name

        print(f"Checking for cached file: {cached_redo_file}")
        if cached_redo_file.exists() and cached_redo_file.stat().st_size != 0:
            print(f"returning cached file: {cached_redo_file}")
            return cached_redo_file

        print(f"Checking for cached file: {cached_rcsb_file}")
        if cached_rcsb_file.exists() and cached_rcsb_file.stat().st_size != 0:
            print(f"returning cached file: {cached_rcsb_file}")
            return cached_rcsb_file

    pdb_files = setup_pdb_files(protein, data_dir, redo, pdb_format)

    print(f'Generating modified cif  file: {pdb_files["output"]}')

    if pdb_files["input"].suffix == ".pdb":
        write_pdb2cif(pdb_files["input"], pdb_files["inp.cif"])

    write_pdb2cif(pdb_files["pqr.pdb"], pdb_files["pqr.cif"])
    write_pdb2cif(pdb_files["sasa.pdb"], pdb_files["sasa.cif"])

    orig_doc = open_cif2doc(pdb_files["inp.cif"])
    pqr_doc = open_cif2doc(pdb_files["pqr.cif"])
    sasa_doc = open_cif2doc(pdb_files["sasa.cif"])

    orig_doc = merge_block_cifDoc(orig_doc, pqr_doc, loopName="_atom_site.")

    pqr_vals = extract_column_cifDoc(
        pqr_doc, loopName="_atom_site.", colName="occupancy"
    )
    sasa_vals = extract_column_cifDoc(
        sasa_doc, loopName="_atom_site.", colName="B_iso_or_equiv"
    )

    add_column_cifDoc(
        orig_doc,
        newFieldName="partial_charge",
        newFieldValues=pqr_vals,
        loopName="_atom_site.",
    )
    add_column_cifDoc(
        orig_doc,
        newFieldName="solvent_accessibility",
        newFieldValues=sasa_vals,
        loopName="_atom_site.",
    )
    # This is needed if using pdb2pqr 2.1 not in 3.x
    out_doc = update_category_from_rcsb_cif(
        orig_doc, pdb_files["code"], ["_entity.", "_entity_poly."]
    )

    write_ciffile(orig_doc, pdb_files["output"])
    print(f'Finished writing modified cif file: {pdb_files["output"]}')

    if not keep_temp_files:
        clean_files(pdb_files)

    return pdb_files["output"]


def generate_modCIF_from_pdb_dataset(
    pdb_file,
    protein_dir,
    redo=False,
    pdb_format=True,
    keep_temp_files=False,
    ignore_cache=False,
):

    import pymp

    with pdb_file.open("rt") as f:
        pdb_codes = [pdb for pdb in f.read().split("\n")]

    with pymp.Parallel() as p:
        for i in p.xrange(len(pdb_codes)):
            try:
                generate_modCIF_from_pdb(
                    pdb_codes[i],
                    protein_dir,
                    redo,
                    pdb_format,
                    keep_temp_files,
                    ignore_cache,
                )
            except Exception as e:
                print(f"Error with {pdb_codes[i]}: {e}")
                continue
