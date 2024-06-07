#!/usr/bin/python3
"""danny305 8/12/21"""

from pathlib import Path
from typing import Callable, Union

import gemmi

from ApChem.cif import find_mmcif_tag_idx, reorder_atoms_and_reset_id


def fix_Se_bug(inp_cif: Union[Path, str]):

    inp_cif = Path(inp_cif)

    doc = gemmi.cif.read_file(str(inp_cif))

    block = doc.sole_block()

    atom_site = block.find_mmcif_category("_atom_site.")

    group_idx = find_mmcif_tag_idx(atom_site, "group_PDB")
    auth_comp_idx = find_mmcif_tag_idx(atom_site, "auth_comp_id")

    for atom in atom_site:
        if atom[group_idx] == "HETATM" and (
            atom[auth_comp_idx] == "MET" or atom[auth_comp_idx] == "CYS"
        ):
            atom[group_idx] = "ATOM"

    reorder_atoms_and_reset_id(block)

    doc.write_file(str(inp_cif), gemmi.cif.Style.Pdbx)

    return inp_cif


def fix_asym_id_bug(inp_cif: Union[Path, str]):

    inp_cif = Path(inp_cif)

    doc = gemmi.cif.read_file(str(inp_cif))

    block = doc.sole_block()

    reorder_atoms_and_reset_id(block)

    doc.write_file(str(inp_cif), gemmi.cif.Style.Pdbx)

    return inp_cif


def patch_cif_file(cif: Path, func: Callable):

    assert isinstance(cif, Path) and cif.exists()
    assert callable(func)

    func(cif)

    if cif.stat().st_size == 0:
        print(f"Patch Failed: {cif.stem}")
        raise RuntimeError(cif.name)

    print(f"Patched: {cif.name}")


def patch_cif_dir(cif_dir: Path, func: Callable, glob: str = "*.cif"):

    import pymp

    assert isinstance(cif_dir, Path)
    assert cif_dir.exists() and cif_dir.is_dir()
    assert callable(func)

    cif_files = list(set(cif_dir.glob(glob)))

    total_cif_files = len(cif_files)

    failed_proteins = pymp.shared.list()

    with pymp.Parallel() as p:
        for idx in p.xrange(total_cif_files):
            cif = cif_files[idx]

            try:
                patch_cif_file(cif, func)

            except Exception:
                failed_proteins.append(cif.stem)
                p.print(
                    f"Added failed protein {cif.stem}. Total count: {len(failed_proteins)}"
                )

            finally:
                p.print(
                    f"Thread: {p.thread_num} Finished protein {cif.stem}: {idx+1} out of {total_cif_files}"
                )

    log_failed_proteins = cif_dir / f"failed_{func.__name__}.txt"

    if len(failed_proteins) > 0:
        with log_failed_proteins.open("wt") as f:
            f.write(f"#Dataset: {cif_dir}\n")
            f.write("\n".join(failed_proteins))
