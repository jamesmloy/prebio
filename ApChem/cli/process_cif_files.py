from pathlib import Path
import shutil
from typing import List, Iterable
from ApChem.cif import run_chimX_fw2_freesasa_pipeline
from argparse import ArgumentParser
import os
import traceback


def get_all_files_with_extensions(
    root_path: Path, extensions: Iterable[str]
) -> List[Path]:
    the_files = []
    folders_stack = [root_path]

    while folders_stack:
        folder = folders_stack.pop()
        for f in folder.iterdir():
            if f.is_file() and (f.suffix in extensions):
                the_files.append(f)
            elif f.is_dir():
                folders_stack.append(f)

    return the_files


def process_cif(
    pdb_file: Path,
    out_dir: Path,
    scratch_dir: Path,
    dry_run: bool,
    keep_intermediate=False,
):
    remove_intermediate = not keep_intermediate
    files_to_delete = []
    try:
        final_cif = out_dir / (pdb_file.stem.lower() + ".fw2.sasa.cif")
        if not final_cif.exists():
            if not dry_run:
                scratch_pdb = Path(shutil.copy(pdb_file, scratch_dir))
                files_to_delete.append(scratch_pdb)
                final_cif = run_chimX_fw2_freesasa_pipeline(
                    init_cif=scratch_pdb, remove_intermediate=remove_intermediate
                )
                files_to_delete.append(final_cif)
                shutil.copy(final_cif, out_dir)
                if remove_intermediate:
                    _ = [f.unlink() for f in files_to_delete if f.exists()]

    except Exception as e:
        print(f"Encountered exception running pipeline: {e}\n{traceback.format_exc()}")
        if remove_intermediate:
            _ = [f.unlink() for f in files_to_delete if f.exists()]
        raise


def process_cif_files(
    cif_dir: Path, out_dir: Path, dry_run: bool, keep_intermediate=False
):
    all_pdbs = get_all_files_with_extensions(cif_dir, [".cif", ".pdb"])

    scratch_folder = out_dir / "scratch"
    scratch_folder.mkdir(exist_ok=True, parents=True)

    import pymp

    failed_proteins = pymp.shared.list()
    total_proteins = len(all_pdbs)
    processed_proteins = pymp.shared.list()

    with pymp.Parallel() as p:
        for idx in p.xrange(total_proteins):
            try:
                process_cif(
                    all_pdbs[idx], out_dir, scratch_folder, dry_run, keep_intermediate
                )
                processed_proteins.append(all_pdbs[idx])
            except Exception as e:
                print(f"{all_pdbs[idx].name} failed with exception: {e}")
                failed_proteins.append(all_pdbs[idx])

    print(f"Processed {len(processed_proteins)} proteins")
    failed_path = out_dir / "failed.txt"
    failed_path.unlink(missing_ok=True)
    with failed_path.open("w") as f:
        _ = [f.write(str(pdb) + "\n") for pdb in failed_proteins]

    print("\n==========\nFINISHED\n==========\n")


def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--in-folder",
        type=Path,
        required=True,
        help="The folder containing all the cif/pdb files (could be nested)",
    )
    parser.add_argument(
        "--out-folder",
        type=Path,
        required=True,
        help="The target folder where all the new cif files will go",
    )
    parser.add_argument(
        "--num-proc", type=int, help="Number of processes to use", default=10
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--keep-intermediate",
        action="store_true",
        help="Keep all intermediate files created (increases storage ~5x)",
    )
    args = parser.parse_args()
    if not args.in_folder.exists():
        raise RuntimeError(f"Folder {args.in_folder} does not exist")
    if not args.out_folder.exists():
        raise RuntimeError(f"Folder {args.out_folder} does not exist")
    return args


if __name__ == "__main__":
    args = get_args()

    os.environ["PYMP_NUM_THREADS"] = str(args.num_proc)
    os.environ["OMP_NUM_THREADS"] = "4"

    process_cif_files(
        args.in_folder, args.out_folder, args.dry_run, args.keep_intermediate
    )
