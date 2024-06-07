from pathlib import Path


from ApChem.exceptions import FreeSASAException
from ApChem.cif import process_pdb, run_chimX_fw2_freesasa_pipeline, run_freesasa



def generate_cif_dataset(
    pdbs: Path,
    cif_dir: Path,
    pdb_redo=True,
    fw2_method="eem",
    fw2_params="EEM_10_Cheminf_b3lyp_aim",
    prefix="_rm_atom_site.",
    use_cache=True,
):

    import pymp

    assert pdbs.is_file()
    assert cif_dir.is_dir()

    with pdbs.open("rt") as f:
        pdb_list = [line.strip() for line in f.readlines() if not "#" in line]

    total_proteins = len(pdb_list)

    failed_proteins = pymp.shared.list()

    with pymp.Parallel() as p:
        for idx in p.xrange(total_proteins):
            pdb = pdb_list[idx]

            cached_fw2_sasa_cif = [cif for cif in cif_dir.glob(f"{pdb[:4]}*.fw2.sasa.cif")]

            if use_cache and len(cached_fw2_sasa_cif) == 1:
                p.print(
                    f"Thread: {p.thread_num} Cached protein {pdb}: {idx+1} out of {total_proteins}"
                )
                continue

            try:
                init_cif = process_pdb(pdb, cif_dir, pdb_redo)

                if not init_cif.exists() or init_cif.stat().st_size == 0:
                    print(f"API Failed: {init_cif.stem}")
                    raise FileNotFoundError(init_cif)

                run_chimX_fw2_freesasa_pipeline(
                    init_cif, fw2_method, fw2_params, prefix
                )

            except Exception:
                failed_proteins.append(pdb)
                p.print(
                    f"Added failed protein {pdb}. Total count: {len(failed_proteins)}"
                )

            finally:
                p.print(
                    f"Thread: {p.thread_num:3d} Finished protein {pdb}: {idx+1} out of {total_proteins}"
                )

    log_failed_proteins = cif_dir / "failed_pdbs.txt"

    if len(failed_proteins) > 0:
        with log_failed_proteins.open("wt") as f:
            f.write(f"#Dataset: {pdbs}\n")
            f.write("\n".join(failed_proteins))


def generate_freesasa_dataset(pdbs: Path, cif_dir: Path, pdb_redo=True):

    import pymp

    assert pdbs.is_file()
    assert cif_dir.is_dir()

    with pdbs.open("rt") as f:
        pdb_list = [line.strip() for line in f.readlines() if not "#" in line]

    total_proteins = len(pdb_list)

    failed_proteins = pymp.shared.list()

    with pymp.Parallel() as p:
        for idx in p.xrange(total_proteins):
            pdb = pdb_list[idx]

            try:
                init_cif = process_pdb(pdb, cif_dir, pdb_redo)

                freesasa_cif = run_freesasa(init_cif)

                if not freesasa_cif.exists() or freesasa_cif.stat().st_size == 0:
                    print(f"Freesasa Failed: {init_cif.stem}")
                    freesasa_cif.unlink(missing_ok=True)
                    raise FreeSASAException(init_cif.name)

            except Exception:
                failed_proteins.append(pdb)
                p.print(
                    f"Added failed protein {pdb}. Total count: {len(failed_proteins)}"
                )

            finally:
                p.print(
                    f"Thread: {p.thread_num} Finished protein {pdb}: {idx+1} out of {total_proteins}"
                )

    log_failed_proteins = cif_dir / "failed_freesasa_pdbs.txt"

    if len(failed_proteins) > 0:
        with log_failed_proteins.open("wt") as f:
            f.write(f"#Dataset: {pdb_list}\n")
            f.write("\n".join(failed_proteins))
