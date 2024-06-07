from pathlib import Path
from datetime import datetime
from shutil import copy

import pandas as pd

from ApChem.cif import (
    reorder_pdb_redo_columns,
    fetch_cif_from_api,
    run_chimX_fw2_freesasa_pipeline,
)


class PiscesCIFDatasetProducer:
    def __init__(self, dataset: Path, redo_dir: Path, out_dir: Path):
        """dataset must be a pisces pdbcull file that is not the fasta one"""

        assert isinstance(dataset, Path)
        assert dataset.exists() and dataset.is_file()

        assert isinstance(redo_dir, Path)
        assert redo_dir.exists() and redo_dir.is_dir()

        self.dataset = dataset
        self.redo_dir = redo_dir

        self.out_dir = out_dir
        out_dir.mkdir(0o777, parents=True, exist_ok=True)

        df = pd.read_csv(self.dataset, delimiter=r"\s+")

        df["PDB_ID"], df["PDB_chain"] = df["PDBchain"].str[:4], df["PDBchain"].str[4:]

        self.inp_df = df[
            ["PDB_ID", "PDB_chain", "len", "method", "resol", "rfac", "freerfac"]
        ]

        self.pdb_ids = set(df.PDB_ID.values)

    def create_dataset(self):

        import pymp

        pdb_ids = list(self.pdb_ids)
        total_proteins = len(pdb_ids)

        assert total_proteins > 0

        failed_proteins = pymp.shared.list()

        with pymp.Parallel() as p:
            for idx in p.xrange(total_proteins):

                pdb = pdb_ids[idx].lower()

                cif = self.redo_dir / f"{pdb}_final.cif"

                new_cif = self.out_dir / cif.name

                rcsb_cif = self.out_dir / f"{pdb}.cif"

                if new_cif.exists() or rcsb_cif.exists():
                    p.print(f"Already copied: {pdb}")

                elif cif.exists():
                    p.print(f"Found REDO: {pdb}")
                    copy(cif, new_cif)
                    reorder_pdb_redo_columns(new_cif)

                else:
                    try:
                        new_cif = fetch_cif_from_api(pdb.lower(), self.out_dir)

                    except FileNotFoundError as e:
                        failed_proteins.append(f"{pdb.lower()},")
                        p.print(
                            f"Failed to download protein: {pdb} Total fail count: {len(failed_proteins)} Error:\n{e}"
                        )

                try:
                    # check if we have already run the pipeline for this
                    processed_file = list(self.out_dir.glob(f"{pdb}*.fw2.sasa.cif"))

                    # if the list is empty, then we have not processed it
                    if not processed_file:
                        run_chimX_fw2_freesasa_pipeline(new_cif)
                    else:
                        p.print(
                            f"Already processed pdb {pdb} -- found {processed_file[0]}"
                        )
                        if len(processed_file) > 1:
                            p.print(
                                f"WARNING: Found multiple final files for pdb {pdb}:\n{processed_file}"
                            )
                except Exception as e:
                    failed_proteins.append(f"{pdb.lower()},")
                    p.print(
                        f"Failed to process cif file: {pdb} Total fail count: {len(failed_proteins)} Error:\n{e}\n===================="
                    )

                cnt = int(idx + 1)
                p.print(
                    f"Thread {p.thread_num:3d} Finished: {cnt:6d} out of {total_proteins}"
                )

        now = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")

        if len(failed_proteins) > 0:
            log_failed_proteins = self.out_dir / "failed_proteins.txt"

            with log_failed_proteins.open("wt") as f:
                f.write(f"# PISCES CIF Producer - proteins that failed downloading\n")
                f.write(f"# Generated: {now}\n")
                f.write(f"# CIF Directory: {self.redo_dir}\n")
                f.write(f"# Number of CIF: {total_proteins}\n")
                f.write("\n".join(failed_proteins))


if __name__ == "__main__":

    pdb_ids = Path(
        "../../data/pdb_datasets/pisces_10_29_21/cullpdb_pc50.0_res0.0-3.0_noBrks_len40-10000_R0.3_Xray_d2021_10_29_chains22759"
    )
    redo_dir = Path("../../data/proteins/cif/ab_5A_renumbered")
    out_dir = Path("./cif")

    producer = PiscesCIFDatasetProducer(pdb_ids, redo_dir, out_dir)

    producer.create_dataset()
