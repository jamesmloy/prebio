from pathlib import Path
import argparse

from ApChem.cif.archive.add_H_ligands_to_cif import add_ligand_hydrogen_rows_to_cif

def cli():

    parser = argparse.ArgumentParser()

    parser.add_argument('-d','--dataset-dir', required=True, type=Path)

    parser.add_argument('--cif-glob', required=False, default='????', type=str)

    parser.add_argument('--mol2-regex', required=False, default='_*_*_*[0-9]_H', type=str)

    args = parser.parse_args()

    return args



def add_ligand_H_atoms_to_cif_dataset(dataset_dir, cif_glob='????', mol2_regex='_*_*_*[0-9]_H'):

    dataset_dir = Path(dataset_dir)

    for cif_file in dataset_dir.rglob(f'{cif_glob}.cif'):
        mol2_dir = cif_file.parent / 'mol2'
        pdb_code = cif_file.stem[:4]

        mol2_files = [mol2_file for mol2_file in mol2_dir.glob(f'{pdb_code}{mol2_regex}.mol2')]

        add_ligand_hydrogen_rows_to_cif(cif_file, mol2_files)



if __name__=="__main__":

    args = cli()

    add_ligand_H_atoms_to_cif_dataset(**vars(args)) 
