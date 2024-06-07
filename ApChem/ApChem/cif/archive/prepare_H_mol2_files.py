import argparse

from pathlib import Path

from ApChem.cif.archive.correct_mol2 import add_H_to_mol2



def cli():

    parser = argparse.ArgumentParser()

    parser.add_argument('-d','--dataset-dir', required=True, type=Path)

    args = parser.parse_args()

    return args


def add_H_to_mol2_dataset(dataset_dir, glob='*/mol2/*[0-9]'):

    dataset_dir = Path(dataset_dir)

    for mol2 in dataset_dir.rglob(f'{glob}.mol2'):
        add_H_to_mol2(mol2)



if __name__=="__main__":

    args = cli()

    add_H_to_mol2_dataset(args.dataset_dir)