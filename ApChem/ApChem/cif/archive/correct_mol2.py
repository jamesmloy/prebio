from pathlib import Path
from openbabel import pybel

import argparse
import re




def cli():

    parser = argparse.ArgumentParser()

    parser.add_argument('-f','--mol2-filename', required=True, type=str)

    parser.add_argument('-d','--data-dir', required=False, default=Path().cwd(), type=Path)

    args = parser.parse_args()

    args.mol2_file = args.data_dir / args.mol2_filename

    if not args.mol2_file.exists():
        raise FileNotFoundError(f'Input mol2 file does not exist: {args.mol2_file}')
    
    if args.mol2_file.suffix != '.mol2':
        raise RuntimeError(f'Input file must be a mol2 file (*.mol2): {args.mol2_file}')

    return args



def convert_mol2_to_cif(mol2_file):
    mol2_f = Path(mol2_file)

    cif_f = mol2_f.with_suffix('.cif')

    mol = next(pybel.readfile("mol2", str(mol2_f)))

    mol.write('cif', str(cif_f), overwrite=True)



def correct_mol2_H_atom_name(line: list, H_idx: int):
    line[1] = f'H{H_idx}'


def correct_mol2_chem_comp_name(line):
    line[7] = line[7][:3]



def add_H_to_mol2(mol2_file):
    mol_f = Path(mol2_file)

    if not mol_f.exists():
        raise FileNotFoundError(f'file: {mol2_file} does not exist.')
    
    if mol_f.suffix != '.mol2':
        raise RuntimeError(f'File: {mol2_file} is not a .mol2 file but a {mol2_file.suffix}')

    mol_h_f = mol_f.parent / (mol_f.stem + "_H.tmp" + mol_f.suffix)

    corr_mol_h_f = mol_f.parent / (mol_f.stem + "_H" + mol_f.suffix)

    mol = next(pybel.readfile("mol2", str(mol_f)))
    mol.addh()
    mol.write('mol2', str(mol_h_f), overwrite=True)


    with mol_h_f.open('rt') as rf, corr_mol_h_f.open('wt') as of:
        atom_rti = re.compile('@<TRIPOS>\s*ATOM')
        edit = False
        H_idx = 1

        for line in rf.readlines():
            if line.startswith('@<TRIPOS>') and atom_rti.search(line):
                edit = True
                of.write(line)
                continue

            elif edit and line.startswith('@<TRIPOS>'):
                edit = False

            line = line.split()


            if edit:
                correct_mol2_chem_comp_name(line)
                if line[5] == 'H':
                    correct_mol2_H_atom_name(line, H_idx)
                    H_idx += 1

            of.write(' '.join(line)+'\n')

    mol_h_f.unlink()

    print(f'Hydrogens added to mol2: {corr_mol_h_f}')




if __name__=="__main__":

    args = cli()

    add_H_to_mol2(args.mol2_file)
