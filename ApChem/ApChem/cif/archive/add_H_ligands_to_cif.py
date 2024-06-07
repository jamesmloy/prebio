from pathlib import Path
from collections import namedtuple
from typing import NamedTuple, Iterable, Union
import argparse
import re

import gemmi


from ApChem.cif.archive.download_protein_ligands import (
    fetch_non_polymer_chem_comp, filter_ligands
)

from ApChem.cif import ChemComp, Ligand, LigandData, Atom, AtomSiteAtom



def cli():

    parser = argparse.ArgumentParser()

    parser.add_argument('-c','--cif-file', required=True, type=str)

    parser.add_argument('-d','--pdb-code-dir', required=False, default=Path().cwd(), type=Path)

    parser.add_argument('-m','--mol2-dir', required=False, default=None, type=Path)

    parser.add_argument('-r','--file-read-regex', required=False, default='_*_*_*[0-9]_H', type=str)

    args = parser.parse_args()

    args.cif_file = args.pdb_code_dir / args.cif_file

    if args.mol2_dir is None: 
        args.mol2_dir = args.pdb_code_dir / 'mol2'

    if not args.cif_file.exists():
        raise FileNotFoundError(f'Input cif file does not exist: {args.cif_file}')
    
    if args.cif_file.suffix != '.cif':
        raise RuntimeError(f'Input file must be a cif file (*.cif): {args.cif_file}')

    pdb_code = args.cif_file.stem[:4]
    file_regex = args.file_read_regex

    args.mol2_files = [mol2_file for mol2_file in args.mol2_dir.glob(f'{pdb_code}{file_regex}.mol2')]

    return args



def prepare_mol2_ligands(mol2_files):

    mol2_ligands = []
    for mol2 in mol2_files:
        _, *ligand, _ = mol2.stem.split('_')

        mol2_ligands.append(Ligand(*ligand))

    return mol2_ligands


def find_mmcif_tag_idx(mmcif_category, tag):

    tags = [tag.split('.')[1] for tag in mmcif_category.tags]

    try: 
        return tags.index(tag)

    except ValueError:
        return None


def find_H_insertion_idx(block, mol2_ligands):

    atom_site = block.find_mmcif_category("_atom_site.")

    comp_id_idx = find_mmcif_tag_idx(atom_site, 'auth_comp_id')
    asym_id_idx = find_mmcif_tag_idx(atom_site, 'label_asym_id')
    seq_id_idx  = find_mmcif_tag_idx(atom_site, 'auth_seq_id')

    insert_idx = []
    for lig in mol2_ligands:
        atom_site_iter = iter(atom_site)

        for atom in atom_site_iter:
            if (atom[comp_id_idx] == lig.comp_id and 
                atom[asym_id_idx] == lig.asym_id and 
                atom[seq_id_idx]  == lig.seq_id):

                for atom in atom_site_iter:
                    if any([atom[comp_id_idx] != lig.comp_id,  
                            atom[asym_id_idx] != lig.asym_id,
                            atom[seq_id_idx]  != lig.seq_id]):

                        insert_idx.append((atom.row_index))

                        break

                break
                
    return insert_idx



def fetch_H_atoms_from_mol2(mol2_file):

    mol2_file = Path(mol2_file)

    if not mol2_file.exists():
        raise FileNotFoundError(f'mol2 file does not exist: {mol2_file}')

    if mol2_file.suffix != '.mol2':
        raise RuntimeError(f'mol2 file does not have a .mol2 suffix: {mol2_file}')

    H_atoms = []
    with mol2_file.open('rt') as f:
        atom_rti = re.compile('@<TRIPOS>\s*ATOM')
        filter_line = False

        for line in f.readlines():
            if line.startswith('@<TRIPOS>') and atom_rti.search(line):
                filter_line = True
                continue

            elif filter_line and line.startswith('@<TRIPOS>'):
                break 

            # line idx : data
            # 1 : atom id
            # 2 : x
            # 3 : y
            # 4 : z
            # 5 : type symbol
            # 7 : chem comp id

            if filter_line and line.split()[1].startswith('H'):
                line = line.split()

                if len(line[5]) != 1 or len(line[7]) > 3 or line[7].isdecimal():
                    raise RuntimeError(f"mol2 file({mol2_file}) is improperly formatted.\nCan't create H atom line: {line}" )

                H_atoms.append(Atom(line[1], line[2], line[3], line[4], line[5], line[7]))

    return H_atoms
        


def same_ligand(lig1, lig2):

    if lig1.comp_id == lig2.comp_id:
        if lig1.asym_id == lig2.asym_id:
            if lig1.seq_id == lig2.seq_id:
                return True

    return False



def H_atom_site_template(block, template_atom_idx):

    atom_site = block.find_mmcif_category("_atom_site.")

    template_atom = atom_site[template_atom_idx]

    group_PDB = 'HETATM'
    # id
    type_symbol = 'H'
    # label_atom_id
    label_alt_id = '.'    
    label_comp_id = template_atom[find_mmcif_tag_idx(atom_site, 'label_comp_id')]
    label_asym_id = template_atom[find_mmcif_tag_idx(atom_site, 'label_asym_id')]
    label_entity_id = template_atom[find_mmcif_tag_idx(atom_site, 'label_entity_id')]
    label_seq_id = template_atom[find_mmcif_tag_idx(atom_site, 'label_seq_id')] # Can be an empty string
    pdbx_PDB_ins_code = '?'
    # Cartn_x
    # Cartn_y
    # Cartn_z
    occupancy = '1.00'
    B_is_or_equiv = '10.00'
    pdbx_formal_charge = '?'
    auth_seq_id = template_atom[find_mmcif_tag_idx(atom_site, 'auth_seq_id')]
    auth_comp_id = template_atom[find_mmcif_tag_idx(atom_site, 'auth_comp_id')]
    auth_asym_id = template_atom[find_mmcif_tag_idx(atom_site, 'auth_asym_id')]
    # auth_atom_id
    pdbx_PDB_model_num = template_atom[find_mmcif_tag_idx(atom_site, 'pdbx_PDB_model_num')]

    H_atom_template = AtomSiteAtom(
            group_PDB, 
            None, # replace
            type_symbol,
            None, # replace
            label_alt_id,
            label_comp_id,
            label_asym_id,
            label_entity_id,
            label_seq_id,
            pdbx_PDB_ins_code,
            None, # replace
            None, # replace
            None, # replace
            occupancy,
            B_is_or_equiv,
            pdbx_formal_charge,
            auth_seq_id,
            auth_comp_id,
            auth_asym_id,
            None, # replace
            pdbx_PDB_model_num
        )

    return H_atom_template



def create_H_atom_site_rows(block, mol2_file, mol2_ligand, ins_H_idx):

    mol2_file = Path(mol2_file)

    _, *ligand, _ = mol2_file.stem.split('_')

    file_ligand = Ligand(*ligand)

    if not same_ligand(file_ligand, mol2_ligand):
        raise RuntimeError(
            f'Mol2 file {mol2_file.stem} is the incorrect ligand for the insertion of {mol2_ligand}'
        )

    atom_site = block.find_mmcif_category("_atom_site.")

    last_lig_atom = atom_site[ins_H_idx -1]

    serial = int(last_lig_atom[find_mmcif_tag_idx(atom_site, 'id')])

    if last_lig_atom[find_mmcif_tag_idx(atom_site, 'label_asym_id')] != mol2_ligand.asym_id:
        print(f"Lig atom:  {last_lig_atom[find_mmcif_tag_idx(atom_site, 'label_asym_id')]}")
        print(f'mol2_atom: {mol2_ligand.asym_id}')
        raise RuntimeError('Label asym id is not equal between cif and mol2')

    if last_lig_atom[find_mmcif_tag_idx(atom_site, 'auth_seq_id')] != mol2_ligand.seq_id:
        print(f"Lig atom:  {last_lig_atom[find_mmcif_tag_idx(atom_site, 'auth_seq_id')]}")
        print(f'mol2_atom: {mol2_ligand.seq_id}')
        raise RuntimeError('Auth seq id is not equal between cif and mol2')

    if last_lig_atom[find_mmcif_tag_idx(atom_site, 'auth_comp_id')] != mol2_ligand.comp_id:
        print(f"Lig atom:  {last_lig_atom[find_mmcif_tag_idx(atom_site, 'auth_comp_id')]}")
        print(f'mol2_atom: {mol2_ligand.comp_id}')
        raise RuntimeError('Auth comp id is not equal between cif and mol2')


    H_atom_template = H_atom_site_template(block, ins_H_idx -1)

    mol2_H_atoms = fetch_H_atoms_from_mol2(mol2_file)

    H_atom_site_rows = [
        list(H_atom_template._replace(
            id=str(serial + count + 1),
            label_atom_id=H_atom.atom_id,
            Cartn_x=H_atom.x,
            Cartn_y=H_atom.y,
            Cartn_z=H_atom.z,
            auth_atom_id=H_atom.atom_id
        ))
        for count, H_atom in enumerate(mol2_H_atoms)
    ]

    return H_atom_site_rows



def add_ligand_H_to_atom_site(block, mol2_file, mol2_ligand, ins_H_idx):

    H_atom_site_rows = create_H_atom_site_rows(block, mol2_file, mol2_ligand, ins_H_idx)

    atom_site = block.find_mmcif_category("_atom_site.")

    atom_site_before = [list(row) for idx, row in enumerate(atom_site) if idx < ins_H_idx]

    atom_site_after = [list(row) for idx, row in enumerate(atom_site) if idx > ins_H_idx - 1]

    updated_atom_site_rows = atom_site_before + H_atom_site_rows + atom_site_after

    tags = [tag.split('.')[1] for tag in atom_site.tags]

    new_atom_site = block.init_mmcif_loop(atom_site.get_prefix(), tags)

    for row in updated_atom_site_rows:

        new_atom_site.add_row(row)



def reset_atom_site_id_column(block):

    atom_site = block.find_mmcif_category("_atom_site.")

    for idx, row in enumerate(atom_site):
        row[find_mmcif_tag_idx(atom_site, 'id')] = str(idx + 1)



def ligands_found_in_cif(block, mol2_ligands: Iterable[NamedTuple], if_all_found=False) \
    -> Union[Iterable, bool]:

    cif_ligand_ids = fetch_non_polymer_chem_comp(block)

    cif_ligands = sorted(filter_ligands(block, cif_ligand_ids), key= lambda x: (x.asym_id, int(x.seq_id)))

    ligand_files_found = []
    for cif in cif_ligands:
        found = any(
            mol2.asym_id == cif.asym_id and mol2.seq_id==cif.seq_id 
            for mol2 in mol2_ligands
        )

        if found:
            ligand_files_found.append(cif)

    if if_all_found:
        all_found = len(cif_ligands) == len(ligand_files_found)
        return all_found

    return ligand_files_found



def mol2_ligands_not_in_cif(block, mol2_ligands: Iterable[NamedTuple]) -> Iterable:

    cif_ligand_ids = fetch_non_polymer_chem_comp(block)

    cif_ligands = filter_ligands(block, cif_ligand_ids)

    random_mol2_files = []
    for mol2 in mol2_ligands:
        absent = all(
            mol2.asym_id != cif.asym_id and mol2.seq_id != cif.seq_id 
            for cif in cif_ligands
        )

        if absent: 
            random_mol2_files.append(mol2)

    return random_mol2_files



def add_ligand_hydrogen_rows_to_cif(cif_file, mol2_files):

    mol2_ligands = prepare_mol2_ligands(mol2_files)

    doc = gemmi.cif.read_file(str(cif_file))
    block = doc.sole_block()

    if mol2_ligands_not_in_cif(block, mol2_ligands):
        raise RuntimeError(
            f"An input mol2 file is not found within cif file's \
                chemp_comp category: {mol2_ligands_not_in_cif()}") 

    mol2_H_ins_idx = find_H_insertion_idx(block, mol2_ligands)

    all_ligand_data = [LigandData(*data) for data in zip(mol2_files, mol2_ligands, mol2_H_ins_idx)]

    # You must insert starting at the bottom of the cif file to not corrupt the H_insert_idx values
    all_ligand_data.sort(key=lambda x: x.H_insert_idx, reverse=True) 

    for ligand in all_ligand_data:
        add_ligand_H_to_atom_site(block, *ligand)

    reset_atom_site_id_column(block)

    doc.write_file(str(cif_file.with_suffix('.H.cif')), gemmi.cif.Style.Pdbx)

    if not ligands_found_in_cif(block, mol2_ligands, if_all_found=True):
        print('Not all ligands present in the cif file were passed as mol2 files.')
        print('Some ligands will be lacking hydrogens.')



if __name__=="__main__":

    args = cli()

    add_ligand_hydrogen_rows_to_cif(args.cif_file, args.mol2_files)