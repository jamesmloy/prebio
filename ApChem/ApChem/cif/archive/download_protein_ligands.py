from pathlib import Path

from wget import download, bar_thermometer
import gemmi

from ApChem.cif import (
    Ligand, fetch_non_polymer_chem_comp
)


def filter_ligands(block, ligand_ids):

    atom_site = block.find("_atom_site.", ["auth_comp_id", "label_asym_id", "auth_seq_id"])

    return {Ligand(*row) for row in atom_site if row[0] in ligand_ids}



def download_ligands(cif_file, data_dir, file_type='mol2'):

    block = gemmi.cif.read_file(str(cif_file)).sole_block()

    pdb_id = block.name.lower()

    ligand_ids = fetch_non_polymer_chem_comp(block)

    all_ligands = filter_ligands(block, ligand_ids)

    url = 'https://models.rcsb.org/v1/{}/ligand?auth_seq_id={}&label_asym_id={}&encoding={}&filename={}' 

    for ligand in all_ligands:

        name = f'{pdb_id.lower()}_{ligand.comp_id}_{ligand.asym_id}_{ligand.seq_id}.{file_type}' 

        print(f'\nDownloading ligand: {ligand}\nWriting to: {data_dir/ name}')

        # TODO refactor with requests?
        download(
            url.format(pdb_id.lower(), ligand.seq_id, ligand.asym_id, file_type, name), 
            out=str(data_dir / name),
            bar=bar_thermometer    
        )  


def download_proteins(pdb_list_file: Path, data_dir: Path, pdb_type='both', api='rcsb'):

    import pymp

    url = '' 
    if api == 'rcsb':
        url = 'http://files.rcsb.org/download/{}'

    elif api == 'pdb-redo':
        url = 'https://pdb-redo.eu/db/{}/{}'

    if pdb_type is 'both':
        pdb_type = ['pdb', 'cif']

    with pdb_list_file.open('rt') as f:
        pdb_list = [line.strip() for line in f.readlines() if not '#' in line]

    total_proteins = len(pdb_list)

    failed_proteins = pymp.shared.list()
    with pymp.Parallel() as p:
        for idx in p.xrange(total_proteins):
            pdb_id = pdb_list[idx]
            try: 
                # TODO refactor with requests?
                if api == 'rcsb':
                    if 'pdb' in pdb_type:
                        download(
                            url.format(pdb_id.upper() + '.pdb'), 
                            out=str(data_dir/ (pdb_id.lower() + '.pdb')),
                            bar=bar_thermometer
                        )

                    if 'cif' in pdb_type:
                        download(
                            url.format(pdb_id.upper() + '.cif'), 
                            out=str(data_dir/ (pdb_id.lower() + '.cif')),
                            bar=bar_thermometer    
                        )

                elif api == 'pdb-redo':
                    if 'pdb' in pdb_type:
                        filename = pdb_id + '_final_tot.pdb'
                        download(
                            url.format(pdb_id, filename), 
                            out=str(data_dir / filename),
                            bar=bar_thermometer    
                        )

                    if 'cif' in pdb_type:
                        filename = pdb_id + '_final.cif'
                        download(
                            url.format(pdb_id, filename), 
                            out=str(data_dir / filename),
                            bar=bar_thermometer    
                        )

            except Exception as e:
                failed_proteins.append(pdb_id)
                p.print(f'Error with protein ({pdb_id}):\n{e}\n')

            finally:
                p.print(f'Thread: {p.thread_num} Finished protein ({pdb_id}): {idx+1} out of {total_proteins}')

    log_failed_proteins = data_dir / 'failed_pdbs.txt'

    if len(failed_proteins) > 0:
        with log_failed_proteins.open('wt') as f:
            f.write(f'#Dataset: {pdb_list_file}\n')
            f.write('\n'.join(failed_proteins))



def setup_pdb_code_dir(pdb_code_dir, overwrite=False):

    pdb_code_dir = Path(pdb_code_dir)
    mol2_dir = pdb_code_dir / "mol2"

    if not overwrite:
        try: 
            mol2_dir.mkdir(0o744, parents=True, exist_ok=False)
        
        except FileExistsError:
            print(f'Directory already exist: {mol2_dir}   Skipping mkdir...')

    else:
        mol2_dir.mkdir(0o744, parents=True, exist_ok=True)



def download_protein_ligand_files(pdb_code, dataset_dir, api, overwrite=False):

    data_dir = dataset_dir / pdb_code
    mol2_dir = data_dir / 'mol2'

    setup_pdb_code_dir(data_dir, overwrite)

    assert(data_dir.exists() and mol2_dir.exists())

    try: 
        download_proteins([pdb_code], data_dir, 'cif', api)

    except Exception as e: 
        print(f'{pdb_code} errored out: {e.msg}')

    try:
        cif_file = [cif for cif in data_dir.glob(f'{pdb_code}*.cif')][0]

        download_ligands(cif_file, mol2_dir)

    except Exception as e:
        print(f'Error in downloading mol2 files for protein: {cif_file}\n{e.msg}')
    
