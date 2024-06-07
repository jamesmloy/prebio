from pathlib import Path
import sys
import time

import chargefw2_python



def time_it(func):

    def inner(*args, **kwargs):
        print(f'Running {func.__name__} function.')
        start = time.time()
        out = func(*args, **kwargs)
        print(f'Done running {func.__name__} {round(time.time() - start, 3)} seconds')
        return out

    return inner


out_dir = Path('../../../data/proteins/fw2/rcsb/')
pdb_id = Path(sys.argv[1]).stem[:4]
out_dir = out_dir / pdb_id

out_dir.mkdir(0o744,parents=True, exist_ok=True)

protein = time_it(chargefw2_python.Molecules)(sys.argv[1])
methods = time_it(chargefw2_python.get_suitable_methods)(protein)
charges = time_it(chargefw2_python.calculate_charges)(protein, "eem", 'EEM_10_Cheminf_b3lyp_aim')

chargefw2_python.write_cif(protein, charges, sys.argv[1], str(out_dir))