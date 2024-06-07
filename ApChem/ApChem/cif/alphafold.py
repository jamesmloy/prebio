from pathlib import Path

import numpy as np
import gemmi


def select_most_populated_block(cif_inp: Path) -> Path:

    doc = gemmi.cif.read(str(cif))

    new_doc = gemmi.cif.Document()

    print(f"Number of blocks: {len(doc)}")

    block_sizes = [len(block.find_mmcif_category("_atom_site.")) for block in doc]

    max_idx = np.argmax(block_sizes)

    new_doc.add_copied_block(doc[max_idx])

    fixed_cif = cif.parent / f"{cif.stem}_fixed.cif"

    new_doc.write_file(str(fixed_cif))

    return fixed_cif



if __name__=="__main__":

    cif = Path("/alphafold/jimmy/A7V3Isolde_real_space_refined_020-coot-1_real_space_refined_026.cif")

    select_most_populated_block(cif)

