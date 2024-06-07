from pathlib import Path
from typing import Callable

from ApChem.exceptions import DirectoryNotFoundError, DuplicatePDBError


def find_missing_pdbs(output_dir: Path, pdb_list: Path, filter: str) -> None:

    if not output_dir.exists() or not output_dir.is_dir():
        raise DirectoryNotFoundError(
            f"{output_dir} was not found or is not a directory."
        )

    if not pdb_list.exists() or not pdb_list.is_file():
        raise FileNotFoundError(f"{pdb_list} was not found or is not a file.")

    dataset = set(pdb_list.read_text().split("\n"))

    for cif in output_dir.glob(filter):

        pdb_id = cif.stem[:4]

        match = [pdb for pdb in dataset if pdb_id in pdb]

        if len(match) == 0:
            print(f"No match found for cif: {cif} using query: {pdb_id}")

        elif len(match) > 1:
            raise DuplicatePDBError(
                f"Multiple matches in dataset for cif: {cif}\nMatches: {match}\n name: {pdb_id}"
            )

        else:
            dataset.remove(match[0])

    return dataset


def mv_files(from_dir: Path, to_dir: Path, glob: str = "*.fw2.sasa.cif") -> None:

    assert isinstance(from_dir, Path) and from_dir.exists()

    assert isinstance(to_dir, Path)

    if not to_dir.exists():
        to_dir.mkdir(0o774, parents=True, exist_ok=True)

    for idx, f in enumerate(from_dir.glob(glob)):

        out_f = to_dir / f.name

        f.replace(out_f)

        print(f"{idx:6d} Moved: {f} to {out_f}")


def rename_cif_by_pdb_id(
    cif: Path, extension=".fw2.sasa.cif", suffix: str = "", **kwargs
) -> Path:

    assert isinstance(cif, Path)
    assert cif.exists() and cif.is_file()

    new_name = cif.parent / (cif.stem[:4].lower() + suffix + extension)

    cif.rename(new_name)

    return new_name


def rename_cif_by_split(
    cif: Path, delimiter="_", extension=".fw2.sasa.cif", suffix: str = "", **kwargs
) -> Path:

    assert isinstance(cif, Path)
    assert cif.exists() and cif.is_file()

    new_name = cif.parent / (
        cif.stem.split(delimiter, maxsplit=1)[0].lower() + suffix + extension
    )

    cif.rename(new_name)

    return new_name


def rename_cif_dir(
    cif_dir: Path,
    filter: str = "*.fw2.sasa.cif",
    rename_func: Callable = rename_cif_by_pdb_id,
    suffix: str = "",
    extension: str = ".fw2.sasa.cif",
    **kwargs,
) -> None:

    assert isinstance(cif_dir, Path)
    assert cif_dir.exists() and cif_dir.is_dir()

    for cif in cif_dir.glob(filter):
        rename_func(cif, suffix=suffix, extension=extension, **kwargs)

