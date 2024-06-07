from pathlib import Path
from sys import stdout
import subprocess


def generate_cif_docker(
    pdb_code: str, cif_dir: Path, num_threads: int = 2, docker_img: str = "apchem", pdb_redo=True
) -> None:

    assert len(pdb_code) == 4
    assert isinstance(cif_dir, Path)

    cif_dir.mkdir(0o777, parents=True, exist_ok=True)

    cmd = f"""
            docker run 
                --rm 
                --env OMP_NUM_THREADS={num_threads} 
                -v {cif_dir.resolve()}:/data
                {docker_img} 
                generate_cif.py {pdb_code.lower()} -d /data 
    """.strip()

    if pdb_redo: 
        cmd += " --pdb-redo"

    print(f"Running: {cmd.strip()}")

    p = subprocess.Popen(cmd.strip().split(), stdout=stdout, stderr=subprocess.STDOUT,)

    statuscode = p.wait()

    if statuscode != 0:
        raise subprocess.SubprocessError("generate_cif.py failed")


if __name__ == "__main__":

    pdb_code = "6ij6"

    cif_dir = Path("../../data/proteins/cif/docker")

    generate_cif_docker(pdb_code, cif_dir, num_threads=4, pdb_redo=True)

