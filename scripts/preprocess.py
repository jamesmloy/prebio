from pathlib import Path
from subprocess import run
from argparse import ArgumentParser


# TODO I dont have the image with process_cif_files.py
def call_process_cif(image: str, in_folder: Path, out_folder: Path, num_proc: int):
    cmd = (
        "docker run "
        f"-v {in_folder.resolve()}:/input "
        f"-v {out_folder.resolve()}:/output "
        f" -t {image} "
        "python process_cif_files.py "
        f"--in-folder /input "
        f"--out-folder /output "
        f"--num-proc {num_proc}"
    )

    run(cmd, shell=True)


def process_cif(image: str, in_folder: Path, num_proc: int):

    cmd = (
            "docker run "
            f"-v {in_folder.resolve()}:/data "
            f"--env PYMP_NUM_THREADS={num_proc} "
            f"-t {image} "
            "python generate_cif_dataset.py "
            f"/data/cifs.txt "
            f"-d /data "
    )

    cifs = [cif.name for cif in in_folder.glob("*.cif")]

    dataset = in_folder / "cifs.txt"

    with dataset.open("w") as f:
        f.write("\n".join(cifs))    
    
    run(cmd, shell=True)


def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--in-folder",
        type=Path,
        help="Folder containing the pdb/cif files to preprocess.",
        required=True,
    )

    parser.add_argument(
        "--out-folder",
        type=Path,
        help="The folder where the cif files will be saved",
        required=False,
    )
    parser.add_argument("-n", "--num-proc", help="The number of processors to use", type=int, default=4)
    parser.add_argument(
        "-i",
        "--image",
        help="Docker image to use",
        default="apchem:latest",
    )

    args = parser.parse_args()

    if args.out_folder is None:
        args.out_folder = args.in_folder

    return parser.parse_args()


def main(
    in_folder: Path,
    out_folder: Path,
    image: str,
    num_proc: int
):
    # call_process_cif(
    #     image=image, in_folder=in_folder, out_folder=out_folder, num_proc=num_proc
    # )

    process_cif(
        image=image, in_folder=in_folder, num_proc=num_proc
    )


if __name__ == "__main__":
    args = vars(get_args())
    main(**args)
