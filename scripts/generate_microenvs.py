from pathlib import Path
from subprocess import run
from argparse import ArgumentParser


def call_generate_boxes(
    image: str, in_folder: Path, out_folder: Path, channel_first: bool, graph: bool
):
    cmd = (
        "docker run "
        f"-v {in_folder.resolve()}:/input "
        f"-v {out_folder.resolve()}:/output "
        f" -t {image} "
        "python /generate_boxes_docker.py "
        "--in-folder /input "
        "--out-folder /output"
    )

    if channel_first:
        cmd += " --channel-first"
    if graph:
        cmd += " --graph"

    run(cmd, shell=True)


def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--in-folder",
        type=Path,
        help="Folder containing the cif files to run predictions on.",
        required=True,
    )

    parser.add_argument(
        "--out-folder",
        type=Path,
        help="The folder where the predictions will be saved",
        required=True,
    )
    parser.add_argument(
        "--image",
        help="Docker image to use",
        default="prebio:latest",
    )

    parser.add_argument(
        "--channel-first",
        action="store_true",
        help="Channel first data ordering (default channel last, ignored for graphs)",
    )

    parser.add_argument(
        "--graph",
        action="store_true",
        help="Generate graph based microenvironment",
    )

    return parser.parse_args()


def main(
    in_folder: Path, out_folder: Path, image: str, channel_first: bool, graph: bool
):
    call_generate_boxes(
        image=image,
        in_folder=in_folder,
        out_folder=out_folder,
        channel_first=channel_first,
        graph=graph,
    )


if __name__ == "__main__":
    args = vars(get_args())
    main(**args)
