#!/bin/bash

# 06-21-24
# test script to generate the cif file and the boxes and graphs from a directory of cif files

cif_dir=../data/cif/
box_dir=../data/box/
graph_dir=../data/graph/

python preprocess.py \
    --image apchem:0424 \
    --in-folder $cif_dir \
    --num-proc 4

python generate_microenvs.py \
    --image prebio:final_build \
    --in-folder $cif_dir \
    --out-folder $box_dir \

python generate_microenvs.py \
    --image prebio:final_build \
    --in-folder $cif_dir \
    --out-folder $graph_dir \
    --graph 

python generate_microenvs.py \
    --image prebio:final_build \
    --in-folder $cif_dir \
    --out-folder $box_dir \
    --channel-first