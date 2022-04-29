#!/bin/bash

. ../conda_init.sh

conda activate liver_snrna_hcc

Rscript add_conda_pkg.R

# add python packages
pip install scanpy

conda deactivate

