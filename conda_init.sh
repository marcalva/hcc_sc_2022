#!/bin/bash

# export HOME=/u/project/pajukant/malvarez/

export conda_start="/u/project/pajukant/malvarez/anaconda3/etc/profile.d/conda.sh"

export init_shell="bash"

if [[ -z $conda_start ]]
then
    conda init $init_shell
else
    . $conda_start
fi

# to run a conda environemnt 'env' for example, type
# . conda_init.sh
# conda activate env
# ...
# conda deactivate

