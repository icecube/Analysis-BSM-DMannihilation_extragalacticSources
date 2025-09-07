#!/bin/bash

set -x

# load python environment
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/setup.sh`

export HDF5_USE_FILE_LOCKING='FALSE'

# point to analysis code repo
#code_path=/data/user/liruohan/dataset/submit_data

# run script
source '/data/user/hmniederhausen/combo_V01-00-00/build/env-shell.sh' python /data/user/liruohan/prepare_Analysis-BSM-DMannihilation_extragalacticSources/data_generation/submit_data/generate_data_files.py $@