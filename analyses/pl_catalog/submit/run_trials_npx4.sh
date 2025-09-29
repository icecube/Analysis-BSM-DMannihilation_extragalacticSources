#!/bin/bash -x

echo $HOSTNAME

#photosplines etcs
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`

# add SkyLLH and i3skyllh to env
export PYTHONPATH="/data/user/liruohan/software/skyllh":$PYTHONPATH
export PYTHONPATH="/data/user/liruohan/software/i3skyllh":$PYTHONPATH
export PYTHONPATH=/cvmfs/icecube.opensciencegrid.org/users/tkontrimas/software/photospline_v2.2.0/:$PYTHONPATH


# add seyfert model code
#export PYTHONPATH="/data/user/$USER/seyfert_galaxies/analysis_code/2022_Seyfert_Galaxies_Analysis/model_code/":$PYTHONPATH

#export HDF5_USE_FILE_LOCKING=FALSE

# code repo.
code_path=/data/user/liruohan/Analysis-BSM-DMannihilation_extragalacticSources/analyses/pl_catalog

python $code_path/generate_trials.py $@