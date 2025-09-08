#!/bin/bash

# load python environment.
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`

# skyllh.
export PYTHONPATH="/data/user/liruohan/software/skyllh":$PYTHONPATH

# point to analysis code repo
code_path=/data/user/$USER/dm_model_stacking

export PYTHONPATH=$PYTHONPATH:$code_path/compute_sens_dp_from_trials

# stacking trials without NGC 1068
## bkg_trials.
bkg_trials_path=/data/user/$USER/dm_model_stacking/trials/1TeV/bkg
## sig_trials.
sig_trials_path=/data/user/$USER/dm_model_stacking/trials/1TeV/sig

################################
# stacking trials with NGC 1068.
## bkg_trials.
#bkg_trials_path=/data/ana/analyses/NuSources/2022_Seyfert_Galaxies_Analysis/trials/stacking/IC2011_IC2021/include_NGC1068/bkg
## sig_trials.
#sig_trials_path=/data/ana/analyses/NuSources/2022_Seyfert_Galaxies_Analysis/trials/stacking/IC2011_IC2021/include_NGC1068/nsbias
################################

# start running.
python $code_path/compute_sens_dp_from_trials/sens_dp_from_trials_stacking.py --bkg_trials_dir $bkg_trials_path --sig_trials_dir $sig_trials_path