#!/bin/bash

# load python environment.
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`

# skyllh.
export PYTHONPATH="/data/user/tkontrimas/software/li/skyllh":$PYTHONPATH

# point to analysis code repo
code_path=/data/user/liruohan/Analysis-BSM-DMannihilation_extragalacticSources

export PYTHONPATH=$PYTHONPATH:$code_path/compute_sens_dp_from_trials

# model flux ps  
## bkg_trials.
bkg_trials_path=/data/user/$USER/powerlaw/trials/bkg
#bkg_trials_path=/data/ana/analyses/NuSources/2022_Seyfert_Galaxies_Analysis/trials/catalog/model/bkg/ 
## sig_trials.
sig_trials_path=/data/user/$USER/powerlaw/trials/sig

# python $code_path/compute_sens_dp_from_trials/sens_dp_from_trials_ps.py --analysis model --trials_location_bkg $bkg_trials_path --trials_location_signal $sig_trials_path

################################
# powerlaw flux ps 
#bkg_trials_path=/data/ana/analyses/NuSources/2022_Seyfert_Galaxies_Analysis/trials/catalog/powerlaw/bkg/ 
## sig_trials.
#sig_trials_path=/data/ana/analyses/NuSources/2022_Seyfert_Galaxies_Analysis/trials/catalog/powerlaw/signal/

python $code_path/compute_sens_dp_from_trials/sens_dp_from_trials_ps.py --analysis powerlaw --trials_location_bkg $bkg_trials_path --trials_location_signal $sig_trials_path
################################