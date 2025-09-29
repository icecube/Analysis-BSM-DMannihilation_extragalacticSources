#!/usr/bin/env python

import os, sys
import logging
import numpy as np
import pandas as pd

import skyllh
from skyllh.core.random import RandomStateService
from skyllh.core.timing import TimeLord
from skyllh.core.source_model import PointLikeSource
from skyllh.core.utils.analysis import create_trial_data_file
from skyllh.core.config import Config
from skyllh.core.minimizers.iminuit import IMinuitMinimizerImpl

from i3skyllh.datasets import data_samples
from i3skyllh.analyses.trad_stacked_ps import analysis as trad_stacking_analysis

from skyllh.core.debugging import setup_logger, setup_console_handler
#setup_logger('skyllh', logging.INFO)
#setup_console_handler('skyllh', logging.INFO)
#setup_logger('i3skyllh', logging.INFO)
#setup_console_handler('i3skyllh')


import argparse
parser=argparse.ArgumentParser()

parser.add_argument("--mean_ns", "-ns", help="how many signal events to inject on average", type=float, default=0.0)
parser.add_argument("--gamma", "-g", help="spectral index for if injecting from powerlaw", type=float, default=2.0)
parser.add_argument("--n_trials", "-nt", help="number of trials per job", type=int, default=10)
parser.add_argument("--outdir", "-od", help="output folder", type=str, default='/data/user/liruohan/powerlaw_stacking/trials/')
parser.add_argument("--outfile", "-of", help="output file", type=str, default='test.npy')
parser.add_argument("--rss_seed", "-rs", help="random_number_seed for mc generator", type=int, default=0)
parser.add_argument("--dset", "-ds", help="name of the dataset to load for powerlaw analysis", type=str, default='NorthernTracks_v005p01')
parser.add_argument("--base_path", "-bp", help="base path to the folder containing datasets", type=str, default='/data/ana/analyses/')
parser.add_argument("--period", "-p", help="name of the data period to load", type=str, default='IC86_2011_2021')
parser.add_argument("--srcs_table_path", "-bstp", help="where to find sources pandas pkl table", type=str, default='/data/user/liruohan/powerlaw_stacking/sources.pkl')
parser.add_argument("--ncpus", "-nc", help="number of CPUs to use", type=int, default=4)
parser.add_argument("--unblind", help="run the unblinded analysis", default=False, action="store_true")
args=parser.parse_args()
print(args)
cfg=Config()
cfg.set_ncpu(4)

dsc = data_samples[args.dset].create_dataset_collection(base_path=args.base_path,cfg=cfg)
datasets = dsc.get_datasets(args.period)

minimizer_impl = IMinuitMinimizerImpl(cfg=cfg)


# create TimeLord
tl = TimeLord()

# apply source selection
df_srcs = pd.read_pickle(args.srcs_table_path)
idx = df_srcs.index
sources = []

for i, (idx, row) in enumerate(df_srcs.iterrows()):
    # Create external seyfert flux object + splinetable.
    src_name = df_srcs.index.values
    #print(src_name)
    ra_arr = df_srcs['ra'].values
    dec_arr = df_srcs['dec'].values
    weight_arr = df_srcs['weight'].values
    #print("source coords (ra,dec):", ra_arr, dec_arr,weight_arr)

    # Generate skyllh inputs.
    sources=[PointLikeSource(ra=src_ra, dec=src_dec, weight=src_weight)
           for (src_ra, src_dec, src_weight) in zip(np.deg2rad(ra_arr), np.deg2rad(dec_arr), weight_arr)]
    #print(sources)

#print("number of sources:", len(sources))

# Generate analysis instance.
ns_seed = 100
evt_sel = 'spatialbox'
optimize_delta_angle=10

#%%scalene
ana_stacking = trad_stacking_analysis.create_analysis(
    cfg=cfg,
    datasets=datasets,
    sources=sources,
    refplflux_E0=1e3,    
    optimize_delta_angle_deg=10,
#    energy_range=(30,200),
    ns_seed=100,
    ns_max=1e4,
    gamma_seed=2.5,
    gamma_min=2.0,
    gamma_max=4.0,
    fit_gamma=True,
    minimizer_impl=minimizer_impl
)


print(tl)

# Define a random state service.
rss_seed = args.rss_seed
rss = RandomStateService(rss_seed)

trials_outfile = os.path.join(args.outdir, args.outfile)

if args.unblind:
    # Run unblinded trial.
    with tl.task_timer('Running unblinded trial.') as tt:
        (TS, fitparam_dict, status) = ana.unblind(rss)

        # Create the structured array data type for the result array.
        result_dtype = [
            ('seed', int),
            ('flux_model', 'U32'),
            ('source_name', 'U32'),
            ('ts', np.float64)
        ] + [(key, np.float64) for key in fitparam_dict.keys()]

        result = np.empty((1,), dtype=result_dtype)
        result['seed'] = rss.seed
        result['flux_model'] = 'seyfert'
        result['source_name'] = 'stacking'
        result['ts'] = TS
        for (key, value) in fitparam_dict.items():
            result[key] = value

    os.makedirs(os.path.dirname(trials_outfile), exist_ok=True)
    np.save(trials_outfile, result)
else:
    # Run trials.
    with tl.task_timer('Running trials.') as tt:
        (seed,mean_ns,mean_ns_null,trial_data)=create_trial_data_file(
            ana=ana_stacking,
            ncpu=args.ncpus,
            rss=rss,
            pathfilename=trials_outfile,
            n_trials=args.n_trials,
            mean_n_sig=args.mean_ns,
            tl=tl)

print(tl)