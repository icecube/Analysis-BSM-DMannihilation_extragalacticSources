#!/usr/bin/env python

import sys
import glob
import os

# Add skyllh and i3skyllh projects to the PYTHONPATH
sys.path.insert(0, '/data/user/liruohan/software/skyllh')
sys.path.insert(0, '/data/user/liruohan/software/i3skyllh')
#sys.path.insert(0, '/home/liruohan/.local/lib/python3.7/site-packages')
#sys.path.insert(0, '/home/cbellenghi/.pyenv/versions/3.8.1/lib/python3.8/site-packages')

# Add missing python packages from cvmfs
#sys.path.insert(0, '/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/RHEL_7_x86_64/lib/python3.7/site-packages')
extra_path = "/cvmfs/icecube.opensciencegrid.org/users/tkontrimas/software/pip/python3.11/site-packages" # whatever individual directory it is
if extra_path not in sys.path:
    sys.path.append(extra_path)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

from skyllh.core.config import Config
from skyllh.core.random import RandomStateService
from skyllh.core.timing import TimeLord
from skyllh.core.source_model import PointLikeSource
from skyllh.core.minimizers.iminuit import IMinuitMinimizerImpl

# Pre-defined datasets
from i3skyllh.datasets import data_samples

#Pre-defined analysis
from i3skyllh.analyses.trad_single_ps import analysis as trad_single_analysis #powerlaw catalog
from i3skyllh.analyses.trad_stacked_ps import analysis as trad_stacking_analysis #powerlaw stacking
from i3skyllh.analyses.trad_single_ps import analysis_dm as trad_single_analysis_dm #dm catalog
from i3skyllh.analyses.trad_stacked_ps import analysis_dm as trad_stacking_analysis_dm #dm stacking

import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--ns", "-ns", help="sens or dp in signal events number", type=float, default="28.37")
parser.add_argument("--model", "-md", help="dark matter dm or powerlaw pl", type=str, default="dm")
parser.add_argument("--channel", "-ch", help="WW or bb channel, only need for dm case", type=str, default="WW")
parser.add_argument("--case", "-c", help="stacking search or catalog search", type=str, default="stacking")
parser.add_argument("--name", "-n", help="only need for catalog search", type=str, default="NGC4151")
parser.add_argument("--energy", "-e", help="reference energy point in GeV, 100,1000 or 10000", type=int, default=10000)
parser.add_argument("--rss_seed", "-rs", help="random_number_seed for mc generator", type=int, default=0)
parser.add_argument("--period", "-p", help="name of the data period to load", type=str, default="IC86_2011_2021")
parser.add_argument("--base_path", "-bp", help="base path to the folder containing datasets", type=str, default='/data/ana/analyses/')
parser.add_argument("--srcs_table_path", "-bstp", help="where to find sources pandas pkl table", type=str, default='/data/user/liruohan/dm_model_stacking/sources.pkl')
parser.add_argument("--ncpus", "-nc", help="number of CPUs to use", type=int, default=4)
parser.add_argument("--unblind", help="run the unblinded analysis", default=False, action="store_true")
args=parser.parse_args()


cfg=Config()
cfg.set_ncpu(4)
minimizer_impl = IMinuitMinimizerImpl(cfg=cfg)

dsc=None
#set other parameters according to input
if args.energy == 100:
    dsc = data_samples['OscNext_v002p04'].create_dataset_collection(base_path=args.base_path,cfg=cfg)
    print('use dataset OscNext')
elif ((args.energy) == 1000 or (args.energy==10000)):
    dsc = data_samples['NorthernTracks_v005p01'].create_dataset_collection(base_path=args.base_path,cfg=cfg)
    print('use dataset NorthernTracks')

datasets = dsc.get_datasets(args.period)


# apply source selection
sources = []
source = None
if args.case == 'stacking':
    df_srcs = pd.read_pickle(args.srcs_table_path)
    idx = df_srcs.index
    
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
    print('stacking analysis loaded sources')
elif args.case == 'catalog':
    single_src= most_significant_src[args.name]
    source = PointLikeSource(ra=np.deg2rad(single_src['ra']),dec=np.deg2rad(single_src['dec']))
    print('catalog analysis loaded source')

# create TimeLord
tl = TimeLord()

flux=None
if ((args.model == 'dm') & (args.case == 'stacking')):
    print('generate analysis object')
    ana=ana_stacking_dm = trad_stacking_analysis_dm.create_analysis(
    cfg=cfg,
    datasets=datasets,
    sources=sources,
    channel=args.channel,
    mass=args.energy,
    refplflux_E0=args.energy,   
#    energy_range=(0,500),
    compress_data=True,
    optimize_delta_angle_deg=10,
    ns_seed=100,
    ns_max=1e4,
    minimizer_impl=minimizer_impl
)
    flux= ana.sig_generator.mu2flux(args.ns)  
elif ((args.model == 'dm') & (args.case == 'catalog')):
    print('generate analysis object')
    ana = trad_single_analysis_dm.create_analysis(
    cfg=cfg,
    datasets=datasets,
    source=source,
    channel=args.channel,
    mass=args.energy,
    refplflux_E0=args.energy,  
    compress_data=True,
    optimize_delta_angle_deg=10,
    ns_seed=100,
    ns_max=1e4,
    minimizer_impl=minimizer_impl
)
    flux= ana.sig_generator.mu2flux(args.ns)  
elif ((args.model == 'pl') & (args.case == 'stacking')):
    print('generate analysis object')
    ana = trad_stacking_analysis.create_analysis(
    cfg=cfg,
    datasets=datasets,
    sources=sources,
#    bkg_event_rate_field_names=['astro', 'conv'],
    compress_data=True,
    refplflux_E0=args.energy,    
    optimize_delta_angle_deg=10,
    ns_seed=100,
    ns_max=1e4,
    gamma_seed=2.5,
    gamma_min=2.0,
    gamma_max=4.0,
    fit_gamma=True,
    minimizer_impl=minimizer_impl
)
    flux= ana.sig_generator.mu2flux(args.ns)  
elif ((args.model == 'pl') & (args.case == 'catalog')):
    print('generate analysis object')
    ana = trad_single_analysis.create_analysis(
    cfg=cfg,
    datasets=datasets,
    source=source,
    compress_data=True,
    refplflux_E0=args.energy,  
    optimize_delta_angle_deg=10,
    ns_seed=100,
    ns_max=1e4,
    gamma_min=2.0,
    gamma_max=4.0,
    fit_gamma=True,
    minimizer_impl=minimizer_impl
)
    flux= ana.sig_generator.mu2flux(args.ns)  
    
print('generate analysis object done')
print('flux is {}'.format(str(flux)))
print(tl)