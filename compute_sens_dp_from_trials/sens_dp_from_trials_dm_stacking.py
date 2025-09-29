#!/usr/bin/env python

import sys
import glob
import os

import numpy as np
from scipy.optimize import root_scalar
from scipy.stats import norm
from scipy.signal import savgol_filter

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

from skyllh.core.utils.analysis import calculate_critical_ts_from_gamma
from photospline import (
    SplineTable,
    glam_fit,
    ndsparse,
)

import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--bkg_trials_dir", "-b", help="folder containing bkg trials",
        type=str,
        default='/data/user/liruohan/dm_model_stacking/trials/10TeV/bkg/')
parser.add_argument("--sig_trials_dir", "-s", help="folder containing sig trials",
        type=str,
        default='/data/user/liruohan/dm_model_stacking/trials/10TeV/sig/')
parser.add_argument("--energy", "-e", help="reference energy in GeV 100, 1000 or 10000",
        type=int,
        default=1000)
args=parser.parse_args()
print(args)


# load bkg trials.
print(os.path.join(args.bkg_trials_dir, '*.npy'))
bkg_trials = np.concatenate([np.load(f) for f in glob.glob(os.path.join(args.bkg_trials_dir, '*.npy'))])

#sensitivity
bkg_p_sens = 0.5
sig_p_sens = 0.9
#3 sigma dp
bkg_p_3dp = norm.sf(3,0,1)
sig_p_3dp = 0.5
# 5 sigma dp
bkg_p_5dp = norm.sf(5,0,1)
sig_p_5dp = 0.5
#
ts_critical_3dp = calculate_critical_ts_from_gamma(bkg_trials['ts'], bkg_p_3dp)
ts_critical_5dp = calculate_critical_ts_from_gamma(bkg_trials['ts'], bkg_p_5dp)
print('critical_ts',ts_critical_3dp,ts_critical_5dp)
# get an idea how many 0 trials. just out of curiosity
# xi = np.sum(bkg_trials['ts']>1.e-5) / len(bkg_trials)
# print("fraction of bkg trials with TS>0:",xi)
# print("")

# load signal trials (with poisson fluctuations)
# between mu=0 and mu=250.
# the trials are loaded into a dictionary with mu as key and the trials as value.
n_max=300
if args.energy == 10000:
    n_max = 100
if (args.energy == 1000) or (args.energy ==100):
    n_max = 300
number_of_signals = np.linspace(1, n_max, n_max)
sig_dir=args.sig_trials_dir
sig_trials_dict = {}

# injected values
percentage_above_bkg = []
percentage_above_critical_ts_3dp = []
percentage_above_critical_ts_5dp = []
mean_ns_inj = []
for ns in number_of_signals:
    if (args.energy == 10000) or (args.energy ==1000):
        fpath = sig_dir + 'mean_ns{}_trials20_rss*.npy'.format(str(ns))
    elif args.energy == 100:
        fpath = sig_dir + 'mean_ns_inj*_mean_ns_fit{}_trials20_rss*.npy'.format(str(round(ns)))
        #mean_ns_inj991.0_mean_ns_fit80_trials20_rss40.npy
    print(fpath)
    sig_trials = np.concatenate([np.load(f) for f in glob.glob(fpath)])
    
    mean_ns_inj.append(ns)
    # Calculate ns for percentile above background median
    percentage_above_bkg.append((
        np.sum(sig_trials['ts'] > np.median(bkg_trials['ts']))
    )/sig_trials.size)

    # Calculate ns for percentile above background median
    percentage_above_critical_ts_3dp.append((
        np.sum(sig_trials['ts'] > ts_critical_3dp)
    )/sig_trials.size)

        # Calculate ns for percentile above background median
    percentage_above_critical_ts_5dp.append((
        np.sum(sig_trials['ts'] > ts_critical_5dp)
    )/sig_trials.size)



# --- Step 2: Extrapolate the data ---
# We will fit a polynomial to the *smoothed* data to capture the trend.
# A degree of 2 (quadratic) is often a good starting point.
x_data=number_of_signals.copy()
y_sens=np.asarray(percentage_above_bkg)
y_3dp=np.asarray(percentage_above_critical_ts_3dp)
y_5dp=np.asarray(percentage_above_critical_ts_5dp)

x_1TeV=[0,200,300,500,1000]
y_1TeV=[0,26,60,150,300]
coef = np.polyfit(x_1TeV,y_1TeV,1)
poly1d_1TeV_corr = np.poly1d(coef)

degree = 3
coeffs_sens = np.polyfit(x_data, y_sens, degree)
poly_func_sens = np.poly1d(coeffs_sens)
coeffs_3dp = np.polyfit(x_data, y_3dp, degree)
poly_func_3dp = np.poly1d(coeffs_3dp)
coeffs_5dp = np.polyfit(x_data, y_5dp, degree)
poly_func_5dp = np.poly1d(coeffs_5dp)

def find_real(numbers):
    positive_reals = [c.real for c in numbers if c.imag == 0 and c.real > 0]
    return min(positive_reals)
    
print((poly_func_sens-0.9).roots)
print((poly_func_3dp-0.5).roots)
print((poly_func_5dp-0.5).roots)

mu_sens=find_real((poly_func_sens-0.9).roots)
if args.energy == 1000:
    mu_sens=poly1d_1TeV_corr(mu_sens)
print(f"sens at {mu_sens:.1f} events number")
print("")
mu_3dp=find_real((poly_func_3dp-0.5).roots)
if args.energy == 1000:
    mu_3dp=poly1d_1TeV_corr(mu_3dp)
print(f"3dp at {mu_3dp:.1f} events number")
print("")
mu_5dp=find_real((poly_func_5dp-0.5).roots)
if args.energy == 1000:
    mu_5dp=poly1d_1TeV_corr(mu_5dp)
print(f"5dp at {mu_5dp:.1f} events number")
print("")