import numpy as np
from scipy.stats import norm
import glob
import sys, os
from lib_sens_dp import get_sens_dp

import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--bkg_trials_dir", "-b", help="folder containing bkg trials",
        type=str,
        default='/data/user/liruohan/dm_model_stacking/trials/1TeV/bkg/')
parser.add_argument("--sig_trials_dir", "-s", help="folder containing sig trials",
        type=str,
        default='/data/user/liruohan/dm_model_stacking/trials/1TeV/sig/')
args=parser.parse_args()
print(args)


# load bkg trials.
print(os.path.join(args.bkg_trials_dir, '*.npy'))
bkg_trials = np.concatenate([np.load(f) for f in glob.glob(os.path.join(args.bkg_trials_dir, '*.npy'))])

# get an idea how many 0 trials. just out of curiosity
xi = np.sum(bkg_trials['ts']>1.e-5) / len(bkg_trials)
print("fraction of bkg trials with TS>0:",xi)
print("")

# load signal trials (with poisson fluctuations)
# between mu=0 and mu=250.
# the trials are loaded into a dictionary with mu as key and the trials as value.
sig_trials_dict = {}

# injected values
#mus = np.linspace(0, 250, 251)
mus = np.linspace(1, 100, 100)
for mu in mus:
#    mu =int(mu) #could delete this line according to format
    print(str(mu))
    fpath = args.sig_trials_dir + '/mean_ns{}_trials20_rss*.npy'.format(str(mu))
    print(mu,fpath)
    #outfile = 'mean_ns{}_rss{:}.npy'.format(tns, rss_seed)
    sig_trials_dict[int(mu)] = np.concatenate([np.load(f) for f in glob.glob(fpath)])


# sensitivity
bkg_p = 0.5
sig_p = 0.9
mu = get_sens_dp(bkg_p, sig_p, bkg_trials, sig_trials_dict, floor=1./len(sig_trials_dict[np.amin(mus)]), plot=True)
print(f"sens at {mu:.1f} events")
print("")

# 3 sigma dp
bkg_p = norm.sf(3,0,1)
sig_p = 0.5
mu = get_sens_dp(bkg_p, sig_p, bkg_trials, sig_trials_dict, floor=1./len(sig_trials_dict[np.amin(mus)]), plot=True)
print(f"3sigma dp at {mu:.1f} events")
print("")

# 5 sigma dp
bkg_p = norm.sf(5,0,1)
sig_p = 0.5
mu = get_sens_dp(bkg_p, sig_p, bkg_trials, sig_trials_dict, floor=1./len(sig_trials_dict[np.amin(mus)]), plot=True)
print(f"5sigma dp at {mu:.1f} events")
print("")