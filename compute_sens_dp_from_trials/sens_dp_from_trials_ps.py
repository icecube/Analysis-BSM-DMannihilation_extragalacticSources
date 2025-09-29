import numpy as np
from scipy.stats import norm
import glob
import sys
import pandas as pd
import argparse

from lib_sens_dp import get_sens_dp

parser=argparse.ArgumentParser()
parser.add_argument("--analysis",  help="powerlaw or model", type=str, default='powerlaw')
parser.add_argument("--bass_table_path", "-tp", help="path to bass tables", type=str, default='/data/user/liruohan/powerlaw/single_sources.pkl')

analysis = parser.parse_known_args()[0].analysis

if analysis == 'model':
    option_parser = argparse.ArgumentParser()
    option_parser.add_argument("--trials_location_bkg", "-pbg" ,help="path to background trials", type=str, default="/data/ana/analyses/NuSources/2022_Seyfert_Galaxies_Analysis/trials/catalog/model/bkg/")
    option_parser.add_argument("--trials_location_signal","-psig",  help="path to signal trials", type=str, default="/data/ana/analyses/NuSources/2022_Seyfert_Galaxies_Analysis/trials/catalog/model/signal/")

elif analysis == 'powerlaw':
    option_parser = argparse.ArgumentParser()
    option_parser.add_argument("--gamma",  help="spectral index for powerlaw analysis", type=float, default=2.0)
    option_parser.add_argument("--trials_location_bkg", "-pbg" ,help="path to background trials", type=str, default="/data/user/liruohan/powerlaw/trials/bkg/")
    option_parser.add_argument("--trials_location_signal","-psig",  help="path to signal trials", type=str, default="/data/user/liruohan/powerlaw/trials/sig/")

args = option_parser.parse_known_args()[0]

print("Computing for {} analysis".format(analysis))

df_srcs = pd.read_pickle('/data/user/liruohan/powerlaw/single_sources.pkl') #args.bass_table_path

# idx1 = np.logical_and(df_bass['DECdeg'] > -5, df_bass['DECdeg'] < 81) # we should move this to -3deg
# idx2 = df_bass['TYPE_105'].str.contains('Sy')
# idx3 = df_bass['F2-10-intr'] > 27
# df_srcs = df_bass[idx1 & idx2 & idx3].sort_values(by='F2-10-intr', ascending=False).copy(deep=True)

names = df_srcs.index.values
decs = df_srcs["dec"].values
ras = df_srcs["ra"].values

if analysis == "model":
    fout = open("sens_dp_catalog_{}_analysis.out".format(analysis),"w")

elif analysis == "powerlaw":
    fout = open("sens_dp_catalog_{}_analysis_gamma_{}.out".format(analysis,args.gamma),"w")

fout.write("#name\tdec\tra\tsens\t3sigma dp\t 5sigma dp\n")

base_path_bkg    = args.trials_location_bkg
base_path_signal = args.trials_location_signal

if analysis == "model":
    mu_range = np.linspace(0, 60, 61)
elif analysis == "powerlaw":
    if args.gamma == 2.0:
        mu_range = np.linspace(0, 20, 21)
    elif args.gamma == 3.0:
        mu_range = np.linspace(0, 20, 21)

for l,source in enumerate(names):
    #print(source)
    #print(base_path_bkg)
    fout.write(source+"\t"+str(decs[l])+'\t'+str(ras[l])+'\t')
    bkg_trials_files = glob.glob(base_path_bkg+'/{}*.npy'.format(source))
    # make sure trials exist
    if len(bkg_trials_files) > 0:
        bkg_trials = np.concatenate([np.load(f) for f in bkg_trials_files])
    else:
        # trials don't exist. move to next source.
        print(f"Can't find background trials for source {source}. Skipping!")
        continue

    # get an idea how many 0 trials. just out of curiosity
    xi = np.sum(bkg_trials['ts']==0.0) / len(bkg_trials)
    print("fraction of bkg trials with TS=0:",xi)
    print("")

    sig_trials_dict = {}

    for mu in mu_range:
        if mu ==0:
            fpath = base_path_bkg + '/{}*.npy'.format(source)
        else:
            if analysis ==  "model":
                fpath = base_path_signal + '/{}_mean_ns{}_rss*.npy'.format(source,int(mu))
            elif analysis == "powerlaw":
                fpath = base_path_signal + '/{}/{}_mean_ns{}_trials100_rss*.npy'.format(source,source,int(mu))
        #print(fpath)
        sig_trials_dict[int(mu)] = np.concatenate([np.load(f) for f in glob.glob(fpath)])

    # sensitivity
    bkg_p = 0.5
    sig_p = 0.9
    mu = get_sens_dp(bkg_p, sig_p, bkg_trials, sig_trials_dict, floor=1./len(sig_trials_dict[np.amin(mu_range)]), plot=False)
    print(f"sens at {mu:.1f} events")
    print("")
    fout.write(str(mu)+"\t")

    # 3 sigma dp
    bkg_p = norm.sf(3,0,1)
    sig_p = 0.5
    mu = get_sens_dp(bkg_p, sig_p, bkg_trials, sig_trials_dict, floor=1./len(sig_trials_dict[np.amin(mu_range)]))
    print(f"3sigma dp at {mu:.1f} events")
    print("")
    fout.write(str(mu)+"\t")


    # 5 sigma dp
    bkg_p = norm.sf(5.0,0,1)
    sig_p = 0.5
    mu = get_sens_dp(bkg_p, sig_p, bkg_trials, sig_trials_dict, floor=1./len(sig_trials_dict[np.amin(mu_range)]))
    print(f"5sigma dp at {mu:.1f} events")
    print("")
    fout.write(str(mu)+"\n")

fout.close()