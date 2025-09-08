#!/usr/bin/env python

import matplotlib as mpl
import matplotlib
import pylab as plt
import numpy as np
import glob
import os

#mpl.rcParams['font.family'] = 'sans-serif'
#mpl.rcParams['font.sans-serif'] = ['Verdana']

import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.WARNING)

import argparse

### change stepsize in xs according to your choice 
### of mean_ns during trials generation
xs = np.linspace(1, 40, 40)

import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--indir", "-i", help="folder containing signal trials",
        type=str,
        default='/data/user/liruohan/dm_model_stacking/trials/1TeV/sig/')
args=parser.parse_args()
print(args)

# import os
# for file in os.listdir("/mydir"):
#     if file.endswith(".txt"):
#         print(os.path.join("/mydir", file))

def get_data(indir, sn_str_list, var='ns'):
    
    data = []
    means = []

    for sn in sn_str_list:
        print(sn)
        vals = []
        inf = os.path.join(indir,  "mean_ns*"+"_trials20_"+"sn"+sn+"_rss*.npy")
        print(inf)
        infs = sorted(glob.glob(inf))
        #print(infs)
        for inf in infs:
            dat = np.load(inf)
            vals += dat[var].tolist()

        data.append(np.percentile(vals, [2.5, 16,50,84, 97.5]))
        means.append(np.mean(vals))
        #print(data)
    return(data, means)


fig = plt.figure(figsize=(8, 6), dpi=200)
ax = plt.axes()

sn_list = [str(int(x)) for x in xs]

fig = plt.figure(figsize=(8, 6), dpi=200)
ax = plt.axes()

indir = args.indir

data, means = get_data(indir, sn_list)
y_low2, y_low, y_mid, y_high, y_high2 = zip(*data)

ax.plot(xs, y_mid, color='tab:blue', linestyle='solid', linewidth=4, label=r"dm_model fit")
ax.fill_between(xs, y_low, y_high, color='tab:blue', alpha=0.3)
ax.fill_between(xs, xs-np.sqrt(xs), xs+np.sqrt(xs), color='tab:orange', alpha=0.2, label='poisson std error')
ax.plot(xs, y_high2, color='tab:blue', linestyle='dotted')
ax.plot(xs, y_low2, color='tab:blue', linestyle='dotted')
ax.plot(xs, xs, "k--")

ax.set_xlabel("$\mu^{inj}_{ns}$", position=(0., 1.), va='top', ha='right', fontsize=24)
ax.set_ylabel('$\mu^{fit}_{ns}$', position=(1., 0.), va='bottom', ha='right', fontsize=24)
ax.yaxis.set_label_coords(-0.12, 1.)
ax.xaxis.set_label_coords(1.0, -0.13)
ax.tick_params(axis='y', which='major', labelsize=18, pad=5)
ax.tick_params(axis='x', which='major', labelsize=18, pad=10)

ax.set_ylim([0.0, 40])
ax.set_xlim([0.0, 40])

leg = ax.legend(bbox_to_anchor=[0.99, 0.01], loc='lower right', prop={'size':12}, ncol=3, fancybox=True)
fig.subplots_adjust(bottom=0.15)

for axis in ['top','bottom','left','right']:
          ax.spines[axis].set_linewidth(1.5)
          ax.spines[axis].set_color('0.0')

ax.tick_params(axis='both', which='both', width=1.5, colors='0.0')
ax.yaxis.set_ticks_position('both')

plt.tight_layout()

plt.savefig("{}_ns_bias.png".format('dm_stacking'))