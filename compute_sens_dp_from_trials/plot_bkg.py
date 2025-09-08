import argparse
import glob
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from iminuit import minimize
from scipy.stats import chi2, gamma

#mpl.rcParams["font.family"] = "sans-serif"
#mpl.rcParams["font.sans-serif"] = ["Verdana"]


# Define helper functions.
def get_confidence_interval(counts):
    """Calculate confidence regions.
    For high statistics: approximate as symmetric gaussian with sqrt(N).
    For low statistics: use asymmetric feldman cousins.
    """
    # Define feldman cousins confidence intervals for the low statistics region
    feldman_cousins = {
        0: (0.00, 1.29),
        1: (0.37, 2.75),
        2: (0.74, 4.25),
        3: (1.10, 5.30),
        4: (2.34, 6.78),
        5: (2.75, 7.81),
        6: (3.82, 9.28),
        7: (4.25, 10.30),
        8: (5.30, 11.32),
        9: (6.33, 12.79),
        10: (6.78, 13.81),
        11: (7.81, 14.82),
        12: (8.83, 16.29),
        13: (9.28, 17.30),
        14: (10.30, 18.32),
        15: (11.32, 19.32),
        16: (12.33, 20.80),
        17: (12.79, 21.81),
        18: (13.81, 22.82),
        19: (14.82, 23.82),
        20: (15.83, 25.30),
    }

    if counts > 20:
        err = np.sqrt(counts)
        return (err, err)
    else:
        fc = feldman_cousins[counts]
        return (counts - fc[0], fc[1] - counts)


def truncated_gamma_logpdf(a, scale, eta, ts_above_eta, N_above_eta):
    """Calculates the -log(likelihood) of a sample of random numbers
    generated from a gamma pdf truncated from below at x=eta.
    Parameters
    ----------
    a : float
        Shape parameter.
    scale : float
        Scale parameter.
    eta : float
        Test-statistic value at which the gamma function is truncated
        from below.
    ts_above_eta : (n_trials,)-shaped 1D ndarray
        The ndarray holding the test-statistic values falling in
        the truncated gamma pdf.
    N_above_eta : int
        Number of test-statistic values falling in the truncated
        gamma pdf.

    Returns
    -------
    -logl : float
    """
    c0 = 1.0 - gamma.cdf(eta, a=a, scale=scale)
    c0 = 1.0 / c0
    logl = N_above_eta * np.log(c0) + np.sum(
        gamma.logpdf(ts_above_eta, a=a, scale=scale)
    )
    return -logl


def get_expectation_from_fit(
    edge_low, edge_high, pars, eta, fract_above_eta, Ntot
):
    """Calculates the expected bin content given the fit parameters from a
    truncated gamma pdf. Fudge the normalization to match above truncation point
    eta.
    """
    norm = fract_above_eta / (1 - gamma.cdf(eta, a=pars[0], scale=pars[1]))
    cdf1 = gamma.cdf(edge_high, a=pars[0], scale=pars[1])
    cdf2 = gamma.cdf(edge_low, a=pars[0], scale=pars[1])
    prob = cdf1 - cdf2
    prob *= norm
    return prob * Ntot


def get_expectation_from_chi2(edge_low, edge_high, Ntot):
    cdf1 = chi2.cdf(edge_high, 1)
    cdf2 = chi2.cdf(edge_low, 1)
    prob = cdf1 - cdf2
    return prob * Ntot * 0.5


def do_fit(ts, label, outdir, eta=1.0, xi=1.0e-2):
    Ntot = len(ts)
    ts_eta = ts[ts > eta]
    N_prime = len(ts_eta)
    fract_above_eta = N_prime / Ntot

    obj = lambda x: truncated_gamma_logpdf(
        x[0], x[1], eta=eta, ts_above_eta=ts_eta, N_above_eta=N_prime
    )
    x0 = [0.75, 1.8]  # Initial values of function parameters.
    bounds = [[0.1, 10], [0.1, 10]]  # Ranges for the minimization fitter.
    r = minimize(obj, x0, bounds=bounds)
    pars = r.x

    # Generate histogram for analysis and visualiztion
    counts_obs, bins = np.histogram(
        ts, bins=np.linspace(0.0, 20, 100), density=False
    )
    binsc = 0.5 * (bins[1:] + bins[:-1])

    # Calculate uncertainties.
    errs = [get_confidence_interval(count) for count in counts_obs]
    lerr, herr = zip(*errs)

    # Calculate expectations.
    counts_exp = np.asarray(
        [
            get_expectation_from_fit(le, he, pars, eta, fract_above_eta, Ntot)
            for le, he in zip(bins[:-1], bins[1:])
        ]
    )
    counts_exp2 = np.asarray(
        [
            get_expectation_from_chi2(le, he, Ntot)
            for le, he in zip(bins[:-1], bins[1:])
        ]
    )

    # Calculate residuals.
    diff = counts_obs - counts_exp
    idxp = diff > 0  # more exp than data
    idxm = np.logical_not(idxp)  # less exp than data

    # Calculate standardized residual depending on direction.
    diff[idxp] = diff[idxp] / np.asarray(lerr)[idxp]
    diff[idxm] = diff[idxm] / np.asarray(herr)[idxm]

    fig, (ax, ax2) = plt.subplots(
        nrows=2,
        ncols=1,
        sharex=True,
        sharey=False,
        figsize=(8, 6),
        gridspec_kw={"height_ratios": [4, 1]},
    )

    # Plot data/fit.
    ax.hist(
        binsc,
        bins=bins,
        weights=counts_obs,
        color="tab:blue",
        histtype="step",
        label=label,
    )
    ax.hist(binsc, bins=bins, weights=counts_obs, color="tab:blue", alpha=0.2)
    ax.errorbar(
        binsc,
        counts_obs,
        yerr=[lerr, herr],
        linewidth=0,
        color="tab:blue",
        elinewidth=1,
    )
    ax.plot(
        binsc,
        counts_exp,
        "k--",
        linewidth=1,
        label="gamma-fit") #($\eta>{},\,\\xi>{}$)")
    ax.plot(
        binsc, counts_exp2, "r--", linewidth=1, label=r"$0.5 \cdot \chi^2(ndof=1)$"
    )

    # Plot residuals.
    ax2.plot([0, 35], [1, 1], "k--", linewidth=0.5)
    ax2.plot([0, 35], [-1, -1], "k--", linewidth=0.5)
    ax2.plot([0, 35], [0, 0], "k-", linewidth=0.5)
    ax2.plot([0, 35], [2, 2], color="k", linestyle="dotted", linewidth=0.5)
    ax2.plot([0, 35], [-2, -2], color="k", linestyle="dotted", linewidth=0.5)

    ax2.plot(binsc, diff, color="tab:blue", marker=".", linewidth=0)

    # Plotting properties.
    ax2.set_xlabel("test-statistic", fontsize=18)
    ax2.set_ylabel(
        r"$\Delta / \sigma$",
        position=(0.0, 0.5),
        va="top",
        ha="right",
        fontsize=13,
    )
    ax.set_ylabel("counts", fontsize=18)
    ax.set_ylim([0.1, 5.0e5])
    ax.set_xlim([0, 20])
    ax.set_yscale("log")
    ax.tick_params(
        axis="both", which="both", width=1.5, colors="0.0", labelsize=16
    )
    ax.yaxis.set_ticks_position("both")
    ax2.set_ylim([-3.5, 3.5])
    ax2.tick_params(axis="y", which="major", labelsize=10)
    ax2.tick_params(axis="x", which="major", labelsize=18)
    ax2.yaxis.set_ticks_position("both")
    ax2.yaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(1.5)
        ax.spines[axis].set_color("0.0")
        ax2.spines[axis].set_linewidth(1.5)
        ax2.spines[axis].set_color("0.0")
    ax2.yaxis.set_label_coords(-0.1, 0.7)
    ax2.set_yticks([-2.0, 0.0, 2.0])
    ax.legend(
        bbox_to_anchor=[0.99, 0.99],
        loc="upper right",
        prop={"size": 12},
        ncol=1,
        fancybox=True,
    )
    plt.subplots_adjust(hspace=0)

    # Use label as a figure name.
    fig_name = "bkg_TS_" + label + ".png"

    os.makedirs(outdir, exist_ok=True)

    plt.savefig(os.path.join(outdir, fig_name), dpi=500)
    plt.close()

    return pars, eta, fract_above_eta


parser = argparse.ArgumentParser()
parser.add_argument("--analysis", help="select analysis type option of: 'catalog_model', 'catalog_powerlaw', 'stacking_model', 'all'", type=str, default='all')
parser.add_argument("--bass_table_path", "-tp", help="path to bass tables", type=str, default='/data/user/liruohan/powerlaw/')
parser.add_argument("--outdir", "-o", help="output folder", type=str, default='./figures/')

# Custom bkg_trials_dir option.
parser.add_argument("--custom_bkg_trials_dir", help="path to custom background trials. If set, the `--analysis` option is ignored", type=str, default=None)
parser.add_argument("--custom_label", help="label of background trials test-statistic histogram", type=str, default="")

args = parser.parse_args()
print(args)

if args.custom_bkg_trials_dir is not None:
    # Ignore analysis option
    infiles = glob.glob(os.path.join(args.custom_bkg_trials_dir, "*.npy"))
    trials = np.concatenate([np.load(infile) for infile in infiles])
    ts = np.asarray(trials["ts"])

    # Fit truncated gamma pdf and plot.
    pars, eta, fract_above_eta = do_fit(ts, args.custom_label, args.outdir)
    print('Plotted {}'.format(args.custom_label))

else:
    # Use analysis option.
    # Background trials directories.
    bkg_dirs = {
        #'catalog_model': '/data/user/liruohan/dm_model/trials/bkg',
        #'catalog_powerlaw': '/data/user/liruohan/powerlaw/trials/bkg/',
        #'stacking_powerlaw': '/data/user/liruohan/powerlaw_stacking/trials/bkg/',
        'stacking_model': '/data/user/liruohan/dm_model_stacking/trials/1TeV/bkg/'
    }

    if args.analysis == 'all':
        analyses = ['stacking_model']
    else:
        analyses = [args.analysis]

    for analysis in analyses:
        if 'catalog' in analysis:
            # Plot all sources.
            # read bass table and keep source properties and flux model
            df_srcs = pd.read_pickle(
                os.path.join(
                    args.bass_table_path, "single_sources.pkl"
                )
            )
            
            source_name = df_srcs.index.values

            for (idx, name) in enumerate(source_name):

                if 'model' in analysis:
                    indir = bkg_dirs['catalog_model']
                    label = name + '_model_fit'
                elif 'powerlaw' in analysis:
                    indir = bkg_dirs['catalog_powerlaw']
                    label = name + '_powerlaw_fit'
                else:
                    raise ValueError("Provided analysis option '{}' is not supported.".format(args.analysis))
                print(indir)
                infiles = glob.glob(os.path.join(indir, '{}*.npy'.format(source_name)))
                trials = np.concatenate([np.load(infile) for infile in infiles])
                ts = np.asarray(trials["ts"])

                # Fit truncated gamma pdf and plot.
                pars, eta, fract_above_eta = do_fit(ts, label, args.outdir)
                print(r"Plotted {analysis} {source_name}.")

        elif 'stacking' in analysis:
            # include_NGC1068 part.
            print(analysis)
            indir = bkg_dirs['stacking_model']
            label = 'dm stacking'

            infiles = glob.glob(os.path.join(indir, "*.npy"))
            trials = np.concatenate([np.load(infile) for infile in infiles])
            ts = np.asarray(trials["ts"])

            # Fit truncated gamma pdf and plot.
            pars, eta, fract_above_eta = do_fit(ts, label, args.outdir)
            print(r"Plotted {analysis}_include_NGC1068.")

#             # exclude_NGC1068 part.
#             indir = bkg_dirs['stacking_model_exclude_NGC1068']
#             label = 'seyfert stacking (N=27, no NGC1068)'

#             infiles = glob.glob(os.path.join(indir, "*.npy"))
#             trials = np.concatenate([np.load(infile) for infile in infiles])
#             ts = np.asarray(trials["ts"])

            # Fit truncated gamma pdf and plot.
#             pars, eta, fract_above_eta = do_fit(ts, label, args.outdir)
#             print("Plotted {}_exclude_NGC1068.".format(analysis))
        else:
            raise ValueError("Provided analysis option '{}' is not supported.".format(args.analysis))