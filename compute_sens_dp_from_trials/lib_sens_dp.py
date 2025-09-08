import numpy as np
from skyllh.core.utils.analysis import calculate_critical_ts_from_gamma
from scipy.stats import binom, norm
from scipy.optimize import root_scalar
import photospline as psp
import matplotlib.pyplot as plt


def get_sens_dp(bkg_ts_tail_prob, signal_ts_tail_prob, bkg_trials, signal_trials_dict, floor=1.e-3, spline_regularization=1, plot=False):
    """
    Computes quantities like sensitivity and discovery potential from existing signal
    and background trials. Strategy: Calculation fraction of signal trials above
    critical TS value as function of injected mean number of ns (make sure the trials
    have poisson fluctuations turned on). Using spline-interpolation, it solves for
    the desired quantity.

    Args
        bkg_ts_tail_prob (float): desired fraction of bkg trials
                that have TS larger than critical value. This is used to determine the
                critical TS value. For sensitivity one would set bkg_ts_tail_prob = 0.5.
                For n-sigma discovery potential one would set bkg_ts_tail_prob = norm.sf(n, 0 ,1).
        signal_ts_tail_prob (float): desired fraction of signal trials that have TS larger
                than critical TS value. For sensitivity use signal_ts_tail_prob = 0.9.
                For n-sigma discovery potential set bkg_ts_tail_prob = 0.5.
        bkg_trials (numpy.ndarray): numpy record array that holds the skyllh trials (output)
        signal_trials_dict (dict): python dictionary that has mu_ns (injected value) as keys
                and corresponding trials (numpy record arrays) as values
        floor (float): smallest possible probability value to avoid problems with log(p) for p-> 0.
                Use 1/N_trials_per_mu as good rule of thumb. default = 1.e-3.
        spline_regularization (float): sets the regularization of the spline that interpolates
                the tail area probability as function of mu (and smooths over statistical errors).
                default = 1
        plot (bool): Whether to plot the interpolation function and corresponding data.
                Note: it's always a good idea to check the plots.
                default = False.

    Returns
        mu_sol (float): solution, i.e. the value of mu that satisfies bkg_ts_tail_prob
                and signal_ts_tail_prob
    """

    # Get critical TS value using skyllh convenience functions.
    # And treat cases where critical TS is below gamma-fit
    # to the test-statistic distribution.
    try:
        crit_ts = calculate_critical_ts_from_gamma(bkg_trials['ts'], bkg_ts_tail_prob)
        # noticed crit_ts = nan in some cases. switch to trials based calculation in such cases
        if np.isnan(crit_ts):
            print("crit_ts evaluated to nan. evaluate from trials.")
            crit_ts = np.percentile(bkg_trials['ts'], [(1.-bkg_ts_tail_prob) * 100])[0]
    except ValueError as error:
        if 'larger than the truncation threshold eta' in str(error):
            print("critical_ts value below range of gamma fit. will compute from trials directly.")
            crit_ts = np.percentile(bkg_trials['ts'], [(1.-bkg_ts_tail_prob) * 100])[0]
        else:
            print(error)
            raise error

    print("found critical TS:", crit_ts)

    # create a dictionary that holds the fraction of trials
    # above the critical TS as function of mu
    signal_prob_dict = {}
    for mu in signal_trials_dict.keys():
        sig_ts = signal_trials_dict[mu]['ts']
        signal_prob_dict[mu] = np.sum(sig_ts > crit_ts) / len(sig_ts)

    # prepare input for photospline interpolation
    xvals = list(signal_prob_dict.keys())
    yvals = [signal_prob_dict[mu] for mu in xvals]
    idx = np.argsort(xvals)
    xvals = np.array(xvals)[idx]
    yvals = np.array(yvals)[idx]
    max_mu = np.amax(xvals)
    # floor yvals
    yvals[yvals < floor] = floor
    yvals[yvals > 1-floor] = 1-floor
    # create logits since creating the spline on logit scale
    # has the beneficial property that it will always remain
    # within [0, 1] (as a probability should!)
    yvals = np.log(yvals / (1-yvals))

    # generate spline on logits
    # an introduce some padding
    knots = np.linspace(xvals[0]-1, xvals[-1]+1, max_mu*3)
    spline = psp.glam_fit(
        *psp.ndsparse.from_data(yvals, np.ones(len(yvals))),
        [xvals],
        [knots],
        [3], # order of spline
        [spline_regularization],
        [2]
    )

    # find root using scipy
    def obj(mu):
        y = spline.evaluate_simple([mu])
        # need to map from logits to probabilities
        y = 1./(1+np.exp(-y)) - signal_ts_tail_prob
        return y

    # search for solution within range of trials
    # i.e. we do not trust the extrapolation of the spline
    # beyond the trials range.
    mu_sol = root_scalar(obj, bracket = [np.amin(xvals), np.amax(xvals)])
	# check if converged
    if(mu_sol.converged):
        mu_sol = mu_sol.root
    else:
        raise RuntimeError('root finding did not converge. Are you sure the solution'
            'is within the range of injected ns?')

    if plot:
        # plot interpolation and solution to allow to visual inspection
        # first plot across entire range
        xvals = np.linspace(np.amin(list(signal_prob_dict.keys())), np.amax(list(signal_prob_dict.keys())), 500)
        yvals = spline.evaluate_simple([xvals])

        # transform spline result back to probabilities
        yvals_p = 1./(1 + np.exp(-yvals))

        plt.figure()
        plt.plot(signal_prob_dict.keys(), [signal_prob_dict[mu] for mu in signal_prob_dict.keys()],
                 lw=0., marker='o', markersize=4)
        plt.plot(xvals, yvals_p, 'r-')
        plt.axvline(mu_sol, color="red")
        plt.ylabel("P(TS>crit_TS)")
        plt.xlabel("signal mean")
        plt.savefig(f"p_mu_{bkg_ts_tail_prob:.2e}_sig_{signal_ts_tail_prob:.2f}.png", dpi=300)

        # second plot zoom around solution
        plt.figure()
        plt.plot(signal_prob_dict.keys(), [signal_prob_dict[mu] for mu in signal_prob_dict.keys()],
                 lw=0., marker='o', markersize=6)
        plt.plot(xvals, yvals_p, 'r-')
        plt.axvline(mu_sol, color="red")
        plt.ylabel("P(TS>crit_TS)")
        plt.xlabel("signal mean")
        plt.xlim([mu_sol-5, mu_sol+5])
        plt.ylim([signal_ts_tail_prob-0.1, signal_ts_tail_prob+0.1])
        plt.savefig(f"p_mu_zoom_bkg_{bkg_ts_tail_prob:.2e}_sig_{signal_ts_tail_prob:.2f}.png", dpi=300)

    return mu_sol