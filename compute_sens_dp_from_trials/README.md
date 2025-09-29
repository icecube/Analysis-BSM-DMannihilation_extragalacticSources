### one need to give its input parameters.  

-`sens_dp_from_trials_stacking.py` and `sens_dp_from_trials_ps.py` compute the sensitivity and detect potential for dark matter analysis. One needs to manually setup the path to trials. One need to manually input the reference energy. The default example is for dm stacking 1TeV.

-`run_ps_sens_dp.sh` and `run_stacking_sens_dp.sh` compute the sensitivity and detect potential for powerlaw stacking and catalog analysis. No need for extra inputs.

-`send_dp_example.ipynb` is an auxiliary notebook illustrate the sens/dp calculation. It is not relevant to the outputs.

-`plot_bkg.py` plots TS-distribution, `plot_nbias.py` plot ns-bias plot. One also needs to manually setup the path to trials. The default example is for dm stacking 10TeV.
In `plot_bkg.py`, one needs to adjust or comment out `bkg_dirs = {}` according to analysis case;
In `plot_nbias.py`one needs to adjust th path to signal.

-`sensitivity_bands_plot.ipynb'` generate sensitivty of cross section for the dm case and the flux sensitivity for powerlaw case. It reads results from results.txt.

The other file are not directly related to results.