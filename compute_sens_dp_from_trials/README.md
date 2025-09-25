-`run_stacking_sens_dp.sh` to compute the sensitivity and detect potential. 
default example case is dark matter stacking at 10TeV '/dm_model_stacking/trials/10TeV'

### one need to give its input parameters. 
1.1 path to anlysis case `dm_model_stacking/trials/$energy` or `powerlaw_stacking//trials/$energy` and 1.2 energy point `$energy`=`100GeV` or `1TeV` or `10TeV`.
 
-`plot_bkg.py` generates TS-distribution, `plot_nbias_corrected.py` generate ns-bias plot. One also need to setup the path to analysis case.

-`sensitivity_bands_plot.ipynb'` generate sensitivty of cross section for dm case and the flux sensitivity for powerlaw case.