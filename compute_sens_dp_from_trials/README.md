### one need to give its input parameters. 
1.1 path to anlysis case `dm_model_stacking/trials/$energy` or `powerlaw_stacking//trials/$energy` and 1.2 energy point `$energy`=`100GeV` or `1TeV` or `10TeV`.
 

-`run_stacking_sens_dp.sh` and `run_ps_sens_dp.sh` to compute the sensitivity and detect potential.




-`plot_bkg.py` plots TS-distribution, `plot_nbias_corrected.py` plot ns-bias plot. One also need manually setup the path to tirals. 

-`sensitivity_bands_plot.ipynb'` generate sensitivty of cross section for dm case and the flux sensitivity for powerlaw case. It reads results from results.txt.