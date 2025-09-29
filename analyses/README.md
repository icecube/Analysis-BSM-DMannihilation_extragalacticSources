## Folder structure

contains the scripts to run analysis example notebook(`example.ipynb`) and to generate trials (`submit` folder + `.py` script) on NPX cluster. Each folder corresponds to the an analysis scene.


please notice:
each analysis scene corresponds to a unique `analysis object` generate by SkyLLH
`from i3skyllh.analyses.trad_single_ps import analysis as trad_single_analysis` powerlaw catalog
`from i3skyllh.analyses.trad_stacked_ps import analysis as trad_stacking_analysis` powerlaw stacking
`from i3skyllh.analyses.trad_stacked_ps import analysis_dm as trad_stacking_analysis_dm` dm stacking
`from i3skyllh.analyses.trad_single_ps import analysis_dm as trad_single_analysis_dm` dm catalog

the `ana.mu2flux(event number)` function converts events number to the flux. It is also only valid for each analysis scene. One can run it either in the notebook or by `get_mu2flux.py` script. The example input could be 'python3 get_mu2flux.py -ns 28.37 -md dm -c stacking -e 10000'

An analysis objects in a notebook could use more than 4GB RAM (cause kernel dead/restarts).