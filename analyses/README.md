## Folder structure

contains the scripts to run analysis example notebook(`example.ipynb`) or to generate trials (`submit` folder + `.py` script) on NPX cluster. Each folder corresponds to the an analysis scene.


please notice:
each analysis scene corresponds to a unique `analysis object` generate by SkyLLH
`from i3skyllh.analyses.trad_single_ps import analysis as trad_single_analysis` powerlaw catalog
`from i3skyllh.analyses.trad_stacked_ps import analysis as trad_stacking_analysis` powerlaw stacking
`from i3skyllh.analyses.trad_stacked_ps import analysis_dm as trad_stacking_analysis_dm` dm stacking
`from i3skyllh.analyses.trad_single_ps import analysis_dm as trad_single_analysis_dm` dm catalog

the function `ana.mu2flux(event number)` function also only valid for each analysis scene.

More than two analysis objects could use more than 4GB RAM (cause kernel dead/restart).