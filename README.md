# Analysis-BSM-DMannihilation_extragalacticSources

[Analysis Wiki Page](https://wiki.icecube.wisc.edu/index.php/DM_annihilation_from_SMBH_as_point_source)

Environment: `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`

Dataset: `/data/ana/analyses/northern_tracks/version-005-p02` and OscNext High Stat data

Software: SkyLLH dm branch  [SkyLLH](https://github.com/icecube/skyllh/tree/dm) and [i3SkyLLH](https://github.com/icecube/i3skyllh/tree/https://github.com/icecube/skyllh/tree/dm)

In the executable script add `skyllh` and `i3skyllh` software and `photospline` package to the `PYTHONPATH`:

```
export PYTHONPATH="/data/user/liruohan/software/skyllh":$PYTHONPATH
export PYTHONPATH="/data/user/liruohan/software/i3skyllh":$PYTHONPATH
export PYTHONPATH=/cvmfs/icecube.opensciencegrid.org/users/tkontrimas/software/photospline_v2.2.0/:$PYTHONPATH
```

## Folder structure

### data_generation
the folder contains scripts to extract OscNext experimental data and apply cuts. NorthernTrack data can be read from nu-source data storage by SkyLLH. The simulation data is not needed because we scrambling the exp data to generated the background.

### spectra
The spectra .dat are taken from [PPPC](http://www.marcocirelli.net/PPPC4DMID.html) . The folder contains all spetra data and the notebook to: 1. extract neutrino spectra at given energy and channel; 2.  generate spectra spline and smoothen the tail. The outpout .fits spline will read by SkyLLH. 

### source list
Source list for two (DM and powerlaw) weighting scene.

### analyses
contains scripts to generate trials on NPX cluster for all 4 analyses cases.

### compute_sens_dp_from_trials
contains scripts check bkg-TS distribution, ns-bias, and compute sensitivity and detect potential.
