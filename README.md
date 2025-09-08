# Analysis-BSM-DMannihilation_extragalacticSources

[Analysis Wiki Page](https://wiki.icecube.wisc.edu/index.php/DM_annihilation_from_SMBH_as_point_source)

Environment: `/cvmfs/icecube.opensciencegrid.org/py3-v4.3.0/setup.sh`

Dataset: `/data/ana/analyses/northern_tracks/version-005-p02` and (add OscNext path)

Software: [SkyLLH](https://github.com/icecube/skyllh/tree/dm) and [i3SkyLLH](https://github.com/icecube/i3skyllh/tree/https://github.com/icecube/skyllh/tree/dm)

In the executable script add `skyllh` and `i3skyllh` software and `photospline` package to the `PYTHONPATH`:

```
export PYTHONPATH="/data/user/liruohan/software/skyllh":$PYTHONPATH
export PYTHONPATH="/data/user/liruohan/software/i3skyllh":$PYTHONPATH
export PYTHONPATH=/cvmfs/icecube.opensciencegrid.org/users/tkontrimas/software/photospline_v2.2.0/:$PYTHONPATH
```

## Folder structure

### data_generation
the folder contains scripts that needed to extract OscNext experimental data and apply cuts. NorthernTrack data can be extracted from nu-source data storage by SkyLLH. The simulation data is not needed because we scrambling the exp data to generated the background.

### spectra
The spectra are taken from [PPPC](http://www.marcocirelli.net/PPPC4DMID.html) . The folder contains all spetra data and the script to extract needed channel and neutrino spectra

### source list
Source list for two (DM and powerlaw) weighting scene.

### analyses
contains script to generate trials on NPX cluster for all 4 analyses cases. And the scripts for sensitivity and detection potential calculations.

### compute_sens_dp_from_trials
contains script check bkg-TS distribution, ns-bias, and compute sensitivity and detect potential.
