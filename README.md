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
contains script needed to generate OscNext exp and simulation data. NorthernTrack data can be extracted from nu-source data storage.

### spectra
contains [PPPC](http://www.marcocirelli.net/PPPC4DMID.html) all spetra data and the script to extract needed neutrino spectra

### source list
Source list for two (DM and powerlaw) weighting scene.

### analyses
contains script to generate trials on NPX cluster for all 4 analyses cases. And the script for following sensitivity and detection potential calculations.
