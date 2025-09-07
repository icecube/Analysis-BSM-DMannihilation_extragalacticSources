## Folder structure


### submit_data ###
to run the data generation script on the NPX cluster, please first change the path in the .sh and .submit files to your own directory.
The script first generate .pkl file with filter conditions saved in it.
Then please run reprocess_data.ipynb notebook to apply the cuts and convert .pkl file into .npy (SkyLLH readable)

The script needs two days to run through all runs. To only check the validation, please use test_submit_script.ipynb and restrict run number about 100 to finish in 10 mins.

### OscNext_data ###
contain .npy data files (output from submit_data) that will be read by SkyLLH.
They are convert to IC86-2011-2021 dataset under `/data/user/liruohan/software/i3skyllh/i3skyllh/datasets/OscNext_v002p04.py`

### test_notebook ###
only contains the notebook to check NorthernTracks and OscNext energy/AngErr distribution, not directly related to the analysis

