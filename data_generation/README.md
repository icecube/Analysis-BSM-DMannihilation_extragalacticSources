## Folder structure


### submit_data ###
To run the `generate_data_file.py` and  script on the NPX cluster, please first change the save directory in `.py`,`.sh`,`.submit` to your own)
The scripts first generate .pkl file with L4-7 filter conditions saved in it.
The`reprocess_data.ipynb` notebook apply the cuts and convert .pkl file into .npy (SkyLLH readable)

The script needs two days to run through all runs. For validation check only, please use `test_submit_script.ipynb` and restrict run number about 100 to finish in 10 mins.

### OscNext_data ###
contain .npy data files (output from submit_data) that will be read by SkyLLH.

### test_notebook ###
only contains the notebook to check NorthernTracks and OscNext energy/AngErr distribution. They are not directly related to the analysis.

