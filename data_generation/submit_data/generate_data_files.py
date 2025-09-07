#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT combo/V01-00-00

import sys
import os

import math
import pickle as pkl
import numpy as np
import tables

from icecube import dataio, dataclasses, astro
import time
import random

sys.path.append("/data/user/liruohan/prepare_Analysis-BSM-DMannihilation_extragalacticSources/data_generation/submit_data/")
#import TimeConverter as TConverter
import os.path

import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--year", "-yr", help="which year of data", type=str, default="12")
args=parser.parse_args()
print(args.year)


def extract_data(year, run_num, subrun_num, data_option, file_version, outputfile):
    
    data_var = dict()

    #dtype([('angErr', '<f8'), ('logE', '<f8'), ('dec', '<f8'), ('ra', '<f8'), 
    #('run', '<i8'), ('event', '<i8'), ('subevent', '<i8'), ('time', '<f8'),
    #('azi', '<f8'), ('zen', '<f8')])
    
    #General infos
    run_id_arr = np.array([])
    event_arr = np.array([])
    subevent_arr = np.array([])
    time_mjd = np.array([])
    accumulated_time_arr = np.array([])
    #Cuts
    L4noise_classifier = np.array([])
    L5nHit_DOMs = np.array([])
    L7osc_next_bool = np.array([])
    L7muon_classifier_all = np.array([])
    L7muon_classifier_up = np.array([])
    L7_ntop15 = np.array([])
    L7_nouter = np.array([])
    L7reco_vertex_z = np.array([])
    L7reco_vertex_rho36 = np.array([])
    L7reco_time = np.array([])   
    #Reco variables
    reco_zenith = np.array([])
    reco_zenith_err = np.array([])
    reco_azimuth = np.array([])
    reco_azimuth_err = np.array([])
    reco_psi = np.array([])
    reco_TotalEnergy = np.array([])
    reco_dec = np.array([])
    reco_ra = np.array([])
    reco_cascade_energy = np.array([])
    reco_track_energy = np.array([])

    NFiles = 0
    
    #print (run_num)
    
    #print(run_num[year][0],run_num[year][1]+1)
    run_range = np.arange(run_num[year][0],run_num[year][1]+1)
    print (run_range)
    
    for run in run_range:
        
        path = "/data/ana/LE/oscNext/pass2/data/level7_v02.00/IC86.{}/Run{}/".format(year,str(run).zfill(8))
        #print(path)
        #In burnsample, only keep run ending by 0
        if (data_option == "burnsample") and (str(run).endswith("0")==False):
            continue
            
        if os.path.isdir(path) == False:
            print ("Run{} does not exist".format(str(run).zfill(8)))
            continue

        else:
            print("Run{}".format(str(run).zfill(8)))
            for subrun in range(subrun_num):

                try:
                    #to_print = "Nothing"
                    filename = "oscNext_data_IC86.{}_level7_v02.00_pass2_Run{}_Subrun{}.hdf5".format(year,str(run).zfill(8),str(subrun).zfill(8))
                    #print('check name',filename)
                    hdf = tables.open_file(path+filename,'r')
                    #to_print = "File opened"
                    

                    ##Cuts##
                    L4noise_classifier = np.append(L4noise_classifier, hdf.root.L4_NoiseClassifier_ProbNu.cols.value[:])
                    L5nHit_DOMs = np.append(L5nHit_DOMs, hdf.root.L5_SANTA_DirectPulsesHitMultiplicity.cols.n_hit_doms[:])
                    L7osc_next_bool = np.append(L7osc_next_bool, hdf.root.L7_oscNext_bool.cols.value[:])
                    L7muon_classifier_all =np.append(L7muon_classifier_all, hdf.root.L7_MuonClassifier_FullSky_ProbNu.cols.value[:])
                    L7muon_classifier_up = np.append(L7muon_classifier_up, hdf.root.L7_MuonClassifier_Upgoing_ProbNu.cols.value[:])
                    L7_ntop15 = np.append(L7_ntop15, hdf.root.L7_CoincidentMuon_Variables.cols.n_top15[:])
                    L7_nouter = np.append(L7_nouter, hdf.root.L7_CoincidentMuon_Variables.cols.n_outer[:])
                    L7reco_vertex_z = np.append(L7reco_vertex_z, hdf.root.L7_reconstructed_vertex_z.cols.value[:])
                    L7reco_vertex_rho36 = np.append(L7reco_vertex_rho36, hdf.root.L7_reconstructed_vertex_rho36.cols.value[:])
                    L7reco_time = np.append(L7reco_time, hdf.root.L7_reconstructed_time.cols.value[:])
                    #print('cut done',filename)
                    ##Reconstructed Informations##
                    reco_zenith = np.append(reco_zenith, hdf.root.retro_crs_prefit__zenith.cols.median[:])
                    reco_zenith_err = np.append(reco_zenith_err,hdf.root.L7_retro_crs_prefit__zenith_sigma_tot.cols.value[:]) #(hdf.root.retro_crs_prefit__zenith.cols.upper_bound[:]-hdf.root.retro_crs_prefit__zenith.cols.lower_bound[:])/2)
                    reco_azimuth = np.append(reco_azimuth, hdf.root.retro_crs_prefit__azimuth.cols.median[:])
                    reco_azimuth_err = np.append(reco_azimuth_err,hdf.root.L7_retro_crs_prefit__azimuth_sigma_tot.cols.value[:])#(hdf.root.retro_crs_prefit__azimuth.cols.upper_bound[:]-hdf.root.retro_crs_prefit__azimuth.cols.lower_bound[:])/2)
                    reco_cascade_energy = np.append(reco_cascade_energy, hdf.root.retro_crs_prefit__cascade_energy.cols.median[:])
                    reco_track_energy = np.append(reco_track_energy, hdf.root.retro_crs_prefit__track_energy.cols.median[:])
                    #print('reco done',filename)
                                        ##General Informations##
                    run_id_arr = np.append(run_id_arr, hdf.root.I3EventHeader.cols.Run[:])
                    event_arr = np.append(event_arr, hdf.root.I3EventHeader.cols.Event[:])
                    subevent_arr = np.append(subevent_arr,hdf.root.I3EventHeader.cols.SubEvent[:])
                    time_mjd = np.append(time_mjd, hdf.root.I3EventHeader.cols.time_start_mjd[:])
                    accumulated_time_arr = np.append(accumulated_time_arr, hdf.root.L4_accumulated_time.cols.value[:])
                    #print('general info done',filename)
                    
                    NFiles+=1
                    hdf.close()
                    
                except Exception as e:
                    hdf.close()
                    #print(e)
                    continue
                    
    #print (NFiles)
    #print(len(reco_zenith),len(reco_azimuth),len(time_mjd))
    #print(reco_zenith,reco_azimuth,time_mjd)
    #Reco
    reco_ra, reco_dec = astro.dir_to_equa(reco_zenith, reco_azimuth, time_mjd)
    #Scramble in RA
    reco_ra = (2*math.pi)*np.random.rand(len(reco_dec))
    reco_psi = astro.angular_distance(reco_ra, reco_dec, np.radians(266.4167), np.radians(-29.0078))
    
    #General informations
    data_var["run"] = run_id_arr               
    data_var["MJD_time"] = time_mjd
    data_var["event"] = event_arr
    data_var["subevent"] = subevent_arr
    data_var["Accumulated_time"] = accumulated_time_arr
    
    #Reconstructed
    data_var["reco_TotalEnergy"] = reco_cascade_energy + reco_track_energy
    data_var["angErr"]= np.sqrt(reco_zenith_err**2 + reco_azimuth_err**2)
    data_var["logE"]= np.log10(reco_cascade_energy + reco_track_energy)
    data_var["zen"] = reco_zenith
    data_var["azi"] = reco_azimuth
    data_var["dec"] = reco_ra
    data_var["ra"] = reco_dec
    data_var["reco_psi"] = reco_psi
     
    #Cuts
    data_var["L4noise_classifier"]=L4noise_classifier
    data_var["L5nHit_DOMs"]=L5nHit_DOMs
    data_var["L7OscNext_bool"] = L7osc_next_bool
    data_var["L7muon_classifier_all"] = L7muon_classifier_all
    data_var["L7muon_classifier_up"] = L7muon_classifier_up
    data_var["L7_ntop15"] = L7_ntop15
    data_var["L7_nouter"] = L7_nouter
    data_var["L7reco_vertex_z"] = L7reco_vertex_z
    data_var["L7reco_vertex_rho36"] = L7reco_vertex_rho36
    data_var["L7reco_time"] = L7reco_time

    #print (NFiles)
    
    pkl.dump(data_var, open(outputfile, "wb"))
    return data_var


data_option = "test_sample" #"data" #
file_version = "v02.00"
years = [args.year]

runs={"12":[121269,122274],"13":[123691,124701],"14":[125370,126376]}
# runs={"12":[120028,121274],"13":[122282,123701],"14":[124551,125376],
#       "15":[126289,126949],"16":[127951,128518],"17":[129523,131059],"18":[131184,132046],
#       "19":[132765,133136],"20":[134064,134928]} #this is first round
#2011 and 2021 was move from previous data folder, but second run will still include them as double check
# runs={"11":[118550,119279],"12":[121275,122274],"13":[123702,124701],"14":[125377,126376],
#       "15":[126950,127949],"16":[128519,129518],"17":[131060,131259],"18":[132047,132846],
#       "19":[133137,134136],"20":[134929,135328],"21":[135251,136220]} # this is the second run
#15,16,17,19 need rerun, where only 17 need change the number
# runs={"11":[118550,119279],"12":[121275,122274],"13":[123702,124701],"14":[125377,126376],
#       "15":[126950,127949],"16":[128519,129518],"17":[129523,131059],"18":[132047,132846],
#       "19":[133137,134136],"20":[134929,135328],"21":[135251,136220]} # this is the third run 
# runs={"11":[118550,119879],"12":[120028,122274],"13":[122282,124701],"14":[124551,126376],
#       "15":[126289,127949],"16":[127951,129518],"17":[129523,131259],"18":[131184,132846],
#       "19":[132765,134136],"20":[134064,135328],"21":[135251,136220]} #this is the back up for most original runs,
subruns = 380 #365

for year in years:
    file_name = "/data/user/liruohan/prepare_Analysis-BSM-DMannihilation_extragalacticSources/data_generation/submit_data/pass2_Level7_v02.00_20{}.pkl".format(year)
    print ("-------------------------------------------")
    print ("### 20{} ###".format(year))
    print ("-------------------------------------------")
    print (file_name)
    data_dict = extract_data(year, runs, subruns, data_option, file_version, file_name)