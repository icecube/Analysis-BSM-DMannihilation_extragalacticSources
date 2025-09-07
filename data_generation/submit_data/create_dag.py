#/usr/bin/env python

import os
USER = os.getenv('USER')


outdir = "/data/user/liruohan/prepare_Analysis-BSM-DMannihilation_extragalacticSources/data_generation/submit_data/"

years = ["12","13","14","15","16","17","18","19","20","21"]


for i,year in enumerate(years):
    with open("run_data{}.dag".format(year), "w") as fout:
        tjob = "job_"+str(i) 
        fout.write("JOB %s run_data.submit\n" %(tjob))
        fout.write("VARS %s YEAR=\"%s\"\n" %(tjob, year))