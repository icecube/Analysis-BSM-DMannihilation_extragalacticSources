#/usr/bin/env python

import numpy as np

import os
USER = os.getenv('USER')

#ana = 'NorthernTracks_v005p01_KDE_PDF_v009_seyfert_model'
#base_path = '/data/ana/analyses/'

njobs_per_ns = 50
ntrials = 20
ncpu =  4
#period  = "IC86_2011-IC86_2021"

outdir = "/data/user/liruohan/dm_model_stacking/trials/"

mean_ns = np.linspace(281, 300, 20)

with open("trials_15.dag", "w") as fout:
    rss_seed = 0 
    for tns in mean_ns:
        for job_id in range(njobs_per_ns):
            rss_seed += 1
            tjob = "job_"+str(rss_seed) 

            outfile = 'mean_ns{}_trials{}_rss{:}.npy'.format(tns, ntrials, rss_seed)

            fout.write("JOB %s trials.submit\n" %(tjob))
            fout.write("VARS %s OUTDIR=\"%s\" OFILE=\"%s\" NTRIALS=\"%s\" MEANNS=\"%s\" RSS_SEED=\"%s\" CPU=\"%s\"\n" %(tjob, outdir, outfile, str(ntrials), str(tns), str(rss_seed), ncpu)) 