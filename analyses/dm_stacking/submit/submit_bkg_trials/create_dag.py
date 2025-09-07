#/usr/bin/env python

import os
USER = os.getenv('USER')

#ana = 'NorthernTracks_v005p01_KDE_PDF_v009_seyfert_model'
#base_path = '/data/ana/analyses/'

njobs = 1000
ntrials = 10
ncpu =  4
#period  = "IC86_2011-IC86_2021"

outdir = "/data/user/liruohan/dm_model_stacking/trials/"

mean_ns = 0.0

with open("trials.dag", "w") as fout:
    rss_seed = 0
    for job_id in range(njobs):
        rss_seed += 1
        tjob = "job_"+str(rss_seed) 

        outfile = 'mean_ns{}_trials{}_rss{:}.npy'.format(mean_ns, ntrials, rss_seed)

        fout.write("JOB %s trials.submit\n" %(tjob))
        fout.write("VARS %s OUTDIR=\"%s\" OFILE=\"%s\" NTRIALS=\"%s\" MEANNS=\"%s\" RSS_SEED=\"%s\" CPU=\"%s\"\n" %(tjob, outdir, outfile, str(ntrials), str(mean_ns), str(rss_seed), ncpu)) 