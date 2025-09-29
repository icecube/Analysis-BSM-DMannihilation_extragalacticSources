#/usr/bin/env python
import pandas as pd
import numpy as np
import os

#USER = os.getenv('USER')

#ana = 'NorthernTracks_v005p01_KDE_PDF_v007'
#base_path = '/data/ana/analyses/'

njobs_per_ns = 50
ntrials = 20
ncpu =  4
period  = "IC86_2011-IC86_2021"

outdir = "/data/user/liruohan/powerlaw/trials/"

df_srcs = pd.read_pickle('/data/user/liruohan/powerlaw/single_sources.pkl')

idx1 = np.logical_and(df_srcs['dec'] > -5, df_srcs['ra'] > -360)
#idx2 = df_bass['TYPE_105'].str.contains('Sy')
#idx3 = df_bass['F2-10-intr'] > 27
df_srcs = df_srcs[idx1].sort_values(by='dec', ascending=False).copy(deep=True)

gamma = 2.0
names = df_srcs.index.values 
mean_ns = np.arange(1,50,1)

for i,src in enumerate(names):
    with open("trials_{}.dag".format(i), "w") as fout:
        rss_seed = 0 
        for tns in mean_ns:
            for job_id in range(njobs_per_ns):
                rss_seed += 1
                tjob = "job"+str(rss_seed) 

                outfile = '{}_mean_ns{}_trials{}_rss{:}.npy'.format(src, tns, ntrials, rss_seed)

                fout.write("JOB %s trials.submit\n" %(tjob))
                fout.write("VARS %s source=\"%s\" OUTDIR=\"%s\" OFILE=\"%s\" GAMMA=\"%s\" CPU=\"%s\" NTRIALS=\"%s\" MEANNS=\"%s\" RSS_SEED=\"%s\"\n" %(tjob,src,outdir,outfile, gamma, str(ncpu), str(ntrials), str(tns), str(rss_seed)))