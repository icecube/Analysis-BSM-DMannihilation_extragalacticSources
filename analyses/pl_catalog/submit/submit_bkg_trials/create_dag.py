#/usr/bin/env python
import pandas as pd
import numpy as np
import os

#ana = 'trials'
#base_path = '/data/user/liruohan/powerlaw/'

njobs_per_ns = 100
ntrials = 100
ncpu =  4
period  = "IC86_2011-IC86_2021"

outdir = "/data/user/liruohan/powerlaw/trials/"

df_srcs = pd.read_pickle('/data/user/liruohan/powerlaw/single_sources.pkl')

idx1 = np.logical_and(df_srcs['dec'] > -5, df_srcs['ra'] > -360)
#idx2 = df_bass['TYPE_105'].str.contains('Sy')
#idx3 = df_bass['F2-10-intr'] > 27
df_srcs = df_srcs[idx1].sort_values(by='dec', ascending=False).copy(deep=True)


#--outdir $outdir --outfile $outfile --ncpus $ncpu --n_trials $ntrials --mean_ns $mean_ns --source_name $source_name
names = df_srcs.index.values 
mean_ns = [0]
tjob = 0
with open("trials.dag", "w") as fout:
    for tns in mean_ns:
        rss_seed = 0
        for job_id in range(njobs_per_ns):
            rss_seed += 1
            for src in names:
                tjob += 1 
    
                outfile = '{}_mean_ns{}_trials{}_rss{:}_.npy'.format(src, tns, ntrials, rss_seed)

                fout.write("JOB %s trials.submit\n" %(tjob))
                fout.write("VARS %s source=\"%s\" OUTDIR=\"%s\" OFILE=\"%s\"  CPU=\"%s\" NTRIALS=\"%s\" MEANNS=\"%s\" RSS_SEED=\"%s\"\n" %(tjob,src,outdir, outfile, str(ncpu), str(ntrials), str(tns), str(rss_seed)))