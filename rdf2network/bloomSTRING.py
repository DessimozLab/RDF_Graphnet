#!/usr/bin/env python
# coding: utf-8

# In[1]:


from dask.distributed import fire_and_forget
from dask.distributed import Client, Variable , Queue , Lock ,LocalCluster
from dask_jobqueue import SLURMCluster
from dask.distributed import  utils_perf
from dask.distributed import Client, LocalCluster
import dask
import redis
from bloom_filter2 import BloomFilter
import lzma
from dask import dataframe as dd


# In[2]:


#!/usr/bin/env python
# coding: utf-8
import sys
sys.path.append('../../..')
#sys.path.append( '/home/cactuskid13/miniconda3/pkgs/')
print(sys.path)
import ete3
import random
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import itertools
import dask
import warnings
from dask import dataframe as dd
import pickle
path = '/work/FAC/FBM/DBC/cdessim2/default/dmoi/'


# In[3]:


from dask.distributed import fire_and_forget
from dask.distributed import Client, Variable , Queue , Lock ,LocalCluster
from dask_jobqueue import SLURMCluster
from dask.distributed import  utils_perf
from dask.distributed import Client, LocalCluster
import dask
import redis
from bloom_filter2 import BloomFilter
import lzma
from dask import dataframe as dd

#using distributed computation on a slurm cluster here. This is my particular config. You will need to alter this: https://distributed.dask.org/en/stable/
NCORE = 4
print('deploying cluster')
cluster = SLURMCluster(
    #change theses settings for your cluster
    walltime='4:00:00',
    n_workers = NCORE,
    cores=NCORE,
    processes = NCORE,
    interface='ib0',
    memory="20GB",
    #job_script_prologue=[
    #' source /work/FAC/FBM/DBC/cdessim2/default/dmoi/miniconda3/etc/profile.d/conda.sh  ' ,
    #'conda activate ML2'
    #],
    #scheduler_options={'interface': 'ens2f0' },
    #if gpu node
    scheduler_options={'interface': 'ens3f0np0' },
    #extra=[ "--lifetime-stagger", "4m"]
)
print(cluster.job_script())
print(cluster)
cluster.scale(jobs = 100)
print(cluster.dashboard_link)
client = Client(cluster , timeout='450s' , set_as_default=True )


# In[10]:


#find which species each of the cogs has an interaction in
link_df = dd.read_csv(path + 'datasets/STRING/protein.links.full.v11.5.txt',  blocksize=15e6 , header = 0, sep = ' ')
print(link_df)


# In[ ]:


#compute bloom filters for protein pairs

import tqdm 

@dask.delayed
def protlinks_species( df ):
    df = df[~df['protein1'].isna()]
    df = df[~df['protein2'].isna()]
    df.protein1 = df.protein1.map(lambda x:str(x))
    df.protein2 = df.protein2.map(lambda x:str(x))
    df['protlinks'] = df.protein1 + '_' + df.protein2 
    ret = set(df.protlinks.unique())
    return ret

@dask.delayed
def return_filter(protlinks, verbose = True):
    if type( protlinks ) == tuple:
        protlinks = protlinks[0]
    b=BloomFilter(max_elements=10**8, error_rate=0.001 ,start_fresh = True)
    for p in protlinks:
        b.add( p )
    retlen = len(protlinks)
    return   b , retlen

@dask.delayed
def sumfilter(f1,f2, total ):
    if type( f1 ) == tuple:
        f1 = f1[0]
    if type( f2 ) == tuple:
        f2 = f2[0]
    f3 = f1.__ior__(f2)
    return f3 , total

def treesum(totalfilter):
    #print(len(totalfilter))
    while len(totalfilter)>1:
        next_round= []
        for i in range(0,len(totalfilter),2):
            if i+1 < len(totalfilter):
                next_round.append( sumfilter( totalfilter[i][0] , totalfilter[i+1][0] , totalfilter[i][1]+totalfilter[i+1][1]  ) )
        if len(totalfilter) % 2 !=0:
            next_round.append(totalfilter[-1])
        totalfilter = next_round
        #print(len(totalfilter))
    return totalfilter

b=BloomFilter(max_elements=10**8, error_rate=0.001 ,start_fresh = True)
partitions  = link_df.to_delayed()
print('map cogs')
res1 = [ protlinks_species(p) for p in partitions ]
print('done')
print('make filters')
res2 = [ return_filter(p) for p in res1 ]


finals =[]

with tqdm.tqdm(total=int(len(res2)/512)+1) as pbar:
    for chunk in range(int(len(res2)/512)+1):
        #print(chunk*1024)
        res3 = res2[chunk*512:(chunk+1)*512]
        res4 = treesum(res3)
        res4 = dask.compute(res4)
        finals.append(res4[0])
        pbar.set_description('processed: %d' % (1 + chunk ))
        pbar.update(1)  

#dask.compute(*finals)

with open('bloomfinal_big.pkl' , 'wb' ) as finalout:
    finalout.write(pickle.dumps(finals))

