#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 10:02:37 2018

@author: Zarathustra
"""
import os, pickle 
SERVER_SCRATCH = '/cluster/scratch/michielk/'

def main(n_parameterizations):
    
    topologies = list(range(34987))

    job_index = int(os.environ.get("LSB_JOBINDEX"))
    number_jobs = int(os.environ.get("LSB_JOBINDEX_END"))
    chunk = topologies[job_index-1::number_jobs]
        
    ### iterate over the topologies in this chunk
    results = dict()
    for j in range(0,len(chunk)): # 35000 / 200 = 175, 25 saved
        
      
        result = list()
        for _ in range(n_parameterizations): # multiple parameterizations per model topology
            
            result += [j]
            
        results[j] = result

    with open(SERVER_SCRATCH+'/part'+str(job_index)+'.p', 'wb') as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL) 

if __name__ == '__main__':
    n_parameterizations = 1#int(sys.argv[2])
    main(n_parameterizations)