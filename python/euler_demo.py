# encoding: utf-8
from __future__ import print_function
import os
import sys
import pickle
import cobra

def main(number_chunks):

	job_index = int(os.environ.get("LSB_JOBINDEX"))
	number_jobs = int(os.environ.get("LSB_JOBINDEX_END"))

	for chunk in range(job_index, number_chunks + 1, number_jobs):

	    result = "job #{} processes chunk {} out of {}".format(job_index, chunk, number_chunks)
	    
	    RESULT_DIR = "/cluster/scratch/michielk"
	    with open(RESULT_DIR+'/part'+str(job_index)+'.p', 'wb') as handle:
	        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL) 

if __name__ == '__main__':
    number_chunks = int(sys.argv[1]) 
    main(number_chunks)


### 100 jobs over 8 works:
### bsub -J "demo_array[1-8]" -R light python demo.py 100