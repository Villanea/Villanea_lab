import os
import sys
import glob
import numpy as np

jobID = sys.argv[1]

try: os.makedirs("data/job"+jobID+"/corrupt")
except FileExistsError: pass
try: os.makedirs("data/job"+jobID+"/corrupt/mat")
except FileExistsError: pass

def check():
    for model_num in range(1, 6):
        count = 0
        for i, file in enumerate(glob.glob('data/job%s/mat/symmetry_matrix_%d*' % (jobID, model_num))):
            corrupt = False
            try: 
                mat = np.loadtxt(file)
                if mat.shape != (64, 64): corrupt = True
            except ValueError: 
                corrupt = True
            if corrupt:
                os.system("mv %s data/job%s/corrupt/mat" %(file, jobID))
            else:
                count += 1
        print("model %d: %d" % (model_num, count))
        
check()
