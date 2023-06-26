import ndmix9_newmsp as ndmix9
import concurrent.futures
import os

if __name__ == '__main__':
    try:
        os.makedirs("data")
    except FileExistsError:
        pass
    try:
        os.makedirs("data/job"+os.environ["SLURM_JOB_ID"])
    except FileExistsError:
        pass
    try:
        os.makedirs("data/job"+os.environ["SLURM_JOB_ID"]+"/mat")
    except FileExistsError:
        pass
    try:
        os.makedirs("data/job"+os.environ["SLURM_JOB_ID"]+"/stat")
    except FileExistsError:
        pass
    try:
        os.makedirs("data/job"+os.environ["SLURM_JOB_ID"]+"/outfile")
    except FileExistsError:
        pass
    for model_num in range(1, 6):
        num_replicates = 40
        # sim_num, model_num
        sim_num = [[i,model_num] for i in range(1,num_replicates+1)]
        with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            results = executor.map(ndmix9.worker, sim_num)

