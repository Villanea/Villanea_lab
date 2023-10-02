import ndmix9_updated as ndmix9
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
    for model_num in [1, 2, 3, 4, 5]:
        num_replicates = 80
        # sim_num, model_num
        sim_num = [[i,model_num] for i in range(1,num_replicates+1)]
        with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            results = executor.map(ndmix9.worker, sim_num)

