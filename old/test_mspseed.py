import ndmix9_oldrand as ndmix9
import concurrent.futures
import os

if __name__ == '__main__':

    num_replicates = 2000
    model_num = 1
    # sim_num, model_num
    sim_num = [[i,model_num] for i in range(1,num_replicates+1)]

    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            results = executor.map(ndmix9.worker, sim_num)