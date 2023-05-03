import ndmix9_oldrand as ndmix9
import concurrent.futures
import os

if __name__ == '__main__':

    num_replicates = 2000
    model_num = 2
    #get pre-generated seed
    with open('random_numbers.txt', 'r') as f:
        num_list = [int(line.strip()) for line in f.readlines()]
    # sim_num, model_num, seed
    sim_num = [[i,model_num,num_list[(i-1)+2000*(model_num-1)]] for i in range(1,num_replicates+1)]

    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
            results = executor.map(ndmix9.worker, sim_num)