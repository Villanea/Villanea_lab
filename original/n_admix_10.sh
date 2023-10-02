#!/bin/bash

# number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80

# Request to run on Villanea Lab node
#SBATCH --qos=blanca-villanea

# Request runtime:
#SBATCH --time=01-00:00:00

# Specify a job name:
#SBATCH --job-name n_admix_10

# Specify an output file
#SBATCH --output n_admix_10_%j.out

# Run a command
module purge
start=`date +%s` 
# using conda env and run program in that env
srun /bin/bash -c "source /curc/sw/anaconda3/latest; conda activate /projects/maea8398/.conda/envs/msprime_original; python n_admix_10.py;" 
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.

