#!/bin/bash

# number of nodes
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1

# Request to run on Villanea Lab node
#SBATCH --qos=blanca-villanea

# Request runtime:
#SBATCH --time=00:01:00

# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
#SBATCH --mem=150G

# Specify a job name:
#SBATCH --job-name template

# Specify an output file
# %j is the job id
#SBATCH --output slurmOut/template%j.out

# Run a command
module purge
start=`date +%s` 
# using conda env and run program in that env
srun /bin/bash -c "source /curc/sw/anaconda3/latest; conda activate /home/maea8398/.conda/envs/msprime-env; python3 template.py;" 
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
