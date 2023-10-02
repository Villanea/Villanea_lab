#!/bin/bash

# number of nodes
#SBATCH --nodes=25
#SBATCH --ntasks=25
#SBATCH --cpus-per-task=40

# Request to run on Villanea Lab node
##SBATCH --qos=blanca-villanea
#SBATCH --qos=preemptable

# Request runtime:
#SBATCH --time=01-00:00:00

# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
##SBATCH --mem=150G

# Specify a job name:
#SBATCH --job-name test_sim

# Specify an output file
# %j is the job id
#SBATCH --output slurmOut/test_sim%j.out

# Run a command
module purge
start=`date +%s` 
# using conda env and run program in that env
srun /bin/bash -c "source /curc/sw/anaconda3/latest; conda activate /projects/maea8398/.conda/envs/msprime_env; python3 test_sim.py;" 
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.

