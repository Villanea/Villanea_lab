#!/bin/bash

# number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80

# Request to run on Villanea Lab node
#SBATCH --qos=blanca-villanea

# Request runtime:
#SBATCH --time=10:00:00

# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
#SBATCH --mem=150G

# Specify a job name:
#SBATCH --job-name test_mod

# Specify an output file
# %j is the job id
#SBATCH --output slurmOut/test_mod%j.out

#----- End of slurm commands ----

# Run a command
module purge
start=`date +%s` 
srun /bin/bash -c "source /curc/sw/anaconda3/latest; conda activate /projects/maea8398/.conda/envs/msprime_env; python3 network_pairwise.py 2000 15735363;" 
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
