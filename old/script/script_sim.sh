#!/bin/bash
# This is an example batch script for slurm on Blanca
# 
# The commands for slurm start with #SBATCH
# All slurm commands need to come before the program 
# you want to run.  In this example, 'echo "Hello World!"
# is the command we are running.
#
# This is a bash script, so any line that starts with # is
# a comment.  If you need to comment out an #SBATCH line 
# use ##SBATCH 
#
# To submit this script to slurm do:
#    sbatch batch.script
#
# Once the job starts you will see a file MySerialJob-****.out
# The **** will be the slurm JobID

# --- Start of slurm commands -----------

#number of nodes

#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=80

# Request to run on Villanea Lab node
##SBATCH --account=pehu5714
#SBATCH --qos=blanca-villanea
##SBATCH --partition=shas
##sinteractive --qos=blanca-villanea --export=NONE --time=02:00:00

# Request an hour of runtime:
#SBATCH --time=04:00:00

# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
#SBATCH --mem=150G

# Specify a job name:
#SBATCH -J test-job

# Specify an output file
# %j is a special variable that is replaced by the JobID when 
# job starts
#SBATCH -o test_job-%j.out
#SBATCH -e test_job-%j.out

#----- End of slurm commands ----

# Run a command
module purge
start=`date +%s` 
# using conda env and run program in that env
srun /bin/bash -c "source /curc/sw/anaconda3/latest; conda activate /projects/pehu5714/popenv; python3 test_preseed.py;" 
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
