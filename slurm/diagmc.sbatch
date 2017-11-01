#!/bin/bash

#SBATCH --job-name=AngulonDiagMC
#SBATCH --output=admc.log

# Number of CPU cores to use within one node
#SBATCH -c 1

# Define the number of hours the job should run. 
# Maximum runtime is limited to 10 days, ie. 240 hours
#SBATCH --time=48:00:00

# Define the amount of RAM used by your job in GigaBytes
# In shared memory applications this is shared among multiple CPUs
#SBATCH --mem=1G

# Send emails when a job starts, it is finished or it exits
#SBATCH --mail-user=giacomo.bighin@ist.ac.at
#SBATCH --mail-type=ALL

# Pick whether you prefer requeue or not. If you use the --requeue
# option, the requeued job script will start from the beginning, 
# potentially overwriting your previous progress, so be careful.
# For some people the --requeue option might be desired if their
# application will continue from the last state.
# Do not requeue the job in the case it fails.
#SBATCH --no-requeue

# Set the number of threads to the SLURM internal variable
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Load the respective software module you intend to use
#module load YourModuleHere

# Run the binary through SLURM's srun

if [ -z "$INIFILE" ] ;
then
        echo "Error: the \$INIFILE environment variable needs to be set!"
        echo "Usage: sbatch --export=INFILE=diagmc.ini diagmc.sbatch"
else
    srun --cpu_bind=verbose ../diagmc --tuning $INIFILE
fi
