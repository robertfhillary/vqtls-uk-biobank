#!/bin/bash -l
 
#SBATCH --job-name=chr22
#SBATCH --cpus-per-task=1
#SBATCH --time 2-00:00:00
#SBATCH --qos=advanced_users_limits
#SBATCH --mail-user=xxxxxx
#SBATCH --mail-type=END,FAIL

#SBATCH --array=1-184
 
# Getting $SLURM_ARRAY_TASK_ID line from parameters file
PER_TASK=8

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM

# Run the loop of runs for this task.
for (( run=$START_NUM; run<=END_NUM; run++ )); do
  echo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run
  PARAM=$(awk "NR==$run" /home/rhillary/disfiles2.txt)
  srun bash osca_batch22.sh $PARAM
done

