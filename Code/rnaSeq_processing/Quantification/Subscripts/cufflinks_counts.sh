#!/bin/bash
#SBATCH -J Cufflinks-Counts                    # '##' means comment, 1 '#SBATCH' can be treated as the setting parameters
#SBATCH -p cpu  # can be compute0, compute2 or gpu, use "sinfo -Nel" to check the crosponding nodes under the spicified partion
###SBATCH -N 1
###SBATCH --ntasks-per-node=1
###SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH -t 100:00:00                # set time limit
###SBATCH -a 1-2%2  # you have 138 jobs to run, 15 jobs can be submitted at the same time.
input_fp=$1
gtf_fp=$2

PARAM1=$(sed -n -e "$SLURM_ARRAY_TASK_ID p" ${input_fp} | awk '{print $1}')
countdir_out="$(basename "$PARAM1" .bam)"

cufflinks \
-o $countdir_out \
-G $gtf_fp \
$PARAM1
