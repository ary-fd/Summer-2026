#!/bin/bash
#$ -cwd
#$ -j y
#$ -N abc_sim
#$ -l h_data=8G,h_rt=08:00:00
#$ -t 1-250
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.$TASK_ID.o$JOB_ID
set -x

cd /u/scratch/a/aryadini/Summer-2026/abc_test
source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env

python csfs_code/abc_simulate.py \
    --task-id  $SGE_TASK_ID \
    --n-sims   200 \
    --out-dir  abc_results