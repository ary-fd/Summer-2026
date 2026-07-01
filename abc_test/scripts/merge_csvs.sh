bash#!/bin/bash
#$ -cwd
#$ -j y
#$ -N merge_csvs
#$ -l h_data=16G,h_rt=01:00:00
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.o$JOB_ID

source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env

python csfs_code/merge_csvs.py \
    --out-dir abc_results_fix_10M \
    --n-tasks 375