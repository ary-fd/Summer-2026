#!/bin/bash
#$ -cwd
#$ -j y
#$ -N abc_test_cases
#$ -l h_data=20G,h_rt=02:00:00
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.o$JOB_ID

source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env

python csfs_code/generate_test_cases.py \
    --n-test 300 \
    --out-dir abc_results_10M \
    --seed 999