#!/bin/bash
#$ -cwd
#$ -j y
#$ -N train_nn
#$ -l h_data=16G,h_rt=02:00:00
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.o$JOB_ID
set -x

cd /u/scratch/a/aryadini/Summer-2026/abc_test
source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env

python csfs_code/train_nn.py \
    --results-dir abc_results_10M \
    --observed    csfs/averaged_csfs.npy \
    --out-dir     abc_results_10M