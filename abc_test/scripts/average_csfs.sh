#!/bin/bash
#$ -cwd
#$ -j y
#$ -N run_csfs
#$ -l h_data=32G,h_rt=03:00:00
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.o$JOB_ID
set -x

cd /u/scratch/a/aryadini/Summer-2026/abc_test
source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env

python average_csfs.py