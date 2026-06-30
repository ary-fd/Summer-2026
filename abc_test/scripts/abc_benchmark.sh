#!/bin/bash
#$ -cwd
#$ -j y
#$ -N abc_benchmark
#$ -l h_data=20G,h_rt=04:00:00
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.o$JOB_ID

. /u/local/Modules/default/init/modules.sh

source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env

module load R/4.2.2
Rscript csfs_code/run_abc_benchmark.R