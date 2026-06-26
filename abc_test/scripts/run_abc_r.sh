#!/bin/bash
#$ -cwd
#$ -j y
#$ -N run_abc_r
#$ -l h_data=32G,h_rt=04:00:00
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.o$JOB_ID
set -x

cd /u/scratch/a/aryadini/Summer-2026/abc_test
source /u/local/Modules/default/init/bash
module load R/4.2.2

Rscript csfs_code/run_abc.R