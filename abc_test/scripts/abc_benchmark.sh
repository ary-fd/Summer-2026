#!/bin/bash
#$ -cwd
#$ -j y
#$ -N abc_benchmark_tol03
#$ -l h_data=25G,h_rt=08:00:00
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.o$JOB_ID

. /u/local/Modules/default/init/modules.sh
module load intel/2020.4
module load R/4.2.2

Rscript csfs_code/run_abc_benchmark.R