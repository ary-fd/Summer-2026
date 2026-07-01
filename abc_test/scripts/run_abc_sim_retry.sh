#!/bin/bash
#$ -cwd
#$ -j y
#$ -N abc_sim_retry
#$ -l h_data=32G,h_rt=08:00:00
#$ -t 1-375
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.$TASK_ID.o$JOB_ID

. /u/local/Modules/default/init/modules.sh

cd /u/scratch/a/aryadini/Summer-2026/abc_test
source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env

if [ ! -f abc_results_fix_10M/csfs_task${SGE_TASK_ID}.csv ]; then
    python csfs_code/abc_simulate.py \
        --task-id $SGE_TASK_ID \
        --n-sims 200 \
        --out-dir abc_results_fix_10M \
        --nsite 10000000
fi