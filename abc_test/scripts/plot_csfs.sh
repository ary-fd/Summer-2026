#!/bin/bash
#$ -cwd
#$ -j y
#$ -N plot_csfs
#$ -l h_data=10G,h_rt=03:00:00
#$ -o /u/scratch/a/aryadini/Summer-2026/abc_test/logs/$JOB_NAME.o$JOB_ID
set -x

cd /u/scratch/a/aryadini/Summer-2026/abc_test
source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env


python plot_csfs.py \
    --csfs  csfs/averaged_csfs.npy \
    --out   csfs/csfs_plot.pdf \
    --title "CSFS: Africa conditioned on Neanderthal" \
    --label "Africa | Neanderthal derived" 2>&1