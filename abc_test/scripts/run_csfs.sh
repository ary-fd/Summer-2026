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

for rep in 0 1 2 3 4; do
    python compute_csfs.py \
        --african-geno eigenstrat/ghost_sim_rep${rep}Africa.geno \
        --african-snp  eigenstrat/ghost_sim_rep${rep}Africa.snp \
        --african-ind  eigenstrat/ghost_sim_rep${rep}Africa.ind \
        --neanderthal-geno eigenstrat/ghost_sim_rep${rep}Neanderthal.geno \
        --neanderthal-snp  eigenstrat/ghost_sim_rep${rep}Neanderthal.snp \
        --neanderthal-ind  eigenstrat/ghost_sim_rep${rep}Neanderthal.ind \
        --save-npy csfs/rep${rep}
done
