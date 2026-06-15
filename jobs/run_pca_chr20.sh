#!/bin/bash
#$ -cwd
#$ -j y
#$ -N run_pca_chr20
#$ -l h_data=1G,h_rt=03:00:00
#$ -o $JOB_NAME.o$JOB_ID

plink \
  --bfile /u/scratch/a/aryadini/PCA/PCA_prune/pruned_chr20 \
  --pca 10 \
  --out /u/scratch/a/aryadini/pca_out_chr20

