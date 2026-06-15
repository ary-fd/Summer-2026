#!/bin/bash
#$ -cwd
#$ -j y
#$ -N prune_step2_chr20
#$ -l h_data=1G,h_rt=03:00:00
#$ -o $JOB_NAME.o$JOB_ID

plink \
  --bfile /u/project/sriram/data/1kg/phase3_plink/combined_final_no_dup \
  --chr 20 \
  --extract /u/scratch/a/aryadini/PCA/PCA_prune/prune.prune.in \
  --make-bed \
  --out /u/scratch/a/aryadini/PCA/PCA_prune/pruned_chr20