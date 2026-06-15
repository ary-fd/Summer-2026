#!/bin/bash
set -ex
#$ -V
#$ -cwd
#$ -j y
#$ -N pca_chr15
#$ -l h_data=16G,h_rt=06:00:00
#$ -o /u/scratch/a/aryadini/logs/$JOB_NAME.o$JOB_ID

DATA=/u/project/sriram/data/1kg/phase3_plink/combined_final_no_dup
SCRATCH=/u/scratch/a/aryadini/PCA

# Step 1 — LD pruning
plink \
  --bfile $DATA \
  --chr 15 \
  --indep-pairwise 50 5 0.2 \
  --out $SCRATCH/prune_chr15

# Step 2 — extract pruned variants
plink \
  --bfile $DATA \
  --chr 15 \
  --extract $SCRATCH/prune_chr15.prune.in \
  --make-bed \
  --out $SCRATCH/pruned_chr15

# Step 3 — PCA
plink \
  --bfile $SCRATCH/pruned_chr15 \
  --pca 10 \
  --out $SCRATCH/PCA_out/pca_out_chr15

# Step 4 — plot
python3 /u/scratch/h/haroldzw/Student_project/PCA/plot_pca.py \
  --eigenvec $SCRATCH/PCA_out/pca_out_chr15.eigenvec \
  --eigenval $SCRATCH/PCA_out/pca_out_chr15.eigenval \
  --panel /u/project/sriram/data/1kg/phase3_plink/integrated_call_samples_v3.20130502.ALL.panel \
  --fam $DATA.fam \
  --out $SCRATCH/PCA_plots/pca_chr15.png