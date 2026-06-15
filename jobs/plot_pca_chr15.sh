#!/bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -N plot_pca_chr15
#$ -l h_data=10G,h_rt=03:00:00
#$ -o logs/$JOB_NAME.o$JOB_ID

python3 /u/scratch/h/haroldzw/Student_project/PCA/plot_pca.py \
  --eigenvec /u/scratch/a/aryadini/PCA/PCA_out/pca_out_chr15.eigenvec \
  --eigenval /u/scratch/a/aryadini/PCA/PCA_out/pca_out_chr15.eigenval \
  --panel /u/project/sriram/data/1kg/phase3_plink/integrated_call_samples_v3.20130502.ALL.panel \
  --fam /u/project/sriram/data/1kg/phase3_plink/combined_final_no_dup.fam \
  --out /u/scratch/a/aryadini/PCA/PCA_plots/pca_chr15.png