#!/bin/bash
set -ex
targetlist=("/u/scratch/h/haroldzw/CSFS/NPE/CSFS/prufer_with_ghost.csfs")
simTAGlist=(
  "/u/scratch/h/haroldzw/CSFS/Pseudo_truth/narrowprior.csv100000"
  "/u/scratch/h/haroldzw/CSFS/Pseudo_truth/wideprior.csv100000"
)
# noramlizelist=("sum" "sum_pre_cut")
# distancelist=("wasserstein")

# noramlizelist=("none" "median" "median_pre_cut" "sum" "sum_pre_cut")
# distancelist=("euclidean")

noramlizelist=("sum")
distancelist=("euclidean")
methodlist=('ridge')
# 'rejection', 'loclinear', 'neuralnet', 'ridge'
# Testing
# targetlist=("/u/scratch/h/haroldzw/CSFS/NPE/CSFS/prufer_with_ghost.csfs")
# simTAGlist=(
#   "/u/scratch/h/haroldzw/CSFS/Pseudo_truth/narrowprior.csv100000"
# )
# noramlizelist=("sum_pre_cut")
# distancelist=("wasserstein")
for target in "${targetlist[@]}"; do
  for simTAG in "${simTAGlist[@]}"; do
    for normalization_method in "${noramlizelist[@]}"; do
      for distance_method in "${distancelist[@]}"; do
        for abc_method in "${methodlist[@]}"; do
        qsub -cwd \
          -e STDERR/ \
          -o STDOUT/ \
          -V \
          -N "run_${normalization_method}_${abc_method}" \
          -l h_data=64G,time=24:00:00 \
          -m e \
          -M haroldzw \
          -b y "bash qc_individual.sh $target $simTAG $normalization_method $distance_method $abc_method"
        done
      done
    done
  done
done