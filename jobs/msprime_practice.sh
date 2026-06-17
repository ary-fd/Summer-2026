#!/bin/bash
#$ -cwd
#$ -j y
#$ -N practice_msprime
#$ -l h_data=1G,h_rt=03:00:00
#$ -o $JOB_NAME.o$JOB_ID

source /u/project/sriram/aryadini/miniforge3/etc/profile.d/conda.sh
conda activate msprime_env

#This script runs the `archiesim.py` simulation tool with various input parameters.
# 
# Parameters:
# - `-r $rep`: Specifies the number of replicates to simulate.
# - `-t $model`: Defines the population history model to use.
# - `-d ${model}.yaml`: Path to the YAML configuration file for the specified model.
# - `-s ${json}.json`: Path to the JSON file containing simulation settings.
# - `-p $pop`: Specifies the population where we're trying to infer introgression ancestry.
# - `-i $introgress`: Indicates the introgression source population, aka Neanderthals
# - `-f $reference`: Indicates the outgroup or reference population with no introgression.
# - `-ns $nsite`: Number of sites to simulate.
# - `--mut_rate $mu`: Mutation rate for the simulation.
# - `--rec_rate $r`: Recombination rate for the simulation.
# - `--target $target`: Specifies the population into which introgression occurs. Can be the same or not as $pop.
#
# Ensure all variables (e.g., $rep, $model, $json, $pop, etc.) are defined in the environment
# or passed appropriately before running this script.

pop=Eurasia
introgress=Neanderthal
reference=Africa
target=Eurasia
model=practice
json=sample
nsite=100000
rep=1
dir=/u/project/sriram/aryadini/msprime_practice

archiesim.py \
    -r $rep \
    -t $model \
    -d /${directory}/${model}.yaml \
    -s /${directory}/${json}.json \
    -p $pop \
    -i $introgress \
    -f $reference \
    -ns $nsite \
    --target $target

    #--mut_rate $mu \
    #--rec_rate $r \