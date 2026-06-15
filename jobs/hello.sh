#!/bin/bash
#$ -cwd
#$ -j y
#$ -N hello_hostname
#$ -l h_data=1G,h_rt=00:05:00
#$ -o $JOB_NAME.o$JOB_ID

echo "Hostname:    $(hostname)"
echo "Date:        $(date)"
echo "Working dir: $(pwd)"
echo "User:        $(whoami)"
echo "Job ID:      $JOB_ID"