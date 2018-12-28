#!/bin/bash

#$ -cwd
#$ -m a
#$ -l h_data=4G,h_rt=0:20:00
#$ -j y
#$ -o ./job_out

export PYTHONNOUSERSITE=True
export PATH=~/anaconda2/bin:$PATH

config=$1
param_i=$2

python ./sim_gwas.py phe --config=${config} --param-i=${param_i}
