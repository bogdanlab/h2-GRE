#!/bin/sh

#$ -cwd
#$ -m a
#$ -l h_data=4G,h_rt=0:20:00
#$ -j y
#$ -o ./job_out

export PATH=~/anaconda2/bin:$PATH
export PYTHONNOUSERSITE=True

config=$1
param_i=$2

python ./sim_gwas.py beta --config=${config} --param-i=${param_i}
