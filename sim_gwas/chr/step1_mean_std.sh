#!/bin/sh

#$ -cwd
#$ -m a
#$ -l h_data=8G,h_rt=0:30:00,highp
#$ -j y
#$ -o ./job_out

export PATH=~/anaconda2/bin:$PATH
export PYTHONNOUSERSITE=True

config=$1

python ./sim_gwas.py mean-std --config=${config}
