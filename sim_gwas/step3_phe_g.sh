#!/bin/sh

#$ -cwd
#$ -m a
#$ -l h_data=12G,h_rt=4:00:00
#$ -j y
#$ -o ./job_out

export PATH=~/anaconda2/bin:$PATH
export PYTHONNOUSERSITE=True

config=$1
param_i=$2
chr_i=$SGE_TASK_ID

python ./sim_gwas.py phe-g --config=${config} --param-i=${param_i} --chr-i=${chr_i}
