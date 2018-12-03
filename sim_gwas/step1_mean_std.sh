#!/bin/sh

#$ -cwd
#$ -m a
#$ -l h_data=8G,h_rt=0:30:00,highp
#$ -j y
#$ -o ./job_out

export PATH=~/anaconda2/bin:$PATH
export PYTHONNOUSERSITE=True

config=$1
chr_i=$SGE_TASK_ID

python ./sim_gwas.py mean-std --config=${config} --chr-i=${chr_i}
