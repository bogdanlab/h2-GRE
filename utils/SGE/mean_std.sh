#!/bin/sh

#$ -cwd
#$ -m a
#$ -l h_data=8G,h_rt=1:00:00,highp
#$ -j y
#$ -o ./job_out

export PATH=~/anaconda2/bin:$PATH
export PYTHONNOUSERSITE=True


chr_i=$SGE_TASK_ID
genotype_dir=/u/project/pasaniuc/kangchen/projects/ukbb_local_h2g/out/genotype
out_dir=/u/project/pasaniuc/kangchen/projects/ukbb_local_h2g/out/ld/${chr_i}
mkdir -p ${out_dir}
python gre.py mean-std --bfile=${genotype_dir}/${chr_i} --chunk-size=250 --ld-dir=${out_dir}
