#!/bin/sh

#$ -cwd
#$ -m a
#$ -l h_data=8G,h_rt=0:30:00,highp
#$ -j y
#$ -o ./job_out

export PATH=~/anaconda2/bin:$PATH
export PYTHONNOUSERSITE=True

chr_i=$SGE_TASK_ID

bfile=$(cat ./metadata.json | python -c "import sys, json; print json.load(sys.stdin)['genotype'].format(${chr_i})")
chunk_size=$(cat ./metadata.json | python -c "import sys, json; print json.load(sys.stdin)['mean_std_chunk_size']")

python ./get_ld.py mean-std --bfile=${bfile} --chunk-size=${chunk_size} --out ./out/mean_std/chr"$chr_i".txt

