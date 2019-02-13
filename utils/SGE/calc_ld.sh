#$ -cwd
#$ -l h_data=16G,h_rt=0:30:00
#$ -j y
#$ -o ./job_out

export PATH=~/anaconda2/bin:$PATH
export PYTHONNOUSERSITE=True

chr_i=$1
chunk_size=$2

out_dir=/u/project/pasaniuc/kangchen/projects/ukbb_local_h2g/out/ld/${chr_i}
bfile=/u/project/pasaniuc/kangchen/projects/ukbb_local_h2g/out/genotype/${chr_i}

for i in $(seq 1 ${chunk_size})
do
    task_id=$((SGE_TASK_ID + i - 1))
    echo "chr_i: ${chr_i}, task_id: ${task_id}"
    python gre.py calc-ld --bfile=${bfile} --part-i=${task_id} --ld-dir=${out_dir}
done
