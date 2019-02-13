#$ -cwd
#$ -l h_data=16G,h_rt=3:00:00
#$ -j y
#$ -o ./job_out

export PATH=~/anaconda2/bin:$PATH
export PYTHONNOUSERSITE=True

for chr_i in $(seq 1 22)
do
    bfile=/u/project/pasaniuc/kangchen/projects/ukbb_local_h2g/out/genotype/${chr_i}
    out_dir=/u/project/pasaniuc/kangchen/projects/ukbb_local_h2g/out/ld/${chr_i}
    part_num=$(wc -l < ${out_dir}/part.info)
    if [ -e ${out_dir}/part.missing ]; then
        rm ${out_dir}/part.missing
    fi
    for part_i in $(seq 1 $part_num)
    do
        if [ ! -e ${out_dir}/part_${part_i}.npy ]; then
            echo ${part_i} >> ${out_dir}/part.missing
        fi
    done
    echo "chromosome: ${chr_i}"
    if [ -e ${out_dir}/part.missing ]; then
        cat ${out_dir}/part.missing | while read task_i
        do
            echo ${task_i} 
            python gre.py calc-ld --bfile=${bfile} --part-i=${task_i} --ld-dir=${out_dir}
        done
    else
        echo "complete"
    fi
done

