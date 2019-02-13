step=20
for chr_i in $(seq 1 22)
do
    out_dir=/u/project/pasaniuc/kangchen/projects/ukbb_local_h2g/out/ld/${chr_i}
    part_num=$(wc -l < ${out_dir}/part.info)
    qsub -t 1-${part_num}:${step} calc_ld.sh ${chr_i} ${step}
done
