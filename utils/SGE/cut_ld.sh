export PATH=~/anaconda2/bin:$PATH
export PYTHONNOUSERSITE=True

chunk_size=1000

for chr_i in $(seq 1 22)
do
    bfile=/u/project/pasaniuc/kangchen/projects/ukbb_local_h2g/out/genotype/${chr_i}
    out_dir=/u/project/pasaniuc/kangchen/projects/ukbb_local_h2g/out/ld/${chr_i}
    python gre.py cut-ld --bfile=${bfile} --chunk-size=${chunk_size} --ld-dir=${out_dir}
done
