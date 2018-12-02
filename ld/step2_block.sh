#$ -cwd
#$ -m a
#$ -l h_data=16G,h_rt=0:06:00
#$ -j y
#$ -o ./job_out

export PYTHONNOUSERSITE=True
export PATH=~/anaconda2/bin:$PATH

chr_i=$1

bfile=$(cat ./metadata.json | python -c "import sys, json; print json.load(sys.stdin)['genotype'].format(${chr_i})")
chunk_size=$(cat ./metadata.json | python -c "import sys, json; print json.load(sys.stdin)['ld_chunk_size']")

snp_num=$(wc -l < ${bfile}.bim)

chunk_num=$(( ($snp_num + $chunk_size - 1) / $chunk_size))

task_i=$(($SGE_TASK_ID-1))

row_i=$((${task_i}/${chunk_num}))
col_i=$((${task_i}%${chunk_num}))

mean_std_file=./out/mean_std/chr${chr_i}.txt

python get_ld.py part --bfile=$bfile --mean-std-file=$mean_std_file --row-i $row_i --col-i $col_i --chunk-size $chunk_size --out-dir ./out/part/chr${chr_i}

