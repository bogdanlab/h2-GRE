mkdir out/part

chunk_size=$(cat ./metadata.json | python -c "import sys, json; print json.load(sys.stdin)['ld_chunk_size']")

for chr_i in $(seq 1 22)
do
    mkdir out/part/chr${chr_i}
    bfile=$(cat ./metadata.json | python -c "import sys, json; print json.load(sys.stdin)['genotype'].format(${chr_i})")
    snp_num=$(wc -l < ${bfile}.bim)
    chunk_num=$(( ($snp_num + $chunk_size - 1) / $chunk_size))
    
    task_num=$(( ${chunk_num} * ${chunk_num} ))
    qsub -t 1-$task_num step2_block.sh $chr_i
done

