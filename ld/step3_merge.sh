#$ -cwd
#$ -m a
#$ -j y
#$ -o ./job_out

export PYTHONNOUSERSITE=True
export PATH=~/anaconda2/bin:$PATH

chr_i=$1

bfile=$(cat ./metadata.json | python -c "import sys, json; print json.load(sys.stdin)['genotype'].format(${chr_i})")
snp_num=$(wc -l < ${bfile}.bim)
chunk_size=$(cat ./metadata.json | python -c "import sys, json; print json.load(sys.stdin)['ld_chunk_size']")

out_dir=./out/part/chr${chr_i}

python get_ld.py merge --snp-num=${snp_num} --chunk-size=${chunk_size} --out-dir=${out_dir}
