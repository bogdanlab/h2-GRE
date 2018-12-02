#$ -cwd
#$ -m a
#$ -j y
#$ -o ./job_out

export PYTHONNOUSERSITE=True
export PATH=~/anaconda2/bin:$PATH

chr_i=$1

python get_ld.py multiply --in-dir=./out/part/chr${chr_i} --out-path=./out/chr${chr_i}

