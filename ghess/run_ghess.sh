#$ -cwd
#$ -m a
#$ -l h_data=24G,h_rt=3:00:00,highp
#$ -j y
#$ -o ./job_out


export PYTHONNOUSERSITE=True
export PATH=~/anaconda2/bin:$PATH

inv_ld_dir=../ld/out
gwas_dir=../sim_gwas
mkdir ./out
python ghess.py estimate --inv-ld-dir=${inv_ld_dir} --gwas-dir=${gwas_dir} --out-dir=./out
