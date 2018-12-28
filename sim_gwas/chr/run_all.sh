export PYTHONNOUSERSITE=True
export PATH=~/anaconda2/bin:$PATH

config=./config.json

mkdir ./out
mkdir ./out/mean_std
mkdir ./out/beta
mkdir ./out/phe_g
mkdir ./out/phe
mkdir ./out/zsc

mean_std_num=$(find ./out/mean_std/ -type f | wc -l)
if [ "$mean_std_num" -ne 1 ]; then
    qsub step1_mean_std.sh ${config}
    echo "computing the mean and stddev right now, please wait for these jobs to finish"
    exit 1
fi

params_num=$(cat ./config.json | python -c "import sys, json; print len(json.load(sys.stdin)['params'])")

for param_i in $(seq 0 $(( ${params_num} - 1 )))
do
    qsub -N sim_gwas_step2_beta_param_${param_i} step2_beta.sh ${config} ${param_i}
    qsub -hold_jid sim_gwas_step2_beta_param_${param_i} -N sim_gwas_step3_phe_g_param_${param_i} step3_phe_g.sh ${config} ${param_i}
    qsub -hold_jid sim_gwas_step3_phe_g_param_${param_i} -N sim_gwas_step4_phe_param_${param_i} step4_phe.sh ${config} ${param_i}
    qsub -hold_jid sim_gwas_step4_phe_param_${param_i} -N sim_gwas_step5_zsc_param_${param_i} step5_zsc.sh ${config} ${param_i}
done
