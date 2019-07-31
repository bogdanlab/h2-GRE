# Perform Association Study, i.e. OLS regression

After getting the data, we are now ready to perform the OLS regression.

To use GRE, be sure to perform the association studies with **the exactly same genotype** you use to compute the LD. We leave GRE with reference panel LD for future work.

```bash
num_cols=$(head -1 $covar_file | awk '{print NF}')
num_cols=$((num_cols-2))

for chr_i in $(seq 1 22)
do
plink \
  --allow-no-sex \
  --bfile ${all_bfile} \
  --chr ${chr_i} \
  --ci 0.95 \
  --covar ${covar_file} \
  --covar-number 1-${num_cols} \
  --linear hide-covar \
  --out ./assoc/chr${chr_i} \ 
  --pheno ${pheno_file}
done
```


You will get 22 association files after this step. Each contains the OLS sumstats for SNPs for 22 chromosomes.