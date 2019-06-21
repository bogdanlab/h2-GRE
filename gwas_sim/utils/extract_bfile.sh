# extract the first 10k individuals
. /u/local/Modules/default/init/modules.sh
module load plink
bfile=/u/project/pasaniuc/pasaniucdata/DATA/UKBiobank/array/allchr.unrelatedbritishqced.mafhwe
awk '{print $1}' ${bfile}.fam | head -n 10000 > ./out/indv_10k.fam
plink --keep-fam ./out/indv_10k.fam --make-bed --bfile ${bfile} --out ./out/10k

for chr_i in $(seq 22)
do
plink --make-bed --bfile ./out/10k --keep-allele-order --out ./out/10k.chr${chr_i} --chr ${chr_i}
done
