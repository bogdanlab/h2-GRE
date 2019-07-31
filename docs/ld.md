# Computing chromosome wide inverse LD
h2-GRE requires the inverse LD matrix for each chromosome. For a small set of SNPs, computing the inverse LD matrix is straightforward; for example, given a genotype matrix in PLINK format,
```python
from pysnptools.snpreader import Bed
import numpy as np
from scipy import linalg

genotype = Bed('/path/to/bfile', count_A1=False).read().val
ld = np.corrcoef(genotype)
inv_ld = linalg.pinv(ld)
```
However, there are typically up to ~60K genotyped SNPs per chromosome and >100K individuals. In order to compute LD in memory, we divide the above into the several steps and parallelize over chunks of SNPs.

The following code computes the LD and inverse LD of a genotype matrix (in PLINK format). We assume that the genotype file only contains the information for one chromosome and the results will be at `ld_dir`.

We will use the following variables below:

- `ld_dir`: location where our results store
- `bfile`: plink genotype file for **only one** chromosome
- `mean_chunk_size`: when computing mean and standard deviation, how many SNPs we compute at a time, change this parameter according to memory.
- `ld_chunk_size`: when computing the LD, how many SNPs we compute at a time, change this parameter according to memory.  

## 1. Computing the mean and standard deviation of each SNP genotype
The mean and standard deviation of each genotype column is needed to standardize the genotype matrix.

```bash
# when number of individuals = ~350K, setting mean_chunk_size=500 takes up about 8G memory
python gre.py mean-std --bfile=${bfile} --chunk-size=${mean_chunk_size} --ld-dir=${ld_dir}
```

## 2. Dividing the LD computations into multiple tasks
Next, we compute the LD matrix which has dimensions number-of-SNPs x number-of-SNPs. 
For example, we can cut the large task of computing a 60k x 60k SNP covariance matrix into 3600 small tasks, each of which computes a 1k x 1k SNP covariance matrix.

```bash
# this step will create a file `part.info` in ${ld_dir}, which contains the information for each sub-task
python gre.py cut-ld --bfile=${bfile} --chunk-size=${ld_chunk_size} --ld-dir=${ld_dir}
python gre.py check-ld ${ld_dir}
```


Now we compute covariance matrix for each part. You can compute the covariance matrix sequentially as follows:

```bash
part_num=$(wc -l < ${ld_dir}/part.info)
for part_i in $(seq 1 ${part_num})
do
    python gre.py calc-ld --bfile=${bfile} --part-i=${part_i} --ld-dir=${ld_dir}
done
```


## 3. Computing chromosome-wide inverse LD
We merge the covariance matrices from the previous step to obtain the chromosome-wide LD matrices. First, check if each task is complete. If a task is not complete, check `${ld_dir}/part.missing` and re-submit these.

For a 60k x 60k covariance matrix, computing the inverse via singular value decomposition (SVD) takes about 60G memory and 4 hours. For a 10k x 10k covariance matrix, this step requires 6G memory and 1 hour. Note that the memory and time complexity are independent of the number of individuals.
```bash
part_num=$(wc -l < ${ld_dir}/part.info)
if [ -e ${ld_dir}/part.missing ]; then
    rm ${ld_dir}/part.missing
fi
for part_i in $(seq 1 ${part_num})
do
    if [ ! -e ${ld_dir}/part_${part_i}.npy ]; then
        echo ${part_i} >> ${ld_dir}/part.missing
    fi
done

if [ -e ${ld_dir}/part.missing ]
then
    echo "the following parts are missing, please rerun calc_ld.sh"
    cat ${ld_dir}/part.missing
else
    python gre.py merge-ld --bfile=${bfile} --ld-dir=${ld_dir}
fi
```
After this step, you will see the inverse of LD matrix and the rank in folder `${ld_dir}`.