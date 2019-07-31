# h2-GRE

h2-GRE is a software package implementing method for accurate estimation of SNP-heritability from biobank-scale data.

Estimating h2-GRE consists of the following steps, please click link for a tutorial for each part. In these tutorial, we are aiming to introduce each step. So for simplicity, the scripts are in sequential version, i.e., these scripts do not utilize the parallel computing on the cluster.

1. [Prepare the data](input_data.md).
2. [Perform GWAS (i.e. OLS regression) using individual level data.](ols.md)
3. [Compute chromosome-wide pseudoinverse of LD](ld.md).
4. [Calculate h2-GRE](estimate_h2gre.md).

We also prepare a [step-by-step tutorial](biobank_demo.md) which demonstrate how we performed analyses on UK Biobank data (#individuals=337k, genome-wide #SNPs=600k, largest chromosome #SNPs=60k).

Details on each step can be also found in our paper.

## Software requirement
- python 3
- python-fire
- pysnptools
- numpy
- pandas
- scipy
- plink 1.9

## References
Kangcheng Hou\*, Kathryn S. Burch\*, Arunabha Majumdar, Huwenbo Shi, Nicholas Mancuso, Yue Wu, Sriram Sankararaman & Bogdan Pasaniuc. "Accurate estimation of SNP-heritability from biobank-scale data irrespective of genetic architecture" [[PDF]](https://rdcu.be/bMhni)

## Contact
- Kangcheng Hou: kangchenghou at gmail dot com