# Input data

## Requirement on #SNPs and #individuals
In order for the desirable properties of the h2-GRE estimator to hold, h2-GRE requires individual-level genotype and phenotype data where the number of individuals (N) is larger than the number of genotyped SNPs per chromosome (m_k for the k-th chromosome). 

As a reference point, in our analyses of the UK Biobank genotyped SNPs (UK Biobank Axiom Array), chromosome 1 (the largest chromosomes) had about 45-55K SNPs depending on how we QC-ed the data, while number of individuals is about 300K.

## Prepare the data
The genotype and phenotypes should be in [PLINK](https://www.cog-genomics.org/plink/2.0/input) format. We will also use PLINK software to perform OLS regression.

!!! note 
    Here we assume that you know the basic about plink, e.g. basic commands, data format. You can read more on the plink website about [covariate file](http://www.cog-genomics.org/plink/1.9/input#covar), [phenotype file](http://www.cog-genomics.org/plink/1.9/input#pheno) and [association analysis](http://www.cog-genomics.org/plink/1.9/assoc).

We will use the following data.

- `all_bfile`: Genotype for the **unrelated** individuals in UK Biobank **after** quality control. This contains genome-wide data in one file.
- `pheno_file`: the phenotype file.
- `covar_file`: the covariates which will be used for association studies.
