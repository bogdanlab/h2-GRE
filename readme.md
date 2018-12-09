# g-HESS scripts

Scripts for manuscript "Accurate estimation of SNP-heritability under minimal
assumptions about genetic architecture".

We show how to apply g-HESS to simulation data. This repository mainly contains three parts:
1. Compute the inverse of the genotype covariance matrix. This is needed for g-HESS estimator.
2. Simulate GWAS using genome-wide data.
3. Apply g-HESS to the simulated data.

**The code is still under construction.**

## Requirement
1. Anaconda Python 2
2. pysnpstools
3. python-fire
4. [Hoffman2 Cluster](https://idre.ucla.edu/hoffman2)

## Scripts
1. [Compute the inverse of the genotype covariance matrix](./ld)
2. [Simulate GWAS](./sim_gwas)
3. [Apply g-HESS](./ghess)
