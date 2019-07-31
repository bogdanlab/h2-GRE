## Estimating heritability using GRE
After we get the inverse of the LD, we are ready to apply GRE to the OLS sumstats data. The GRE method requires that the sumstats and LD computed using the same genotype reference panel.

The sumstats should be in [LDSC sumstats format](https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format).

Then we can apply the GRE estimator to the sumstats.
```bash
legend=/path/to/bim_file
ld_dir=./ld
sumstats=/path/to/sumstats
python gre.py estimate --legend=${legend} --ld-dir=${ld_dir} --sumstats=${sumstats}
```

In the above we demonstrate applying GRE to the sumstats on one chromosome. To get the genome-wide \(h^2_{\text{GRE}}\) estimate, we can add up the estimates and the variance in each chromosome.

