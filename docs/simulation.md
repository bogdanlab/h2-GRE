We also provide code for simulating OLS sumstats with different genetic architecture.

Simulating OLS sumstats can be divided into the following steps:

1. computing the mean and standard deviation for each SNP in population and standardize the genotype matrix
2. draw the effect sizes for each SNP
3. simulating the phenotypes
4. simulating the OLS sumstats


```python
geno = Bed('/path/to/bfile', count_A1=False).read().val
# impute nan SNPs
nanidx = np.where(np.isnan(genotype))
mean = np.nanmean(genotype, axis=0)
genotype[nanidx] = mean_geno[nanidx[1]] 
# standardize genotypes
geno -= geno.mean(axis=0)
geno /= geno.std(axis=0)
num_indv, num_snps = geno.shape
# here we draw from infinitesimal effect sizes model with heritability 0.1 for an example
beta = np.random.normal(loc=0.0, scale=1.0, size=num_snps)
# simulating phenotypes
phe_g = np.dot(geno, beta)
phe_g = phe_g * np.sqrt(hsq / np.var(phe_g))
phe_e = np.random.normal(loc=0.0, scale=1.0, size=num_indv)
phe_e = phe_e * np.sqrt((1 - hsq) / np.var(phe_e))
phe = phe_g + phe_e
# simulating OLS sumstats z-scores
zsc = np.dot(geno.T, phe) / np.sqrt(num_indv)
```

This code should work fine for small dataset. But in some cases, we may want to
1. Simulate using large dataset.
2. Specify how the effect sizes are drawn (different genetic architecture)
   
We can use the provided code located at `sim_gwas/chr` in this repository to achieve these two easily. By Changing `option: chunk_size` in `config.json`, we can fit the computation into memory. By specifies `params`, we can change how the effect sizes are drawn.

In order to simulate using the provided code, we can input the command as follows:
### 1. Mean and standard deviation
```shell 
config=/path/to/config
python ./sim_gwas.py mean-std --config=${config}
```

### 2. Effect sizes
```shell
config=/path/to/config
python ./sim_gwas.py beta --config=${config} --param-i=${param_i}
```

### 3. phenotype from genetic component
```shell
config=/path/to/config
python ./sim_gwas.py phe-g --config=${config} --param-i=${param_i}
```

### 4. phenotype (with environmental component added)
```shell
config=/path/to/config
python ./sim_gwas.py phe --config=${config} --param-i=${param_i}
```

### 5. OLS sumstats
```shell
config=/path/to/config
python ./sim_gwas.py zsc --config=${config} --param-i=${param_i}
```