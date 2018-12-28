# Simulating Genome-wide association study

## Configuration
Configure the simulating settings through `config.json`
- `path/genotype`: The genotype bed file for 22 chromosomes.
- `path/ldak_weights`: A single file generated from LDAK software containing the LDAK weights for each SNP.
- `meta/num_sim`: The number of simulations
- `meta/chunk_size`: The chunk of SNP when performing the simulation. In our own simulation, we set this to 250 when using machines with 8G memory.

## Simulation
After configuring, type `bash run_all.sh` to run the simulation. The first time typing this command will first submit the jobs for computing the mean and standard deviation. After it's done, type `bash run_all.sh` again to simulate the GWAS.
