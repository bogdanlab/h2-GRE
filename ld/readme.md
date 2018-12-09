# Computing the covariance matrix

## Configuration
Configuration through `metadata.json`.

`mean_std_chunk_size`: Chunk size for SNPswhen computing the mean and standardization
`ld_chunk_size`: Chunk size for SNPs when estimating the LD matrix.

## Steps
1. Computing the mean and standard deviation for each chromosome. `bash run_step1_mean_std.sh`
2. Computing covariance matrix of each block. `bash run_step2_block.sh`
3. Merge the covariance matrix from blocks and perform singular value decomposition. `bash run_step3_merge.sh`
4. Calculate the inverse of the matrix using the matrix computed from singular value decomposition. `bash run_step4_multiply.sh`

