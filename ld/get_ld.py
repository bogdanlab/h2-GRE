from pysnptools.snpreader import Bed
import numpy as np
import pandas as pd
import fire
from os.path import join
from scipy import linalg

def mean_std(bfile, chunk_size, out):
    
    genotype = Bed(bfile, count_A1=False)

    row_count = genotype.row_count
    col_count = genotype.col_count
    mean = np.zeros(col_count)
    std = np.zeros(col_count)

    for i in range(0, col_count, chunk_size):
        sub_geno = genotype[:, i : i + chunk_size].read().val
        sub_mean = np.nanmean(sub_geno, axis=0)
        mean[i : i + chunk_size] = sub_mean
        nanidx = np.where(np.isnan(sub_geno))
        sub_geno[nanidx] = sub_mean[nanidx[1]]
        std[i : i + chunk_size] = np.std(sub_geno, axis=0)

    df = pd.DataFrame({'mean': mean, 'std': std})
    df.to_csv(out, sep='\t', index=False)

# calculate the LD for one block

def std_geno(geno, mean, std):
    nanidx = np.where(np.isnan(geno))
    geno[nanidx] = mean[nanidx[1]]
    geno_std = np.nan_to_num((geno - mean) / std)
    return geno_std

def part(bfile, mean_std_file, row_i, col_i, chunk_size, out_dir):

    geno = Bed(bfile, count_A1=False)
    mean_std = pd.read_table(mean_std_file)

    row_start = row_i * chunk_size
    row_end = (row_i + 1) * chunk_size
    col_start = col_i * chunk_size
    col_end = (col_i + 1) * chunk_size
    print '[{}:{}, {}:{}]'.format(row_start, row_end, col_start, col_end)
    std_row_geno = std_geno(geno[:, row_start : row_end].read().val, mean_std['mean'][row_start : row_end].values, mean_std['std'][row_start : row_end].values)
    std_col_geno = std_geno(geno[:, col_start : col_end].read().val, mean_std['mean'][col_start : col_end].values, mean_std['std'][col_start : col_end].values)

    assert(std_row_geno.shape[0] == std_col_geno.shape[0])
    cov = np.dot(std_row_geno.T, std_col_geno) / (std_row_geno.shape[0])

    np.save(join(out_dir, 'cov_row{}_col{}.npy'.format(row_i, col_i)), cov)

def merge(snp_num, chunk_size, out_dir):
    chunk_num = int( (snp_num + chunk_size - 1) / chunk_size)
    cov = np.zeros([snp_num, snp_num])
    
    for row_i in range(chunk_num):
        for col_i in range(chunk_num):
            row_start = row_i * chunk_size
            row_end = (row_i + 1) * chunk_size
            col_start = col_i * chunk_size
            col_end = (col_i + 1) * chunk_size
            cov[row_start : row_end, col_start : col_end] = np.load(join(out_dir, 'cov_row{}_col{}.npy'.format(row_i, col_i)))
    
    stddev = np.sqrt(np.diag(cov))
    cov /= stddev[:, None]
    cov /= stddev[None, :]

    eig_w, eig_v = linalg.eigh(cov)
    
    np.save(join(out_dir, 'eig_w.npy'), eig_w)
    np.save(join(out_dir, 'eig_v.npy'), eig_v)

def multiply(in_dir, out_path):

    eig_w = np.load(join(in_dir, 'eig_w.npy'))
    eig_v = np.load(join(in_dir, 'eig_v.npy'))

    print('before cutoff: {}'.format(eig_v.shape[0]))
    cond = 1e6 * np.finfo('d').eps
    above_cutoff = (abs(eig_w) > cond * np.max(abs(eig_w)))
    psigma_diag = 1.0 / eig_w[above_cutoff]
    if np.all(above_cutoff) == False:
        print('rank deficent')
        eig_v = eig_v[:, above_cutoff]
    print('after cutoff: {}'.format(eig_v.shape[0]))

    num_snps = eig_v.shape[0]
    inv_cov = np.zeros_like(eig_v)
    chunk_size=5000
    for i in range(0, num_snps, chunk_size):
        for j in range(0, num_snps, chunk_size):
            print i,j
            inv_cov[i : i + chunk_size, j : j + chunk_size] = \
                np.dot(eig_v[i : i + chunk_size, :] * psigma_diag, eig_v[j : j + chunk_size, :].T)
    # save
    np.save(out_path + '.npy', inv_cov)
    with open(out_path + '.txt', 'w') as f:
        f.write(str(len(psigma_diag)))

if __name__ == '__main__':
    fire.Fire()

