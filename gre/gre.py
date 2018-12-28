from pysnptools.snpreader import Bed
import numpy as np
import pandas as pd
import fire
from os.path import join
from scipy import linalg
import logging

def mean_std(bfile, chunk_size, ld_dir):

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
    df.to_csv(join(ld_dir, 'mean_std.txt'), sep='\t', index=False)

def std_geno(geno, mean, std):
    nanidx = np.where(np.isnan(geno))
    geno[nanidx] = mean[nanidx[1]]
    geno_std = np.nan_to_num((geno - mean) / std)
    return geno_std

def cut_ld(bfile, chunk_size, ld_dir):
    geno = Bed(bfile, count_A1=False)
    snp_num = geno.col_count
    # number of chunks snps devides into
    chunk_num = int(np.ceil(float(snp_num) / chunk_size))
    # number of parts of the covariance matrix that needs to be computed
    part_num = chunk_num * chunk_num
    part_list = []
    row_list = []
    col_list = []
    for part_i in range(1, part_num + 1):
        row_i = int((part_i - 1) / chunk_num)
        col_i = int((part_i - 1) % chunk_num)
        row_start = row_i * chunk_size
        row_end = (row_i + 1) * chunk_size
        col_start = col_i * chunk_size
        col_end = (col_i + 1) * chunk_size
        part_list.append(part_i)
        row_list.append('{}-{}'.format(row_start, row_end))
        col_list.append('{}-{}'.format(col_start, col_end))
    df = pd.DataFrame({'row': row_list, 'col': col_list}, columns=['row', 'col'])
    df.to_csv(join(ld_dir, 'part.info'), index=False, header=False, sep='\t')

def calc_ld(bfile, part_i, ld_dir):
    # INFO: part_i counting from 1
    geno = Bed(bfile, count_A1=False)
    mean_std = pd.read_table(join(ld_dir, 'mean_std.txt'))
    part_info = pd.read_table(join(ld_dir, 'part.info'), header=None, sep='\t')
    part_info.columns = ['row', 'col']

    row_start, row_end = [int(i) for i in part_info['row'][part_i - 1].split('-')]
    col_start, col_end = [int(i) for i in part_info['col'][part_i - 1].split('-')]

    std_row_geno = std_geno(geno[:, row_start : row_end].read().val, mean_std['mean'][row_start : row_end].values, mean_std['std'][row_start : row_end].values)
    std_col_geno = std_geno(geno[:, col_start : col_end].read().val, mean_std['mean'][col_start : col_end].values, mean_std['std'][col_start : col_end].values)

    cov = np.dot(std_row_geno.T, std_col_geno) / (std_row_geno.shape[0])
    np.save(join(ld_dir, 'part_{}.npy'.format(part_i)), cov)

def merge_ld(bfile, ld_dir):

    geno = Bed(bfile, count_A1=False)
    snp_num = geno.col_count
    cov = np.zeros([snp_num, snp_num])
    part_info = pd.read_table(join(ld_dir, 'part.info'), header=None, sep='\t', names=['row', 'col'])

    for part_i, part in part_info.iterrows():
        row_start, row_end = [int(i) for i in part_info['row'][part_i].split('-')]
        col_start, col_end = [int(i) for i in part_info['col'][part_i].split('-')]
        cov[row_start : row_end, col_start : col_end] = np.load(join(ld_dir, 'part_{}.npy'.format(part_i + 1)))

    stddev = np.sqrt(np.diag(cov))
    cov /= stddev[:, None]
    cov /= stddev[None, :]
    inv_cov, rank = linalg.pinvh(cov, return_rank=True)

    np.save(join(ld_dir, 'inv_ld.npy'), inv_cov)
    with open(join(ld_dir, 'rank.txt'), 'w') as f:
        f.write(str(rank))

def estimate(legend, ld_dir, sumstats):
    # 1. read legend & check if match
    # 2. estimate
    legend = pd.read_table(legend, header=None, sep='\t', names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
    sumstats = pd.read_table(sumstats)
    SNP_same = (legend['SNP'] == sumstats['SNP']).all()
    A1_same = (legend['A1'] == sumstats['A1']).all()
    A2_same = (legend['A2'] == sumstats['A2']).all()
    if not (SNP_same & A1_same & A2_same):
        print('Reference panel of sumstats and genotype does not match')
        return

    # if match, start to estimate
    inv_ld = np.load(join(ld_dir, 'inv_ld.npy'))
    rank = np.loadtxt(join(ld_dir, 'rank.txt'))
    num_indv = np.mean(sumstats.N.values)

    beta_gwas = sumstats.Z.values / np.sqrt(num_indv)
    quad_form_func = (lambda x, A : np.dot(np.dot(x.T, A), x))
    quad_form = quad_form_func(beta_gwas, inv_ld)
    est = (num_indv * quad_form - rank) / (num_indv - rank)
    var = ((num_indv / (num_indv - rank)) ** 2) * \
            (2 * rank * (1 - est) / num_indv + 4 * est) * \
            (1 - est) / num_indv
    print('h2-GRE estimate: {:.6f}, standard deviation: {:.6f}'.format(est, np.sqrt(var)))

if __name__ == '__main__':
    fire.Fire()
