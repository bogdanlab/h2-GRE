from pysnptools.snpreader import Bed
import numpy as np
import pandas as pd
import fire
from os.path import join
from scipy import linalg
import json

def mean_std(config, chr_i):
    with open(config) as f:
        config = json.load(f)
    
    bfile = config['path']['genotype'].format(chr_i)
    chunk_size = config['meta']['chunk_size']
    mean_std_path = config['path']['mean_std'].format(chr_i)

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
    df.to_csv(mean_std_path, sep='\t', index=False)

def beta(config, param_i):

    # utility function
    def read_maf(config):
        maf = []
        for chr_i in range(1, 23):
            mean_std = pd.read_table(config['path']['mean_std'].format(chr_i))
            maf.append(1 - mean_std['mean'].values / 2)
        return np.concatenate(maf)
            
    def read_ldak_weights(config):
        weights = pd.read_table(config['path']['ldak_weights'], sep=' ')
        return weights.Weight.values

    with open(config) as f:
        config = json.load(f)

    snps_num_list = []
    for chr_i in range(1, 23):
        bfile = config['path']['genotype'].format(chr_i)
        snps_num_list.append(Bed(bfile, count_A1=False).col_count)
    num_snps = np.sum(snps_num_list)
    num_sim = config['meta']['num_sim']
    betas = np.zeros([num_snps, num_sim])
    
    # for `typical` simulation setting, we need `maf` and `ldak_weights`.
    maf = read_maf(config)
    ldak_weights = read_ldak_weights(config)

    param = config["params"][param_i]
    maf_range = param['maf_range']
    ncau = int(num_snps * param['cau_ratio'])
    
    # genome-wide coefficients
    maf_coef = (maf * (1 - maf)) ** param['maf_exponent']
    ld_coef = ldak_weights ** param['ld_exponent']
    maf_ld_coef = np.sqrt(maf_coef * ld_coef)
    
    
    for i in range(num_sim):
        candidate_list = np.where((maf >= maf_range[0]) & (maf <= maf_range[1]))[0]
        cau = sorted(np.random.choice(candidate_list, size=ncau, replace=False))
        betas[cau, i] = np.random.normal(loc=0.0, scale=1.0, size=len(cau)) * maf_ld_coef[cau]
    
    # output the betas per chromosome
    betas_split = np.vsplit(betas, np.cumsum(snps_num_list))[:-1]
    for chr_i in range(1, 23):
        np.save(config['path']['beta'].format(param_i, chr_i), betas_split[chr_i - 1])

def phe_g(config, param_i, chr_i):
    with open(config) as f:
        config = json.load(f)

    genotype = Bed(config['path']['genotype'].format(chr_i), count_A1=False)
    mean_std = pd.read_table(config['path']['mean_std'].format(chr_i))
    beta = np.load(config['path']['beta'].format(param_i, chr_i))
    num_indv = genotype.row_count
    num_snps = genotype.col_count
    num_sim = config['meta']['num_sim']
    chunk_size = config['meta']["chunk_size"]
    phe_g = np.zeros([num_indv, num_sim])
    
    for i in range(0, num_snps, chunk_size):
        sub_geno = genotype[:, i : i + chunk_size].read().val
        sub_mean = mean_std['mean'][i : i + chunk_size].values
        sub_std = mean_std['std'][i : i + chunk_size].values
        nanidx = np.where(np.isnan(sub_geno))
        sub_geno[nanidx] = sub_mean[nanidx[1]]
        std_sub_geno = np.nan_to_num((sub_geno - sub_mean) / sub_std)
        phe_g += np.dot(std_sub_geno, beta[i : i + chunk_size, :]) 
    
    np.save(config['path']['phe_g'].format(param_i, chr_i), phe_g)

def phe(config, param_i):
    # add up the contribution from 22 chromosome
    # first get the shape of the phe 
    with open(config) as f:
        config = json.load(f)
    total_phe_g = np.load(config['path']['phe_g'].format(param_i, 1))
    total_phe_g = np.zeros_like(total_phe_g)
    phe_e = np.zeros_like(total_phe_g)
    num_indv = total_phe_g.shape[0]
    num_sim = total_phe_g.shape[1]
    for chr_i in range(1, 23):
        phe_g = np.load(config['path']['phe_g'].format(param_i, chr_i))
        total_phe_g = total_phe_g + phe_g
    hsq = config['params'][param_i]['hsq']
    # normalize to get the phe
    for i in range(num_sim):
        total_phe_g[:, i] = total_phe_g[:, i] * np.sqrt(hsq / np.var(total_phe_g[:, i]))
        phe_e[:, i] = np.random.normal(loc=0.0, scale=1.0, size=num_indv)
        phe_e[:, i] = phe_e[:, i] * np.sqrt((1 - hsq) / np.var(phe_e[:, i]))
    phe = total_phe_g + phe_e
    np.save(config['path']['phe'].format(param_i), phe)

def zsc(config, param_i, chr_i):
    with open(config) as f:
        config = json.load(f)
    phe = np.load(config['path']['phe'].format(param_i))
    chunk_size = config['meta']['chunk_size']
    mean_std = pd.read_table(config['path']['mean_std'].format(chr_i))
    genotype = Bed(config['path']['genotype'].format(chr_i), count_A1=False)
    num_indv = genotype.row_count
    num_snps = genotype.col_count
    num_sim = config['meta']['num_sim']
    zsc = np.zeros([num_snps, num_sim])
    
    for i in range(0, num_snps, chunk_size):
        sub_geno = genotype[:, i : i + chunk_size].read().val
        sub_mean = mean_std['mean'][i : i + chunk_size].values
        sub_std = mean_std['std'][i : i + chunk_size].values
        nanidx = np.where(np.isnan(sub_geno))
        sub_geno[nanidx] = sub_mean[nanidx[1]]
        std_sub_geno = np.nan_to_num((sub_geno - sub_mean) / sub_std)
        zsc[i : i + chunk_size, :] = np.dot(std_sub_geno.T, phe) / np.sqrt(num_indv)
    # write zsc
    print zsc.shape
    np.save(config['path']['zsc'].format(param_i, chr_i), zsc)

if __name__ == '__main__':
    fire.Fire()

