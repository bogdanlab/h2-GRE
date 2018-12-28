from pysnptools.snpreader import Bed
import numpy as np
import pandas as pd
import fire
from os.path import join
from scipy import linalg
import json

class Simulation(object): 
    
    def __init__(self, config):
        with open(config) as f:
            config = json.load(f)
        # options
        self.chr_i = config['option']['chr']
        self.chunk_size = config['option']['chunk_size']
        self.num_sim = config['option']['num_sim']
        
        # paths
        self.paths = {'bfile': config['path']['bfile'],
                      'mean_std': config['path']['mean_std'],
                      'beta': config['path']['beta'],
                      'phe': config['path']['phe'],
                      'phe_g': config['path']['phe_g'],
                      'zsc': config['path']['zsc'],
                      'ldak_weights': config['path']['ldak_weights']}
        
        # variables
        self.genotype = Bed(self.paths['bfile'].format(self.chr_i), count_A1=False)
        self.num_indv = self.genotype.row_count
        self.num_snps = self.genotype.col_count

        # parameter settings
        self.params = config['params']

    # utility function
    def read_maf(self):
        mean_std = pd.read_table(self.paths['mean_std'].format(self.chr_i))
        maf = 1 - mean_std['mean'].values / 2
        return maf
            
    def read_ldak_weights(self):
        weights = pd.read_table(self.paths['ldak_weights'], sep=' ')
        weights = weights[weights.chr == self.chr_i].Weight.values
        return weights

    # step1
    def mean_std(self):
        mean = np.zeros(self.num_snps)
        std = np.zeros(self.num_snps)

        for i in range(0, self.num_snps, self.chunk_size):
            sub_geno = self.genotype[:, i : i + self.chunk_size].read().val
            sub_mean = np.nanmean(sub_geno, axis=0)
            mean[i : i + self.chunk_size] = sub_mean
            nanidx = np.where(np.isnan(sub_geno))
            sub_geno[nanidx] = sub_mean[nanidx[1]]
            std[i : i + self.chunk_size] = np.std(sub_geno, axis=0)

        df = pd.DataFrame({'mean': mean, 'std': std})
        df.to_csv(self.paths['mean_std'].format(self.chr_i), sep='\t', index=False)

    # step2 
    def beta(self, param_i):

        betas = np.zeros([self.num_snps, self.num_sim])
        # for `typical` simulation setting, we need `maf` and `ldak_weights`.
        maf = self.read_maf()
        ldak_weights = self.read_ldak_weights()
        param = self.params[param_i]
        maf_range = param['maf_range']
        ncau = int(self.num_snps * param['cau_ratio'])
        
        # genome-wide coefficients
        maf_coef = (maf * (1 - maf)) ** param['maf_exponent']
        ld_coef = ldak_weights ** param['ld_exponent']
        maf_ld_coef = np.sqrt(maf_coef * ld_coef)
        
        for i in range(self.num_sim):
            candidate_list = np.where((maf >= maf_range[0]) & (maf <= maf_range[1]))[0]
            cau = sorted(np.random.choice(candidate_list, size=ncau, replace=False))
            betas[cau, i] = np.random.normal(loc=0.0, scale=1.0, size=len(cau)) * maf_ld_coef[cau]
        np.save(self.paths['beta'].format(param_i, self.chr_i), betas)


    # step3
    def phe_g(self, param_i):
        
        mean_std = pd.read_table(self.paths['mean_std'].format(self.chr_i))
        beta = np.load(self.paths['beta'].format(param_i, self.chr_i))
        phe_g = np.zeros([self.num_indv, self.num_sim])
        
        for i in range(0, self.num_snps, self.chunk_size):
            sub_geno = self.genotype[:, i : i + self.chunk_size].read().val
            sub_mean = mean_std['mean'][i : i + self.chunk_size].values
            sub_std = mean_std['std'][i : i + self.chunk_size].values
            nanidx = np.where(np.isnan(sub_geno))
            sub_geno[nanidx] = sub_mean[nanidx[1]]
            std_sub_geno = np.nan_to_num((sub_geno - sub_mean) / sub_std)
            phe_g += np.dot(std_sub_geno, beta[i : i + self.chunk_size, :]) 
        
        np.save(self.paths['phe_g'].format(param_i, self.chr_i), phe_g)
    
    # step4
    def phe(self, param_i):
        phe_g = np.load(self.paths['phe_g'].format(param_i, self.chr_i))

        # normalize phe_g and add the environment component
        phe_e = np.zeros_like(phe_g)
        hsq = self.params[param_i]['hsq']
        # normalize to get the phe
        for i in range(self.num_sim):
            phe_g[:, i] = phe_g[:, i] * np.sqrt(hsq / np.var(phe_g[:, i]))
            phe_e[:, i] = np.random.normal(loc=0.0, scale=1.0, size=self.num_indv)
            phe_e[:, i] = phe_e[:, i] * np.sqrt((1 - hsq) / np.var(phe_e[:, i]))
        phe = phe_g + phe_e
        np.save(self.paths['phe'].format(param_i, self.chr_i), phe)

    # step5
    def zsc(self, param_i):
        
        phe = np.load(self.paths['phe'].format(param_i, self.chr_i))
        mean_std = pd.read_table(self.paths['mean_std'].format(self.chr_i))
        zsc = np.zeros([self.num_snps, self.num_sim])
        
        for i in range(0, self.num_snps, self.chunk_size):
            sub_geno = self.genotype[:, i : i + self.chunk_size].read().val
            sub_mean = mean_std['mean'][i : i + self.chunk_size].values
            sub_std = mean_std['std'][i : i + self.chunk_size].values
            nanidx = np.where(np.isnan(sub_geno))
            sub_geno[nanidx] = sub_mean[nanidx[1]]
            std_sub_geno = np.nan_to_num((sub_geno - sub_mean) / sub_std)
            zsc[i : i + self.chunk_size, :] = np.dot(std_sub_geno.T, phe) / np.sqrt(self.num_indv)
        # write zsc
        np.save(self.paths['zsc'].format(param_i, self.chr_i), zsc)

if __name__ == '__main__':
    fire.Fire(Simulation)
