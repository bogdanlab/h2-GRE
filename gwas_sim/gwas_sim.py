import numpy as np
from pysnptools.snpreader import Bed
import pandas as pd
import yaml
from os.path import join
import os
import sys
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

class BaseModel(object):
    
    def __init__(self, root_dir):
        """
            read the config and initialize the model
        """
        # load configuration
        with open(join(root_dir, 'config.yaml')) as f:
            self.config = yaml.safe_load(f)
        np.random.seed(self.config['seed'])
        self.root_dir = root_dir
        self.geno = Bed(join(root_dir, self.config['bfile']), count_A1=False)
        self.num_indv, self.num_snps = self.geno.shape
        self.snp_list = np.linspace(0, self.num_snps, num=self.config['num_part'] + 1).astype(int)
        
    def submit(self, task):
        assert(task in ['mean_std', 'beta', 'phe_g', 'phe', 'zsc'])
        # qsub 
        if task == 'mean_std':
            for i in range(self.config['num_part']):
                self.compute_mean_std(i)
        elif task == 'beta':
            self.compute_beta()
        elif task == 'phe_g':
            for i in range(self.config['num_part']):
                self.compute_phe_g(i)
        elif task == 'phe':
            self.compute_phe()
        elif task == 'zsc':
            for i in range(self.config['num_part']):
                self.compute_zsc(i)
        else:
            raise NotImplementedError
        

    def compute_mean_std(self, part_i):
        """
        Compute the mean and standard deviation for a genotype file
        """
        
        chunk_size = self.config['chunk_size']
        start_snp = self.snp_list[part_i]
        stop_snp = self.snp_list[part_i + 1]
        geno_part = self.geno[:, start_snp : stop_snp]
        mean = np.zeros(geno_part.shape[1])
        std = np.zeros(geno_part.shape[1])
        
        for i in range(0, geno_part.shape[1], chunk_size):
            sub_geno = geno_part[:, i : i + chunk_size].read().val
            sub_mean = np.nanmean(sub_geno, axis=0)
            mean[i : i + chunk_size] = sub_mean
            # fill in the missing value
            nanidx = np.where(np.isnan(sub_geno))
            sub_geno[nanidx] = sub_mean[nanidx[1]]
            std[i : i + chunk_size] = np.std(sub_geno, axis=0)
    
        df = pd.DataFrame({'mean': mean, 'std': std})
        
        save_dir = join(self.root_dir, self.config['mean_std_dir'])
        
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        df.to_csv(join(save_dir, '{}.txt'.format(part_i)), index=False)
        

    def compute_beta(self):
        """
        model: the generative model
        """
        betas = np.zeros([self.num_snps, self.config['num_sim']])
        
        for sim_i in range(self.config['num_sim']):
            betas[:, sim_i] = np.random.normal(loc=0.0, scale=1.0, size=self.num_snps)
        save_dir = join(self.root_dir, self.config['beta_dir'])
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        for part_i in range(self.config['num_part']):
            start_snp = self.snp_list[part_i]
            stop_snp = self.snp_list[part_i + 1]
            betas_part = betas[start_snp : stop_snp, :]
            np.save(join(save_dir, '{}.npy'.format(part_i)), betas_part)


    def compute_phe_g(self, part_i):
        """
            compute the 'unnormalized' genetic component of phenotypes
        """
        betas = np.load(join(self.root_dir, 
                        self.config['beta_dir'], 
                        '{}.npy'.format(part_i)))
        
        mean_std = pd.read_csv(join(self.root_dir,
                        self.config['mean_std_dir'],
                        '{}.txt'.format(part_i)))
        
        mean, std = mean_std['mean'].values, mean_std['std'].values
        chunk_size = self.config['chunk_size']
        start_snp = self.snp_list[part_i]
        stop_snp = self.snp_list[part_i + 1]
        geno_part = self.geno[:, start_snp : stop_snp]
        num_snps = geno_part.shape[1]
        phe_g = np.zeros([self.num_indv, self.config['num_sim']])
        for i in range(0, num_snps, chunk_size):
            sub_geno = geno_part[:, i : i + chunk_size].read().val
            sub_mean = mean[i : i + chunk_size]
            sub_std = std[i : i + chunk_size]
            nanidx = np.where(np.isnan(sub_geno))
            sub_geno[nanidx] = sub_mean[nanidx[1]]
            std_sub_geno = (sub_geno - sub_mean) / sub_std
            
            phe_g += np.dot(std_sub_geno, betas[i : i + chunk_size, :])

        save_dir = join(self.root_dir, self.config['phe_g_dir'])

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        np.save(join(save_dir, '{}.npy'.format(part_i)), phe_g)


    def compute_phe(self):
        phe_g = np.zeros([self.num_indv, self.config['num_sim']])
        for part_i in range(self.config['num_part']):
            phe_g += np.load(join(self.root_dir,
                            self.config['phe_g_dir'],
                            '{}.npy'.format(part_i)))
        phe_e = np.zeros_like(phe_g)

        h2g = self.config['h2g']

        for sim_i in range(self.config['num_sim']):
            phe_g[:, sim_i] = phe_g[:, sim_i] * np.sqrt(h2g / np.var(phe_g[:, sim_i]))
            phe_e[:, sim_i] = np.random.normal(loc=0.0, scale=1.0, size=self.num_indv)
            phe_e[:, sim_i] = phe_e[:, sim_i] * np.sqrt((1 - h2g) / np.var(phe_e[:, sim_i]))
        phe = phe_g + phe_e

        save_dir = join(self.root_dir, self.config['phe_dir'])

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        np.save(join(save_dir, 'phe.npy'), phe)

    
    def compute_zsc(self, part_i):
        phe = np.load(join(self.root_dir, 
                    self.config['phe_dir'],
                    'phe.npy'))
        
        mean_std = pd.read_csv(join(self.root_dir,
                        self.config['mean_std_dir'],
                        '{}.txt'.format(part_i)))
        mean, std = mean_std['mean'].values, mean_std['std'].values
        chunk_size = self.config['chunk_size']
        start_snp = self.snp_list[part_i]
        stop_snp = self.snp_list[part_i + 1]
        geno_part = self.geno[:, start_snp : stop_snp]
        num_snps = geno_part.shape[1]
        zsc = np.zeros([num_snps, self.config['num_sim']])
        for i in range(0, num_snps, chunk_size):
            sub_geno = geno_part[:, i : i + chunk_size].read().val
            sub_mean = mean[i : i + chunk_size]
            sub_std = std[i : i + chunk_size]
            # impute missing data using sample average
            nanidx = np.where(np.isnan(sub_geno))
            sub_geno[nanidx] = sub_mean[nanidx[1]]
            # standardized genotypes of this subpart of causal genotypes
            std_sub_geno = (sub_geno - sub_mean) / sub_std
            zsc[i : i + chunk_size, :] = np.dot(std_sub_geno.T, phe) / np.sqrt(self.num_indv)
        
        save_dir = join(self.root_dir, self.config['zsc_dir'])

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        np.save(join(save_dir, '{}.npy'.format(part_i)), zsc)


class SumherModel(BaseModel):
    def __init__(self, root_dir):
        BaseModel.__init__(self, root_dir)
    def compute_beta(self):
        """
        implment SumHer model
        """
        betas = np.zeros([self.num_snps, self.config['num_sim']])
        
        # read maf
        maf = []
        for part_i in range(self.config['num_part']):
            mean_std = pd.read_csv(join(self.root_dir,
                self.config['mean_std_dir'],
                '{}.txt'.format(part_i)))
            maf.append(1 - mean_std['mean'].values / 2)
        maf = np.concatenate(maf)
        param = self.config['sumher']
        # read ldak weights
        ldak_weights = pd.read_table(param['ldak_weights'], sep=' ').Weight.values
        maf_coef = (maf * (1 - maf)) ** param['maf_exponent']
        ld_coef = ldak_weights ** param['ld_exponent']
        maf_ld_coef = np.sqrt(maf_coef * ld_coef)
        ncau = int(self.num_snps * param['cau_ratio'])


        for sim_i in range(self.config['num_sim']):
            candidate_list = np.where((maf >= param['maf_low']) & (maf <= param['maf_high']))[0]
            cau = sorted(np.random.choice(candidate_list, size=ncau, replace=False))
            betas[cau, sim_i] = np.random.normal(loc=0.0, scale=1.0, size=len(cau)) * maf_ld_coef[cau]
        print(cau)    
        save_dir = join(self.root_dir, self.config['beta_dir'])
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        for part_i in range(self.config['num_part']):
            start_snp = self.snp_list[part_i]
            stop_snp = self.snp_list[part_i + 1]
            betas_part = betas[start_snp : stop_snp, :]
            np.save(join(save_dir, '{}.npy'.format(part_i)), betas_part)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='GWAS simulation')
    parser.add_argument('--root_dir', type=str,
                    help="root directory of simulation model")
    parser.add_argument('--step', type=str, 
                choices=['mean_std', 'beta', 'phe_g', 'phe', 'zsc'],
                help='simulation step')
    parser.add_argument('--part_i', type=int, default=None,
                help='which part')
    
    args = parser.parse_args()
    model = SumherModel(args.root_dir)
    if args.step == 'mean_std':
        model.compute_mean_std(args.part_i)
    elif args.step == 'beta':
        model.compute_beta()
    elif args.step == 'phe_g':
        model.compute_phe_g(args.part_i)
    elif args.step == 'phe':
        model.compute_phe()
    elif args.step == 'zsc':
        model.compute_zsc(args.part_i)
    else:
        raise NotImplementedError
