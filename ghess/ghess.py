import numpy as np
from pysnptools.snpreader import Bed
import fire
from os.path import join
import pandas as pd
import json

# proof of concept code for g-HESS
# only use it with the accompanying simulation code
def get_quad_form(zsc, inv_ld, num_indv):
    def quad_form(x, A):
        return np.dot(np.dot(x.T, A), x)
    beta_gwas = zsc / np.sqrt(num_indv)
    return quad_form(beta_gwas, inv_ld)

def ghess_estimator(quad_form, num_indv, rank_list):
    # return total SNP-heritability and variance
    # hsq estimation
    quad_form = np.array(quad_form)
    rank_list = np.array(rank_list)
    local_hsq = (num_indv * quad_form - rank_list) / (num_indv - rank_list)
    total_hsq = np.sum(local_hsq)
    # hsq variance
    var_part1 = (num_indv / (num_indv - rank_list)) ** 2
    var_part2 = (2 * rank_list * (1 - local_hsq) / num_indv + 4 * local_hsq)
    var_part3 = (1 - local_hsq) / num_indv
    total_hsq_var = np.sum(var_part1 * var_part2 * var_part3)
    return local_hsq, total_hsq, total_hsq_var

def estimate(inv_ld_dir, gwas_dir, out_dir):
    with open(join(gwas_dir, 'config.json')) as f:
        config = json.load(f)
    
    num_sim = config['meta']['num_sim']
    genotype = Bed(config['path']['genotype'].format(1), count_A1=False)
    num_indv = genotype.row_count
    params = config['params']
    
    rank_list = []
    # use quad_form_list[param_i][sim_i] to access
    quad_form_list = [[[] for sim_i in range(num_sim)] for param_i in range(len(params))]
    
    for chr_i in range(1, 23):
        inv_ld = None # clear memory
        inv_ld = np.load(join(inv_ld_dir, '{}.npy'.format(chr_i)))
        rank_list.append(np.loadtxt(join(inv_ld_dir, '{}.txt'.format(chr_i))))
        for param_i in range(len(params)):
            # zsc: num_snps x num_sim
            zsc = np.load(join(gwas_dir, config['path']['zsc'].format(param_i, chr_i)))
            for sim_i in range(num_sim):
                quad_form_list[param_i][sim_i].append(get_quad_form(zsc[:, sim_i], inv_ld, num_indv))

    # now compute the total heritability estimation and analytical variance
    for param_i, param in enumerate(params):
        print(param)
        for sim_i in range(num_sim):
            local_hsq_estimate, total_hsq_estimate, hsq_variance = ghess_estimator(quad_form_list[param_i][sim_i], num_indv, rank_list)
            with open(join(out_dir, 'param{}_sim{}.txt'.format(param_i, sim_i)), 'w') as f:
                f.write('local heritability estimation (one for each chromosome): ' + ' '.join([str(i) for i in local_hsq_estimate]) + '\n')
                f.write('total heritability estimation: {}\n'.format(total_hsq_estimate))
                f.write('total heritability variance: {}\n'.format(hsq_variance))
                
if __name__ == '__main__':
    fire.Fire()

