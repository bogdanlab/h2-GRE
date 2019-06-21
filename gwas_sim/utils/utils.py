from os.path import join
import yaml
import os.path

def generate_config(template, out_dir):

    def write_config(data, out_dir, cnt):
        """
        helper function to write yaml configuration
        """
        save_dir = join(out_dir, 'model_{}'.format(cnt))
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        with open(join(save_dir, 'config.yaml'), 'w') as f:
            yaml.dump(data, f, default_flow_style=False)

    with open(template, 'r') as f:
        data = yaml.safe_load(f)
    data['bfile'] = '/u/scratch/k/kangchen/gwas_sim/utils/out/10k'
    data['num_part'] = 3
    data['sumher']['ldak_weights'] = '/u/scratch/k/kangchen/sim_gwas_typed_10k/step0_metadata/ldak_weights.txt'
    data['h2g'] = 0.25
    cnt = 0
    for maf in [0.0, 0.75]:
        for ld in [0.0, 1.0]:
            data['sumher']['maf_exponent'] = maf
            data['sumher']['ld_exponent'] = ld
            # vary causal snps position
            data['sumher']['cau_ratio'] = 1.0
            data['sumher']['maf_low'] = 0.0
            data['sumher']['maf_high'] = 1.0
            write_config(data, out_dir, cnt)
            cnt += 1
            for maf_low, maf_high in [(0.01, 0.05), (0.05, 0.5), (0.0, 1.0)]:
                
                data['sumher']['cau_ratio'] = 0.01
                data['sumher']['maf_low'] = maf_low
                data['sumher']['maf_high'] = maf_high
                write_config(data, out_dir, cnt)
                cnt += 1

def ldak_weights_add_chr(ldak_weights, bim, out):
    weights = pd.read_csv(ldak_weights, sep=' ', dtype=str)
    bim = pd.read_csv(bim, sep='\t', header=None)
    weights['chr'] = bim[0]
    weights = weights[['chr'] + weights.columns.tolist()[:-1]]
    weights.to_csv(out, sep=' ', index=False)

def npy2plink(npy_dir, params, fam, plink_dir):
    """
    convert phenotypes from .npy to plink format
    """
    params = pd.read_table(params, sep=' ', skiprows=2, dtype=str)
    fam = pd.read_table(fam, sep=' ', header=None, dtype=str)
    fam = fam[[0,1]]
    params = params.reset_index(drop=True)
    for param_i, param in params.iterrows():
        print('cau_ratio: {}, hsq: {}, maf_exp: {}, ld_exp: {}, low: {}, high: {}'.format(param.cau_ratio, param.hsq, param.maf, param.ld, param.maf_low, param.maf_high))
        phe = np.load(join(npy_dir, 'cau_ratio_{}_hsq_{}_maf_{}_ld_{}_range_{}_{}.npy'.format(param.cau_ratio, param.hsq, param.maf, param.ld, param.maf_low, param.maf_high)))

        num_sim = phe.shape[1]
        for sim_i in range(num_sim):
            phen = pd.DataFrame(fam)
            phen[2] = phe[:, sim_i]
            phen.to_csv(join(plink_dir, 'cau_ratio_{}_hsq_{}_maf_{}_ld_{}_range_{}_{}_sim_{}.phen'.format(param.cau_ratio, param.hsq, param.maf, param.ld, param.maf_low, param.maf_high, sim_i + 1)), index=False, header=None, sep='\t')

if __name__ == '__main__':
    fire.Fire()

