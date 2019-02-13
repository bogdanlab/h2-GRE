import pandas as pd
import fire


def sumstats_plink2ldsc(plink_sumstats, legend, out):
    sumstats = {}
    legend = pd.read_csv(legend, header=None, sep='\t')
    with open(plink_sumstats) as f:
        lines = f.readlines()[1:]
        lines = [line.strip().split() for line in lines]
        sumstats['SNP'] = [line[1] for line in lines]
        sumstats['N'] = [int(line[5]) for line in lines]
        sumstats['Z'] = [float(line[6]) / float(line[7]) for line in lines]
        sumstats['A1'] = [line[3] for line in lines]
        sumstats['A2'] = legend[5]
    sumstats = pd.DataFrame(sumstats, columns=['SNP', 'N', 'Z', 'A1', 'A2'])
    sumstats.to_csv(out, sep='\t', index=False)

if __name__ == '__main__':
    fire.Fire()
