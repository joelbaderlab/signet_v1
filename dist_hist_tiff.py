import pandas as pd
import numpy as np
import csv
import os
import pdb
import matplotlib.pyplot as plt
import utils_plot
import pickle
import logging
import argparse

logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger('gene_selection_simple')
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Distance histogram')
parser.add_argument('--results_dir', help='full path to results directory', required=True)
args = parser.parse_args()

def get_gwas_distance_hist(summary_file, init_file, snp_file, dist_file):
    logger.info(f'{summary_file} + {init_file} + {snp_file} -> {dist_file}')

    summary_df = pd.read_csv(summary_file, sep='\t')
    gwas_df = summary_df[~summary_df['locus_name'].str.startswith('loc')].copy()

    # the problem is that the distance in the summary genestats is abs(distance)
    # so ... have to look up the SNP, get the distance, and then recalculate

    snp_bp_df = pd.read_csv(snp_file)
    print(f'{snp_file} -> {snp_bp_df.shape}')
    snp2bp = dict()
    for (snp, bp) in zip(snp_bp_df['rsid'], snp_bp_df['rsid_start']):
        snp2bp[snp] = bp

    signeddist_list = [ ]
    snp_pos = [ ]
    for (gene_tss, gene_strand, snp) in zip(gwas_df['gene_tss'], gwas_df['gene_strand'], gwas_df['gene_closestsnp']):
        assert(snp in snp2bp), f'could not find snp {snp}'
        snpbp = snp2bp[snp]
        snp_pos.append(snpbp)
        mydist = int(snpbp) - int(gene_tss)
        if gene_strand == '-':
            mydist = -mydist
        mydist = float(mydist) / 1000.0
        signeddist_list.append(mydist)

    gwas_df['gene_signedlocusdist_kbp'] = signeddist_list
    gwas_df['snp_pos'] = snp_pos

    init_df = pd.read_csv(init_file, sep=',', index_col=False)
    init_df.columns = ['col1', 'col2']
    init_df.rename(columns={'col1': 'locus_name', 'col2': 'gene_initactive'}, inplace=True)
    init_df['gene_initactive'] = init_df['gene_initactive'].str.strip("{}' ")
    init_df = init_df[~init_df['locus_name'].str.startswith('loc')]
    init_activegenes = set(init_df['gene_initactive'].values)
    gwas_df['gene_initactive'] = '.'
    gwas_df.loc[gwas_df['gene_name'].isin(init_activegenes), 'gene_init_active'] = 'init_active'
    gwas_df = gwas_df[['locus_name', 'locus_chr', 'gene_name', 'gene_strand', 'gene_closestsnp', 'gene_tss', 'snp_pos', 'gene_locusdist_kbp', 'gene_signedlocusdist_kbp', 'gene_init_active', 'gene_mindist', 'gene_signet']]
    dist_file = os.path.join(results_dir, 'gwas_distance_hist.txt')
    gwas_df.to_csv(dist_file, sep='\t')

    logger.info(f'done writing {dist_file}')
    return None

# start of main program

results_dir = args.results_dir

logger.info(f'creating distance histogram table')
summary_file = os.path.join(results_dir, 'summary_genestats_manuscript.txt')
init_file = os.path.join(results_dir, 'init_activeset.csv')
snp_file = os.path.join(results_dir, 'snp_bp.csv')
dist_file = os.path.join(results_dir, 'gwas_distance_hist.txt')
get_gwas_distance_hist(summary_file, init_file, snp_file, dist_file)

gwas_df = pd.read_csv(dist_file, sep='\t')
signet_dist_values = gwas_df.loc[gwas_df['gene_signet']=='signet', 'gene_signedlocusdist_kbp']
mindist_dist_values = gwas_df.loc[gwas_df['gene_mindist']=='mindist', 'gene_signedlocusdist_kbp']
init_dist_values = gwas_df.loc[gwas_df['gene_init_active']=='init_active', 'gene_signedlocusdist_kbp']

bin_width = 20
minval = 250
bins1 = np.arange(-minval, minval+1, bin_width) # 0 is the center of the middle bin, -10 to 10
# bins2 = np.arange(-260, 261, bin_width) # 0 is an edge, bins are at -20,0 and 0,20
bins = bins1
densebins = np.arange(-minval,minval+1,0.1)

colors=['limegreen', 'cornflowerblue', 'sandybrown']
label = ['SigNet', 'BestGuess', 'MinDist']
values_array =  np.array([signet_dist_values, init_dist_values, mindist_dist_values])
values_array = np.transpose(values_array)

distance = densebins
param = 161.3
exp_model = (0.5/param) * np.exp(-abs(distance)/param)
total_counts = len(signet_dist_values)
exp_model_scaled = exp_model * bin_width * total_counts

plt.hist(values_array , bins, histtype='bar', color=colors, label=label, stacked=False, fill=True, density=False)
plt.plot(distance, exp_model_scaled, 'k', linewidth=1,  label='Learned distribution')
plt.legend(loc='upper right')
plt.xlabel('Distance from SNP to TSS, kb')
plt.ylabel('Number of genes')

mypdf = os.path.join(results_dir, 'fig3-dist-hist.pdf')
mytiff = os.path.join(results_dir, 'fig3-dist-hist.tiff')
plt.savefig(mypdf)
plt.savefig(mytiff, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
plt.close()

logger.info(f'done writing {mypdf} {mytiff}')

