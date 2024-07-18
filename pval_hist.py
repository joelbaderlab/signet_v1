import matplotlib.pyplot as plt
import pandas as pd
import pdb
import numpy as np
import os
import logging
import argparse

logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger('gene_selection_simple')
logger.setLevel(logging.INFO)

parser = argparse.ArgumentParser(description='Select genes')
parser.add_argument('--results_dir', help='full path to results directory', required=True)
args = parser.parse_args()

results_dir = args.results_dir
filename = os.path.join(results_dir,  'summary_genestats_manuscript.txt')
summary_df = pd.read_csv(filename, sep='\t')
gwas_df = summary_df[~summary_df['locus_name'].str.startswith('loc')]

pvalues = gwas_df['gene_prob_selected'].values
expvalues = gwas_df['exp_score'].values

binwidth1 = 0.05
bins1 = np.arange(0.0, 1.01, binwidth1)

binwidth2 = 0.1
bins2 = np.arange(0.0, 1.01, binwidth2)

colors=['limegreen', 'cornflowerblue']
label = ['Selection frequency', 'Gene weight']


plt.hist([pvalues, expvalues], bins2, label = label, color=colors)
plt.legend(loc='upper right')
plt.xlabel('Value')
plt.ylabel('Number of genes')
plot_filename = os.path.join(results_dir,  'pval_hist_all.pdf')
plt.savefig(plot_filename)
#plt.show()
plt.close()

logger.info(f'{filename} -> {plot_filename}')



