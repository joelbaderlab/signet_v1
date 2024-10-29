import matplotlib.pyplot as plt
import pandas as pd
import pdb
import numpy as np
import os
import logging
import argparse

logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger('comparison_signet_pops_tiff')
logger.setLevel(logging.INFO)

def main():
    parser = argparse.ArgumentParser(description='compare genes selected and scored by signet and pops')
    parser.add_argument('--resultsdir', help='full path to results directory', required=True)
    parser.add_argument('--summaryfile', help='file with combined results of signet and pops', required=True)
    args = parser.parse_args()
    
    resultsdir = args.resultsdir
    summaryfile = args.summaryfile
    summary_df = pd.read_csv(summaryfile, sep='\t')
    
    gwas_df = summary_df[~summary_df['locus_name'].str.startswith('loc')]
    
    pvalues = gwas_df['gene_prob_selected'].values
    expvalues = gwas_df['exp_score'].values
    
    colors=['limegreen', 'cornflowerblue']
    label = ['Selection frequency', 'Gene weight']
    
    binwidth2 = 0.1
    bins2 = np.arange(0.0, 1.01, binwidth2)
    
    plt.hist([pvalues, expvalues], bins2, label = label, color=colors)
    plt.legend(loc='upper right')
    plt.xlabel('Value')
    plt.ylabel('Number of genes')
    pdf_filename = os.path.join(resultsdir, 'fig7-comparison.pdf')
    tiff_filename = os.path.join(resultsdir, 'fig7-comparison.tiff')
    plt.savefig(pdf_filename)
    plt.savefig(tiff_filename, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close()

    logger.info(f'{summaryfile} -> {pdf_filename} {tiff_filename}')

    return None


if __name__ == "__main__":
    main()
