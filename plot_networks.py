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

# start of main program

results_dir = args.results_dir

tfpickle = os.path.join(results_dir, 'TFsource_target.pickle')
networkpickle = os.path.join(results_dir, 'network_ppi.pickle')
genepickle = os.path.join(results_dir, 'gene_special.pickle')

with open(tfpickle, 'rb') as handle:
    TFsource_target = pickle.load(handle)

with open(networkpickle, 'rb') as handle:
    network_ppi = pickle.load(handle)

with open(genepickle, 'rb') as handle:
    gene_special = pickle.load(handle)

summary_filename = 'summary_genestats_manuscript.txt'
summary_file = os.path.join(results_dir, summary_filename)

summary_df = pd.read_csv(summary_file, sep='\t')
selected_geneset = set(
    summary_df[summary_df['gene_signet'] == 'signet']['gene_name'].values)

for anchor_gene in selected_geneset:
    utils_plot.plot_InteractionNetwork_v1(
        anchor_gene, results_dir, summary_filename, network_ppi, gene_special, TFsource_target)

logger.info(f'done creating network plots')

