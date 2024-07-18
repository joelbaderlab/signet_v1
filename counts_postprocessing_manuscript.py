import argparse
import yaml
import logging
import pandas as pd
from collections import defaultdict
import numpy as np
from scipy.stats import gmean
import os as os
import glob as glob
import csv
import pdb
import pickle
import matplotlib.pyplot as plt
import utils_ppi
import utils_data
import utils_plot
import random


logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger('gene_selection_simple')
logger.setLevel(logging.INFO)


def get_args():
    # ,epilog='small call: see run.sh', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser = argparse.ArgumentParser(description='Select genes')
    parser.add_argument('--results_dir', help='full path to results directory', required=True)
    parser.add_argument('--seed_start', help='integer bottom of python range, usually 0', type=int, required=True)
    parser.add_argument('--seed_end', help='integer top of python range, last seed is seed_end-1', type=int, required=True)
    args = parser.parse_args()
    return(args)




def finalnetworkscores(results_dir, seed_start, seed_end):
    df = pd.DataFrame(columns=['seed', 'finalnetworkscore'])

    for seed in range(seed_start, seed_end):
        base_filename = f'seed{seed:02d}'

        file = os.path.join(results_dir, base_filename +
                            '_networkscore_log.txt')

        seed_df = pd.read_csv(file, delimiter='\t')
        seed = seed_df['seed'].values[-1] # [-1] to only consider final value
        seed_finalnetworkscore = seed_df['network_score'].values[-1]
        df.loc[len(df.index)] = [seed, seed_finalnetworkscore]

    max_scores = np.max(df['finalnetworkscore'])
    min_scsores = np.min(df['finalnetworkscore'])
    mean_scsores = np.mean(df['finalnetworkscore'])
    std_scsores = np.std(df['finalnetworkscore'])

    finalscores_file = os.path.join(results_dir, 'final_networkscores.txt')
    df.to_csv(finalscores_file, index=False, sep='\t')

    maxnetworkscore_seed = int(df.loc[df['finalnetworkscore']
                                      == max_scores, 'seed'].values[0])
    return maxnetworkscore_seed



def write_aggregated_statsfiles(seedless_base_filename, results_dir, seed_start, seed_end):

    maxnetworkscore_seed = finalnetworkscores(results_dir, seed_start=seed_start, seed_end=seed_end)


    base_filename = f'seed{maxnetworkscore_seed:02d}'

    genestats_file = os.path.join(results_dir, base_filename + '_extended_finalpass_scores.txt')

    genestats_df = pd.read_csv(genestats_file, sep='\t', dtype=str)

    genestats_df = genestats_df[['locus_name', 'gene_name', 'gene_pattern', 'gene_mindist',
                                 'gene_locusdist_kbp', 'distance_score', 'gene_special', 'gene_nonzero', 'special_score', 'ppi_score', 'exp_score', 'ppi_deg', 'ppi_deg_active_cnt', 'ppi_deg_active_names', 'TF_score', 'TF_indeg', 'TF_in_active_cnt', 'TF_in_active_names', 'TF_outdeg', 'TF_out_active_cnt', 'TF_out_active_names']]

    counts_filename = os.path.join(
        results_dir, seedless_base_filename + '_counts_allruns_prob.txt')

    counts_df = pd.read_csv(counts_filename, sep='\t', dtype=str)
    counts_df = counts_df[['locus_name', 'gene_name', 'gene_prob_selected']]
    final_df = pd.merge(counts_df, genestats_df, on=[
                        'locus_name', 'gene_name'], validate='one_to_one')

    gene_closest_filename = os.path.join(results_dir, 'gene_closest.txt')
    gene_df = pd.read_csv(gene_closest_filename, sep='\t', dtype=str)
    gene_df = gene_df[['gene_name', 'gene_closestsnp']]
    final_df = final_df.merge(gene_df, how='left', on='gene_name')

    locus_to_genes_filename = os.path.join(results_dir, 'locus_to_genes.txt')
    locus_df = pd.read_csv(locus_to_genes_filename, sep='\t', dtype=str)

    final_df = final_df.merge(locus_df, how='left', on=[
                              'locus_name', 'gene_name', 'gene_locusdist_kbp'])


    # add a column for locus_mindist_gene

    final_df['gene_signet'] = '.'
    final_df['gene_signet+'] = '.'
    for locus in final_df['locus_name'].unique():

        locus_df = final_df[final_df['locus_name'] == locus].copy()

        mindist_genestr = locus_df.loc[locus_df['gene_mindist']
                                       == 'mindist']['gene_name'].values[0]

        locus_nspecialgene = sum(locus_df['gene_special'] != '.')

        locus_ninteractiongene = sum(locus_df['ppi_deg'].astype(
            int) + locus_df['TF_indeg'].astype(int) + locus_df['TF_outdeg'].astype(int) != 0)

        max_prob = np.max(locus_df['gene_prob_selected'])
        highest_genes = locus_df.loc[locus_df['gene_prob_selected']
                                     == max_prob]['gene_name'].values

        selected_gene = random.choice(sorted(highest_genes))

        final_df.loc[final_df['locus_name'] == locus,
                     'locus_mindistgene'] = mindist_genestr
        final_df.loc[final_df['locus_name'] == locus,
                     'locus_nspecialgene'] = locus_nspecialgene
        final_df.loc[final_df['locus_name'] == locus,
                     'locus_ninteractiongene'] = locus_ninteractiongene

        final_df.loc[(final_df['locus_name'] == locus) & (
            final_df['gene_name'] == selected_gene), 'gene_signet'] = 'signet'
        final_df.loc[(final_df['locus_name'] == locus) & ((final_df['gene_name'] == selected_gene) | (
            final_df['gene_special'] != '.')), 'gene_signet+'] = 'signet+'


    final_df['gene_interesting'] = '.'
    final_df.loc[(final_df['gene_nonzero'] == 'nonzero') | (final_df['gene_mindist'] == 'mindist') | (
        final_df['gene_special'] != '.'), 'gene_interesting'] = 'interesting'

    final_df = final_df.astype({'locus_chr': 'int', 'locus_start': 'int',
                                'locus_end': 'int', 'gene_start': 'int', 'gene_end': 'int', 'gene_tss': 'int', 'gene_strand': 'str'})


    final_df = final_df.sort_values(
        by=['locus_chr', 'locus_start', 'locus_end', 'gene_tss', 'gene_end'])


    summary_df = final_df

    summary_df['locus_hisig_type'] = '.'
    summary_df['gene_hisig_type'] = '.'
    locus_df = pd.DataFrame()
    for locus in summary_df['locus_name'].unique():
        if locus.startswith('loc'):
            continue
        locus_df = summary_df[summary_df['locus_name'] == locus].copy()
        locus_special_str = "".join(locus_df['gene_special'])
        if 'omim' in locus_special_str:
            summary_df.loc[summary_df['locus_name']
                           == locus, 'locus_hisig_type'] = 'omim'
        elif 'exome' in locus_special_str:
            summary_df.loc[summary_df['locus_name'] ==
                           locus, 'locus_hisig_type'] = 'exome'
        elif 'coloc' in locus_special_str:
            summary_df.loc[summary_df['locus_name'] ==
                           locus, 'locus_hisig_type'] = 'coloc'
        else:
            summary_df.loc[summary_df['locus_name'] ==
                           locus, 'locus_hisig_type'] = 'poorsig'

        locus_geneset = set(locus_df['gene_name'].values)
        locus_mindist_gene = locus_df[locus_df['gene_mindist']
                                      != '.']['gene_name'].values[0]

        for gene in locus_geneset:
            gene_special_str = "".join(
                locus_df.loc[locus_df['gene_name'] == gene]['gene_special'])

            if 'omim' in gene_special_str:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'omim'
            elif 'exome' in gene_special_str:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'exome'
            elif 'coloc' in gene_special_str:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'coloc'
            elif gene == locus_mindist_gene:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'mindist'
            else:
                summary_df.loc[(summary_df['locus_name'] == locus) & (
                    summary_df['gene_name'] == gene), 'gene_hisig_type'] = 'poorsig'

    summary_df = summary_df[['locus_chr', 'locus_start', 'locus_end', 'locus_name', 'locus_ngene', 'locus_width', 'locus_pheno', 'locus_mindistgene', 'locus_nspecialgene', 'locus_ninteractiongene', 'locus_hisig_type', 'gene_start', 'gene_end', 'gene_tss', 'gene_strand', 'gene_locusdist_kbp', 'gene_closestsnp', 'gene_name', 'gene_type', 'gene_prob_selected',
                             'gene_signet', 'gene_mindist', 'gene_special',  'gene_signet+', 'distance_score', 'special_score', 'ppi_score', 'TF_score', 'exp_score', 'gene_hisig_type', 'ppi_deg', 'ppi_deg_active_cnt', 'ppi_deg_active_names', 'TF_indeg', 'TF_in_active_cnt', 'TF_in_active_names', 'TF_outdeg', 'TF_out_active_cnt', 'TF_out_active_names']]

    summary_df_filename = os.path.join(
        results_dir, 'summary_genestats_manuscript.txt')
    summary_df.to_csv(summary_df_filename, index=False, sep='\t')

def write_feature_weight_finalval_CV(CV_results_dir):
    feature_weights = defaultdict(list)
    feature_list = ['distfac_active', 'omim', 'exome', 'coloc']
    cvseed_feature_finalval_df = pd.DataFrame()

    for cv_seed in range(100, 199, 1):
        small_df = pd.DataFrame()
        results_dir = os.path.join(CV_results_dir, '_'.join([project_dirname, phenostr, flankstr, 'NewDist', 'CV', str(cv_seed)]))

        for feature in feature_list:
            df = pd.DataFrame()
            for file in glob.glob(results_dir + '/seed*' + feature + '_log.txt'):
                df = pd.read_csv(file, delimiter='\t')
                feature_weights[feature].append(df[feature].values[-1])

            small_df[feature] = feature_weights[feature]

        small_df['cv_seed'] = np.repeat( cv_seed, len(feature_weights[feature]))

        if cvseed_feature_finalval_df.empty:
            cvseed_feature_finalval_df = small_df
        else:
            cvseed_feature_finalval_df = pd.concat([cvseed_feature_finalval_df, small_df], ignore_index=True)

    file = os.path.join(CV_results_dir, 'CV_features_weights_final.txt')
    cvseed_feature_finalval_df.to_csv(file, index=False, sep='\t')

    for feature in feature_list:
        print('feature: ', feature, 'mean: ', np.mean(
            feature_weights[feature]), 'std: ', np.std(feature_weights[feature]))

    feature_summary_df = pd.DataFrame(columns=['feature', 'mean', 'std'])


    for feature in feature_list:
        df2 = {'feature': feature, 'mean': np.mean(
            feature_weights[feature]), 'std': np.std(feature_weights[feature])}
        feature_summary_df.loc[len(feature_summary_df)] = df2


    file = os.path.join(CV_results_dir, 'CV_features_weights_final_summary.txt')
    feature_summary_df.to_csv(file, index=False, sep='\t')


# start of main program
    
args = get_args()
logger.info('args %s', str(args))
seedless_base_filename = 'SIGNET'
results_dir = args.results_dir
seed_start = args.seed_start
seed_end = args.seed_end


if True:
    logger.info(f'calling write_aggregated_statsfile on main results')
    #main signet
    write_aggregated_statsfiles(seedless_base_filename, args.results_dir, args.seed_start, args.seed_end)

# special logic to combine results for shuffled replica networks
if False:
    #shuffled signet
    shuffledtop_results_dir = config.results_dir + '/Shuffled'
    for nx_seed in range(100, 200):
        results_subdir = os.path.join(shuffledtop_results_dir, str(nx_seed))
        seed_start = 0
        seed_end = 100
        write_aggregated_statsfiles(seedless_base_filename, results_subdir, seed_start = seed_start, seed_end=seed_end)
        selected_df = pd.DataFrame()
        summary_filename = 'summary_genestats_manuscript.txt'
        summary_file = os.path.join(results_subdir, summary_filename)
        summary_df = pd.read_csv(summary_file, sep='\t')
        selected_geneset = sorted(summary_df[summary_df['gene_signet'] == 'signet']['gene_name'].values)
        selected_df = selected_df._append(selected_geneset)
        result_filename = os.path.join(results_subdir, 'signet_genes_nxseed_'+str(nx_seed)+'.txt')
        selected_df.to_csv(result_filename, index=False, header=None, sep='\t')

if False:
    #shuffled signet combine pathway enrich results
    shuffledtop_results_dir = config.results_dir + '/Shuffled'
    enrichr_overlap_df = pd.DataFrame()
    enrichr_pval_df = pd.DataFrame()
    pathways = ['Adrenergic signaling in cardiomyocytes', 'Arrhythmogenic right ventricular cardiomyopathy','Cardiac muscle contraction', 'Cholinergic synapse','Circadian entrainment','Dilated cardiomyopathy','GnRH signaling pathway','Hypertrophic cardiomyopathy','Oxytocin signaling pathway']
    enrichr_overlap_df['Term'] = pathways
    enrichr_pval_df['Term'] = pathways
    for nx_seed in range(100, 200):
        results_subdir = os.path.join(shuffledtop_results_dir, str(nx_seed))
        one_enrichr_filename = 'KEGG_2021_nx' + str(nx_seed) + '.txt'
        one_enrichr_file = os.path.join(results_dir, one_enrichr_filename)
        df = pd.read_csv(one_enrichr_file, sep='\t')
        overlap_filtered_df = df[df['Term'].isin(pathways)][['Term', 'Overlap']]
        overlap_filtered_df['Overlap'] = overlap_filtered_df['Overlap'].str.split('/').str[0]
        overlap_filtered_df = overlap_filtered_df.rename(columns={'Overlap': ('nxseed_'+str(nx_seed))})
        enrichr_overlap_df = enrichr_overlap_df.merge(overlap_filtered_df, on='Term')

        pval_filtered_df = df[df['Term'].isin(pathways)][['Term', 'Adjusted P-value']]
        pval_filtered_df = pval_filtered_df.rename(columns={'Adjusted P-value': ('nxseed_'+str(nx_seed))})
        enrichr_pval_df = enrichr_pval_df.merge(pval_filtered_df, on='Term')


    enrichr_overlap_df = enrichr_overlap_df.iloc[:, 1:].apply(pd.to_numeric)
    enrichr_overlap_df['arith_mean'] = enrichr_overlap_df.iloc[:, 1:].mean(axis=1)
    enrichr_overlap_df['std'] = enrichr_overlap_df.iloc[:, 1:].std(axis=1)
    enrichr_overlap_filename = os.path.join(shuffledtop_results_dir, 'enrichr_overlap_shuffled.txt')
    enrichr_overlap_df.to_csv(enrichr_overlap_filename, index=False, header=True, sep='\t')

    enrichr_pval_df['geom_mean'] = enrichr_pval_df.iloc[:, 1:].apply(lambda x: gmean(x), axis=1)
    enrichr_pval_filename = os.path.join(shuffledtop_results_dir, 'enrichr_pval_shuffled.txt')
    enrichr_pval_df.to_csv(enrichr_pval_filename, index=False, header=True, sep='\t')

if False:
    # 80% data cross-validation
    CV_results_dir = results_dir + '/CV'
    write_feature_weight_finalval_CV(CV_results_dir)


