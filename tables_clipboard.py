import pandas as pd
import numpy as np
import csv
import os
import pdb
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
from collections import defaultdict
import glob as glob
import utils_plot
import pickle


results_dir = './results/TopMed_QT+HR+QRS+PR+JT_250kb_NewDist'
data_dir = './data/TopMed'
seedless_base_filename = 'SIGNET'


def write_finalstats_selected(results_dir):
    filename = os.path.join(results_dir,  'summary_genestats_manuscript.txt')
    summary_df = pd.read_csv(filename, sep='\t')
    gwas_df = summary_df[~summary_df['locus_name'].str.startswith('loc')]
    n_genes = len(gwas_df)

    selected_df = gwas_df.loc[gwas_df['gene_signet'] == 'signet'] #take locus once
    n_loci = len(selected_df)
    locus_width = selected_df['locus_width'].values
    median_locus_width = np.median(locus_width)
    ngenes_locus = selected_df['locus_ngene'].values
    median_ngenes_locus = np.median(ngenes_locus)


    loci_stats = {'n_genes': n_genes,
                  'n_loci': n_loci, 'median_locus_width': median_locus_width, 'median_ngenes_locus': median_ngenes_locus}

    df = pd.DataFrame([loci_stats])
    filename = os.path.join(results_dir, 'selected_dist_stats.csv')

    df.to_csv(filename, index=False, header=True)



def write_feature_weight_dict_final(results_dir):
    feature_weights = defaultdict(list)
    feature_list = ['distfac_active', 'omim', 'exome', 'coloc']
    feature_finalvalues_df = pd.DataFrame()

    for feature in feature_list:
        feature_weights[feature] = []
        for file in glob.glob(results_dir + '/seed*' + feature + '_log.txt'):
            df = pd.read_csv(file, delimiter='\t')
            # ADDED the [-1] to only consider final value
            feature_weights[feature].append(df[feature].values[-1])

        feature_finalvalues_df[feature] = feature_weights[feature]

    file = os.path.join(
        results_dir, 'features_weights_final.txt')

    feature_finalvalues_df.to_csv(file, index=False, sep='\t')

    for feature in feature_list:
        print('feature: ', feature, 'mean: ', np.mean(
            feature_weights[feature]), 'std: ', np.std(feature_weights[feature]))

    feature_summary_df = pd.DataFrame(columns=['feature', 'mean', 'std'])

    for feature in feature_list:
        df = {'feature': feature, 'mean': np.mean(
            feature_weights[feature]), 'std': np.std(feature_weights[feature])}
        feature_summary_df = feature_summary_df._append(df, ignore_index=True)

    file = os.path.join(results_dir, 'features_weights_final_summary.txt')
    feature_summary_df.to_csv(file, index=False, sep='\t')

def write_feature_weight_dict_init(results_dir): #used for init rand, otherwise init is all the same, but here is the code
    feature_weights = defaultdict(list)
    feature_list = ['distfac_active', 'omim', 'exome', 'coloc']
    feature_finalvalues_df = pd.DataFrame()

    for feature in feature_list:
        feature_weights[feature] = []
        for file in glob.glob(results_dir + '/seed*' + feature + '_log.txt'):
            df = pd.read_csv(file, delimiter='\t')
            # ADDED the [0] to only consider init value
            feature_weights[feature].append(df[feature].values[0])

        feature_finalvalues_df[feature] = feature_weights[feature]

    file = os.path.join(
        results_dir, 'features_weights_init.txt')
    feature_finalvalues_df.to_csv(file, index=False, sep='\t')

    for feature in feature_list:
        print('feature: ', feature, 'mean: ', np.mean(
            feature_weights[feature]), 'std: ', np.std(feature_weights[feature]))

    feature_summary_df = pd.DataFrame(columns=['feature', 'mean', 'std'])

    for feature in feature_list:
        df = {'feature': feature, 'mean': np.mean(
            feature_weights[feature]), 'std': np.std(feature_weights[feature])}
        feature_summary_df = feature_summary_df._append(df, ignore_index=True)

    file = os.path.join(results_dir, 'features_weights_init_summary.txt')
    feature_summary_df.to_csv(file, index=False, sep='\t')




def kegg_pathway_comparison(file1path, file1name, file2path, file2name, results_dir, seedless_base_filename):

    gene_pathwaysdf_file1 = pd.read_csv(file1path, sep='\t')
    gene_pathwaysdf_file2 = pd.read_csv(file2path, sep='\t')

    pathways_df = gene_pathwaysdf_file1.merge(
        right=gene_pathwaysdf_file2, on="Term", validate='one_to_one', suffixes=('_' + file1name, '_' + file2name))

    pathways_df[('-log10(Adjsuted_pval)_' + file1name)] = - \
        np.log10(pathways_df[('Adjusted P-value_' + file1name)])
    pathways_df[('-log10(Adjsuted_pval)_' + file2name)] = - \
        np.log10(pathways_df[('Adjusted P-value_' + file2name)])

    pathways_df['pval_diff'] = pathways_df['-log10(Adjsuted_pval)_' + file1name] - \
        pathways_df['-log10(Adjsuted_pval)_' + file2name]

    pathways_df[file1name + '_genes'] = [sorted(x.split(';'))
                                         for x in pathways_df['Genes_' + file1name].values]
    pathways_df[file2name + '_genes'] = [sorted(x.split(';'))
                                         for x in pathways_df['Genes_' + file2name].values]

    pathways_df['CommonGenes'] = sorted([set(pathways_df[file1name + '_genes'][x]).intersection(
        set(pathways_df[file2name + '_genes'][x])) for x in range(len(pathways_df))])

    pathways_df[file1name + '_only'] = [sorted(set(pathways_df[file1name + '_genes'][x]) - set(
        pathways_df[file2name + '_genes'][x])) for x in range(len(pathways_df))]

    pathways_df[file2name + '_only'] = [sorted(set(pathways_df[file2name + '_genes'][x]) -
                                               set(pathways_df[file1name + '_genes'][x])) for x in range(len(pathways_df))]

    file = os.path.join(results_dir, (seedless_base_filename + '_' +
                                      file1name + '_vs_' + file2name + '_pathways.txt'))
    pathways_df.to_csv(file, index=False, sep='\t')


def count_npass(results_dir):
    feature = 'distfac_active' #doesn't matter which feature we pick; just need the number
    npass_list = []
    for file in glob.glob(results_dir + '/seed*' + feature + '_log.txt'):
        df = pd.read_csv(file, delimiter='\t')
        npass_list.append(len(df))
    print(np.median(npass_list))
    pdb.set_trace()


##########
# TABLE 2
##########
if 0:
    write_finalstats_selected(results_dir)
    pdb.set_trace()

##########
# TABLE 5
##########
if 0:
    write_feature_weight_dict_final(results_dir)
    write_feature_weight_dict_init(results_dir)
    pdb.set_trace()


##########
# Table 6
##########
if 0:
    results_dir = './results/TopMed_QT+HR+QRS+PR+JT_250kb_NewDist'

    summary_file = os.path.join(results_dir, 'summary_genestats_manuscript.txt')
    summary_df = pd.read_csv(summary_file, sep='\t')

    omim_loci = set()
    exome_loci = set()
    coloc_loci = set()
    poorsig_loci = set()

    omim_loci_omim_selected = set()
    omim_loci_exome_selected = set()
    omim_loci_coloc_selected = set()
    omim_loci_poorsig_selected = set()
    omim_loci_mindist_selected = set()

    exome_loci_omim_selected = set()
    exome_loci_exome_selected = set()
    exome_loci_coloc_selected = set()
    exome_loci_poorsig_selected = set()
    exome_loci_mindist_selected = set()

    coloc_loci_omim_selected = set()
    coloc_loci_exome_selected = set()
    coloc_loci_coloc_selected = set()
    coloc_loci_poorsig_selected = set()
    coloc_loci_mindist_selected = set()

    poorsig_loci_omim_selected = set()
    poorsig_loci_exome_selected = set()
    poorsig_loci_coloc_selected = set()
    poorsig_loci_poorsig_selected = set()
    poorsig_loci_mindist_selected = set()

    locus_df = pd.DataFrame()
    for locus in summary_df['locus_name'].unique():
        if locus.startswith('loc'):
            continue

        locus_df = summary_df[summary_df['locus_name'] == locus].copy()
        locus_special_str = "".join(locus_df['gene_special'])
        locus_geneset = set(locus_df['gene_name'].values)

        selected_gene = locus_df[locus_df['gene_signet']
                                 == 'signet']['gene_name'].values[0]

        selected_gene_special_str = "".join(
            locus_df.loc[locus_df['gene_name'] == selected_gene]['gene_special'])
        locus_mindist_gene = locus_df[locus_df['gene_mindist']
                                      != '.']['gene_name'].values[0]

        if 'omim' in locus_special_str:
            omim_loci.add(locus)

            # now keep count of the type of selected gene:

            if 'omim' in selected_gene_special_str:
                omim_loci_omim_selected.add(locus)
            elif 'exome' in selected_gene_special_str:  # should be 0
                omim_loci_exome_selected.add(locus)
            elif 'coloc' in selected_gene_special_str:
                omim_loci_coloc_selected.add(locus)
            elif selected_gene == locus_mindist_gene:
                omim_loci_mindist_selected.add(locus)
            else:
                omim_loci_poorsig_selected.add(locus)

            continue

        if 'exome' in locus_special_str:
            exome_loci.add(locus)
            if 'omim' in selected_gene_special_str:
                exome_loci_omim_selected.add(locus)
            elif 'exome' in selected_gene_special_str:
                exome_loci_exome_selected.add(locus)
            elif 'coloc' in selected_gene_special_str:
                exome_loci_coloc_selected.add(locus)
            elif selected_gene == locus_mindist_gene:
                exome_loci_mindist_selected.add(locus)
            else:
                exome_loci_poorsig_selected.add(locus)

            continue
        if 'coloc' in locus_special_str:
            coloc_loci.add(locus)
            if 'omim' in selected_gene_special_str:
                coloc_loci_omim_selected.add(locus)
            elif 'exome' in selected_gene_special_str:
                coloc_loci_exome_selected.add(locus)
            elif 'coloc' in selected_gene_special_str:
                coloc_loci_coloc_selected.add(locus)
            elif selected_gene == locus_mindist_gene:
                coloc_loci_mindist_selected.add(locus)
            else:
                coloc_loci_poorsig_selected.add(locus)

            continue

        else:
            poorsig_loci.add(locus)
            if 'omim' in selected_gene_special_str:
                poorsig_loci_omim_selected.add(locus)
            elif 'exome' in selected_gene_special_str:
                poorsig_loci_exome_selected.add(locus)
            elif 'coloc' in selected_gene_special_str:
                poorsig_loci_coloc_selected.add(locus)
            elif selected_gene == locus_mindist_gene:
                poorsig_loci_mindist_selected.add(locus)
            else:
                poorsig_loci_poorsig_selected.add(locus)

            continue


    print('omim_loci:  %d' % len(omim_loci))
    print(len(omim_loci_omim_selected))
    print(len(omim_loci_exome_selected))
    print(len(omim_loci_coloc_selected))
    print(len(omim_loci_mindist_selected))
    print(len(omim_loci_poorsig_selected))


    print('exome_loci: %d' % len(exome_loci))
    print(len(exome_loci_omim_selected))
    print(len(exome_loci_exome_selected))
    print(len(exome_loci_coloc_selected))
    print(len(exome_loci_mindist_selected))
    print(len(exome_loci_poorsig_selected))


    print('coloc_loci: %d' % len(coloc_loci))
    print(len(coloc_loci_omim_selected))
    print(len(coloc_loci_exome_selected))
    print(len(coloc_loci_coloc_selected))
    print(len(coloc_loci_mindist_selected))
    print(len(coloc_loci_poorsig_selected))


    print('poorsig_loci: %d' % len(poorsig_loci))
    print(len(poorsig_loci_omim_selected))
    print(len(poorsig_loci_exome_selected))
    print(len(poorsig_loci_coloc_selected))
    print(len(poorsig_loci_mindist_selected))
    print(len(poorsig_loci_poorsig_selected))


    omim_df = summary_df.loc[summary_df['locus_name'].isin(omim_loci)]
    omim_file = os.path.join(results_dir, 'summary_genestats_omimloci.txt')
    omim_df.to_csv(omim_file, index=False, sep='\t')

    exome_df = summary_df.loc[summary_df['locus_name'].isin(exome_loci)]
    exome_file = os.path.join(results_dir, 'summary_genestats_exomeloci.txt')
    exome_df.to_csv(exome_file, index=False, sep='\t')

    coloc_df = summary_df.loc[summary_df['locus_name'].isin(coloc_loci)]
    coloc_file = os.path.join(results_dir, 'summary_genestats_colocloci.txt')
    coloc_df.to_csv(coloc_file, index=False, sep='\t')

    poorsig_df = summary_df.loc[summary_df['locus_name'].isin(poorsig_loci)]
    poorsig_file = os.path.join(results_dir, 'summary_genestats_poorsigloci.txt')
    poorsig_df.to_csv(poorsig_file, index=False, sep='\t')


    pdb.set_trace()


    ##########
    # Table 7
    ##########
    #check summary_df filepath top of Table6 script
    hiomim_geneset = set()
    hiomim_selected = set()
    hiomim_notselected = set()
    hiomim_gene_selected_omim = set()
    hiomim_gene_selected_exome = set()
    hiomim_gene_selected_coloc = set()
    hiomim_gene_selected_mindist = set()
    hiomim_gene_selected_poorsig = set()

    hiexome_geneset = set()
    hiexome_selected = set()
    hiexome_notselected = set()
    hiexome_gene_selected_omim = set()
    hiexome_gene_selected_exome = set()
    hiexome_gene_selected_coloc = set()
    hiexome_gene_selected_mindist = set()
    hiexome_gene_selected_poorsig = set()

    hicoloc_geneset = set()
    hicoloc_selected = set()
    hicoloc_notselected = set()
    hicoloc_gene_selected_omim = set()
    hicoloc_gene_selected_exome = set()
    hicoloc_gene_selected_coloc = set()
    hicoloc_gene_selected_mindist = set()
    hicoloc_gene_selected_poorsig = set()

    himindist_geneset = set()
    himindist_selected = set()
    himindist_notselected = set()
    himindist_gene_selected_omim = set()
    himindist_gene_selected_exome = set()
    himindist_gene_selected_coloc = set()
    himindist_gene_selected_mindist = set()
    himindist_gene_selected_poorsig = set()

    hipoorsig_geneset = set()
    hipoorsig_selected = set()
    hipoorsig_notselected = set()
    hipoorsig_gene_selected_omim = set()
    hipoorsig_gene_selected_exome = set()
    hipoorsig_gene_selected_coloc = set()
    hipoorsig_gene_selected_mindist = set()
    hipoorsig_gene_selected_poorsig = set()


    locus_df = pd.DataFrame()
    for locus in summary_df['locus_name'].unique():
        if locus.startswith('loc'):
            continue
        locus_df = summary_df[summary_df['locus_name'] == locus].copy()

        locus_geneset = set(locus_df['gene_name'].values)
        selected_gene = locus_df[locus_df['gene_signet']
                                 == 'signet']['gene_name'].values[0]
        selected_gene_special_str = "".join(
            locus_df.loc[locus_df['gene_name'] == selected_gene]['gene_special'])

        for gene in locus_geneset:
            gene_special_str = "".join(
                locus_df.loc[locus_df['gene_name'] == gene]['gene_special'])
            locus_mindist_gene = locus_df[locus_df['gene_mindist']
                                          != '.']['gene_name'].values[0]

            if 'omim' in gene_special_str:
                hiomim_geneset.add(gene)
                if gene == selected_gene:
                    hiomim_selected.add(gene)
                else:
                    hiomim_notselected.add(gene)
                    if 'omim' in selected_gene_special_str:
                        hiomim_gene_selected_omim.add(gene)
                    elif 'exome' in selected_gene_special_str:
                        hiomim_gene_selected_exome.add(gene)
                    elif 'coloc' in selected_gene_special_str:
                        hiomim_gene_selected_coloc.add(gene)
                    elif selected_gene == locus_mindist_gene:
                        hiomim_gene_selected_mindist.add(gene)
                    else:
                        hiomim_gene_selected_poorsig.add(gene)

            elif 'exome' in gene_special_str:
                hiexome_geneset.add(gene)
                if gene == selected_gene:
                    hiexome_selected.add(gene)
                else:
                    hiexome_notselected.add(gene)
                    if 'omim' in selected_gene_special_str:
                        hiexome_gene_selected_omim.add(gene)
                    elif 'exome' in selected_gene_special_str:
                        hiexome_gene_selected_exome.add(gene)
                    elif 'coloc' in selected_gene_special_str:
                        hiexome_gene_selected_coloc.add(gene)
                    elif selected_gene == locus_mindist_gene:
                        hiexome_gene_selected_mindist.add(gene)
                    else:
                        hiexome_gene_selected_poorsig.add(gene)

            elif 'coloc' in gene_special_str:
                hicoloc_geneset.add(gene)
                if gene == selected_gene:
                    hicoloc_selected.add(gene)
                else:
                    hicoloc_notselected.add(gene)
                    if 'omim' in selected_gene_special_str:
                        hicoloc_gene_selected_omim.add(gene)
                    elif 'exome' in selected_gene_special_str:
                        hicoloc_gene_selected_exome.add(gene)
                    elif 'coloc' in selected_gene_special_str:
                        hicoloc_gene_selected_coloc.add(gene)
                    elif selected_gene == locus_mindist_gene:
                        hicoloc_gene_selected_mindist.add(gene)
                    else:
                        hicoloc_gene_selected_poorsig.add(gene)

            elif gene == locus_mindist_gene:  # selected_gene == locus_mindist_gene:
                himindist_geneset.add(gene)
                if gene == selected_gene:
                    himindist_selected.add(gene)
                else:
                    himindist_notselected.add(gene)
                    if 'omim' in selected_gene_special_str:
                        himindist_gene_selected_omim.add(gene)
                    elif 'exome' in selected_gene_special_str:
                        himindist_gene_selected_exome.add(gene)
                    elif 'coloc' in selected_gene_special_str:
                        himindist_gene_selected_coloc.add(gene)
                    elif selected_gene == locus_mindist_gene:
                        himindist_gene_selected_mindist.add(gene)
                    else:
                        himindist_gene_selected_poorsig.add(gene)

            else:
                hipoorsig_geneset.add(gene)
                if gene == selected_gene:
                    hipoorsig_selected.add(gene)
                else:
                    hipoorsig_notselected.add(gene)
                    if 'omim' in selected_gene_special_str:
                        hipoorsig_gene_selected_omim.add(gene)
                    elif 'exome' in selected_gene_special_str:
                        hipoorsig_gene_selected_exome.add(gene)
                    elif 'coloc' in selected_gene_special_str:
                        hipoorsig_gene_selected_coloc.add(gene)
                    elif selected_gene == locus_mindist_gene:
                        hipoorsig_gene_selected_mindist.add(gene)
                    else:
                        hipoorsig_gene_selected_poorsig.add(gene)


    print(len(hiomim_geneset))
    print(len(hiomim_selected))
    print(len(hiomim_gene_selected_omim))
    print(len(hiomim_gene_selected_exome))
    print(len(hiomim_gene_selected_coloc))
    print(len(hiomim_gene_selected_mindist))
    print(len(hiomim_gene_selected_poorsig))

    print('exome')
    print(len(hiexome_geneset))
    print(len(hiexome_selected))
    print(len(hiexome_gene_selected_omim))
    print(len(hiexome_gene_selected_exome))
    print(len(hiexome_gene_selected_coloc))
    print(len(hiexome_gene_selected_mindist))
    print(len(hiexome_gene_selected_poorsig))

    print('coloc')
    print(len(hicoloc_geneset))
    print(len(hicoloc_selected))
    print(len(hicoloc_gene_selected_omim))
    print(len(hicoloc_gene_selected_exome))
    print(len(hicoloc_gene_selected_coloc))
    print(len(hicoloc_gene_selected_mindist))
    print(len(hicoloc_gene_selected_poorsig))

    print('mindist')
    print(len(himindist_geneset))
    print(len(himindist_selected))
    print(len(himindist_gene_selected_omim))
    print(len(himindist_gene_selected_exome))
    print(len(himindist_gene_selected_coloc))
    print(len(himindist_gene_selected_mindist))
    print(len(himindist_gene_selected_poorsig))

    print('poorsig')
    print(len(hipoorsig_geneset))
    print(len(hipoorsig_selected))
    print(len(hipoorsig_gene_selected_omim))
    print(len(hipoorsig_gene_selected_exome))
    print(len(hipoorsig_gene_selected_coloc))
    print(len(hipoorsig_gene_selected_mindist))
    print(len(hipoorsig_gene_selected_poorsig))

    pdb.set_trace()

##########
# TABLE 8
##########
if 0:
    import counts_postprocessing_manuscript
    #list of the special gene from loci with 1 special genes
    filename = os.path.join(results_dir,  'summary_genestats_manuscript.txt')
    summary_df = pd.read_csv(filename, sep='\t')
    gwas_df = summary_df[~summary_df['locus_name'].str.startswith('loc')]
    specialgene_loci_onespecial = gwas_df.loc[(gwas_df['locus_nspecialgene'] == 1) & (gwas_df['gene_special'] != '.')]['gene_name'].values
    #print(sorted(specialgene_loci_onespecial))
    #run_blockedspecial.sh

    with open(results_dir + '/gene_special.pickle', 'rb') as handle:
        gene_special = pickle.load(handle)

    blockedgene_df = pd.DataFrame(columns=['gene_name', 'gene_special', 'gene_recovered'])

    #summary_file for the blocked genes:
    top_results_dir = './results/blocked_genes'
    for blocked_gene in sorted(specialgene_loci_onespecial):
        blocked_gene_dir = os.path.join(top_results_dir, '_'.join(['blocked', blocked_gene]))
        blocked_dir_results = os.path.join(blocked_gene_dir, 'TopMed_QT+HR+QRS+PR+JT_250kb_NewDist')

        if 0: #need to run once
            seed_start = 0
            seed_end = 100
            counts_postprocessing_manuscript.write_aggregated_statsfiles(seedless_base_filename, blocked_dir_results, seed_start = seed_start, seed_end=seed_end)

        blockedgene_summary_file = os.path.join(blocked_dir_results, 'summary_genestats_manuscript.txt')
        blockedgene_summary_df = pd.read_csv(blockedgene_summary_file, sep='\t')

        blockedgene_recovery = blockedgene_summary_df.loc[blockedgene_summary_df['gene_name']==blocked_gene, 'gene_signet'].values[0]

        df = {'gene_name': blocked_gene, 'gene_special': gene_special[blocked_gene], 'gene_recovered': blockedgene_recovery}
        blockedgene_df  = blockedgene_df._append(df, ignore_index = True)

    result_filename = os.path.join(top_results_dir, 'blockedgene_summary.txt')
    blockedgene_df.to_csv(result_filename, index=False, header=True, sep='\t')

    pdb.set_trace()



##########
# TABLE 9
##########
if 0:
    #pathway analysis
    file1path = os.path.join(results_dir, 'KEGG_2021_signet.txt')
    file2path = os.path.join(results_dir,  'KEGG_2021_mindist.txt')
    file3path = os.path.join(results_dir,  'KEGG_2021_signetplus.txt')
    file1name = 'signet'
    file2name = 'mindist'
    file3name = 'signetplus'
    kegg_pathway_comparison(file1path, file1name, file2path,
                            file2name, results_dir, seedless_base_filename)
    kegg_pathway_comparison(file1path, file1name, file3path,
                            file3name, results_dir, seedless_base_filename)

    pdb.set_trace()



##########
# number passes to convergence
##########
if 0:
    results_dir = './results/TopMed_QT+HR+QRS+PR+JT_250kb_NewDist'
    count_npass(results_dir)
    pdb.set_trace()



##########
# compare rand-init with bestguess-init
##########
if 0:
    results_dir_bestinit = './results/TopMed_QT+HR+QRS+PR+JT_250kb_NewDist'
    filename_bestinit = os.path.join(results_dir_bestinit,  'summary_genestats_manuscript.txt')
    summary_df_bestinit = pd.read_csv(filename_bestinit, sep='\t')[['locus_name', 'gene_name', 'gene_prob_selected']]

    results_dir_randinit = './results/TopMed_QT+HR+QRS+PR+JT_250kb_initrandTrue_NewDist'
    filename_randinit = os.path.join(results_dir_randinit,  'summary_genestats_manuscript.txt')
    summary_df_randinit = pd.read_csv(filename_randinit, sep='\t')[['locus_name', 'gene_name', 'gene_prob_selected']]
    new_df = pd.merge(summary_df_bestinit, summary_df_randinit, on=['locus_name','gene_name'], suffixes = ('_bestinit', '_randinit'))

    new_file = os.path.join(results_dir_bestinit, 'init_best_rand.txt')
    new_df.to_csv(new_file, sep='\t')


    pdb.set_trace()
