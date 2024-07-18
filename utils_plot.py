import pandas as pd
import numpy as np
import os as os
import logging
import matplotlib.pyplot as plt
import graphviz as gv
import utils_ppi
import pdb

logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger('gene_selection_simple')
logger.setLevel(logging.INFO)

SKIPPLOT = False

def plot_InteractionNetwork(results_dir, summary_file, network_ppi, gene_special,  TFsource_target, prob_threshold_min=0, prob_threshold_max=1):
    if SKIPPLOT:
        return None
    logger.info('plotting network')

    graph_name =  'network_' + \
        str(prob_threshold_min) + '_maxprob_' + \
        str(prob_threshold_max)
    img_file = os.path.join(results_dir, graph_name)

    # selectedgenes_file = os.path.join(
    #     results_dir, seedless_base_filename + selectedgenes_filename)

    # finalprob_genestats_df = pd.read_csv(selectedgenes_file, sep='\t')

    summary_df = pd.read_csv(summary_file, sep='\t') #'summary_genestats_manuscript.txt'

    selected_df = summary_df.loc[(summary_df['gene_prob_selected'] >= prob_threshold_min) & (
        summary_df['gene_prob_selected'] <= prob_threshold_max)]

    selected_genes = selected_df['gene_name'].values

   # mindist_genes = set()
    # for locus in locus_mindsistgeneset:
    #     mindist_genes = mindist_genes.union(locus_mindsistgeneset[locus])

    mindist_genes = set(summary_df['locus_mindistgene'])

    G = gv.Digraph(graph_name, filename=img_file,
                   engine='fdp', strict=True)  # 'neato' format='png'
    G.attr(bgcolor="lightgrey")
    # PPI  pairs
    selected_ppi = dict()
    for gene in selected_genes:
        if gene in network_ppi:
            paired = network_ppi[gene].intersection(selected_genes)
            if len(paired) > 0:
                selected_ppi[gene] = paired

    selected_TF = dict()
    for gene in selected_genes:
        if gene in TFsource_target:
            targets = TFsource_target[gene].intersection(selected_genes)
            selected_TF[gene] = targets

    reverse_TF = utils_ppi.reverse_network(TFsource_target)
    target_genes = set(reverse_TF.keys()).intersection(selected_genes)

    interaction_genes = set(selected_ppi.keys()).union(
        set(selected_TF.keys())).union(target_genes)

    selected_interactive_df = selected_df.loc[selected_df['gene_name'].isin(
        interaction_genes)]

    # ppi Genes - special
    G.attr('node', fixedsize='true')
    for gene in interaction_genes:
        color = 'white'
        special = gene_special[gene]
        if gene in mindist_genes:
            color = 'lightgrey'
        elif 'omim' in special:
            color = 'pink'
        elif ('exome' in special or 'lof' in special):
            color = 'orange'
        elif 'coloc' in special:
            color = 'lightyellow'

        G.node(gene, style='filled', fillcolor=color, shape='rectangle',
               width='0.6',  height='0.2', fontsize='8')

        # PPI edges
        if gene in selected_ppi:
            for pair in selected_ppi[gene]:
                if pair < gene:
                    continue
                G.edge(gene, pair, color='dimgray', arrowhead='none')

        # TF edges
        if gene in selected_TF:
            for target in selected_TF[gene]:
                G.edge(gene, target, color='darkgreen',
                       arrowhead='normal', arrowsize='0.5')

    G.render()

    return None


def neighbor_genes(innerlevel_genes, selected_geneset, network_ppi, TFsource_target, TFtarget_source):
    if SKIPPLOT:
        return None
    neighbors = set()
    for gene in innerlevel_genes:
        if gene in network_ppi:
            paired = network_ppi[gene].intersection(selected_geneset)
            neighbors = neighbors.union(paired)

        if gene in TFsource_target:
            TFtargets = TFsource_target[gene].intersection(selected_geneset)
            neighbors = neighbors.union(TFtargets)

        if gene in TFtarget_source:
            TFsources = TFtarget_source[gene].intersection(selected_geneset)
            neighbors = neighbors.union(TFsources)

    return neighbors



def plot_InteractionNetwork_v1(anchor_gene, results_dir, summary_filename,  network_ppi, gene_special, TFsource_target):
    logger.info('plotting network')

    TFtarget_source = utils_ppi.reverse_network(TFsource_target)

    summary_file = os.path.join(results_dir, summary_filename)
    summary_df = pd.read_csv(summary_file, sep='\t')
    selected_df = summary_df[summary_df['gene_signet']=='signet']
    selected_geneset = set(selected_df['gene_name'])

    mindist_geneset = set(
        summary_df.loc[summary_df['gene_mindist'] == 'mindist', 'gene_name'])

    anchor_locus = summary_df.loc[summary_df['gene_name'] == anchor_gene, 'locus_name'].values[0]
    anchor_geneset = set(summary_df.loc[summary_df['locus_name'] == anchor_locus]['gene_name'])

    anchor_interesting_geneset = {anchor_gene}
    for gene in anchor_geneset:

        if gene in gene_special and len(gene_special[gene]) > 0:
            anchor_interesting_geneset.add(gene)
        if gene in mindist_geneset:
            anchor_interesting_geneset.add(gene)

    firstlevel_neighbors = neighbor_genes(anchor_interesting_geneset, selected_geneset, network_ppi, TFsource_target, TFtarget_source)
    secondlevel_neighbors = neighbor_genes(firstlevel_neighbors, selected_geneset, network_ppi, TFsource_target, TFtarget_source)
    thirdlevel_neighbors = neighbor_genes(secondlevel_neighbors, selected_geneset, network_ppi, TFsource_target, TFtarget_source)

    # for level in [1, 2, 3]:
    for level in [1]:
        graphresults_dir = results_dir + '/network_plots'
        graph_name = os.path.join(graphresults_dir, anchor_gene + '_level' +str(level))


        if level == 1:
            neighbors = firstlevel_neighbors
        if level == 2:
            neighbors = firstlevel_neighbors.union(secondlevel_neighbors)
        if level == 3:
            neighbors = secondlevel_neighbors.union(thirdlevel_neighbors)

        # G=gv.Digraph(graph_name, filename = graph_name,
        #                engine = 'fdp', strict = True)  # 'neato', 'circo', 'twopi' format='png'
        G=gv.Digraph(graph_name, filename = graph_name, engine= 'fdp', strict = True)  # 'neato' format='png'
        G.attr(bgcolor = 'white')

        for gene in neighbors:
            color='white'
            special=gene_special[gene]
            # if gene in mindist_geneset:
            #     color = 'lightgrey'
            if 'omim' in special:
                color='pink'
            elif ('exome' in special or 'lof' in special):
                color='orange'
            elif 'coloc' in special:
                color='lightyellow'

            print(gene, special, color)

            G.node(gene, style = 'filled', fillcolor = color, shape = 'ellipse',
                   width='0.6',  height='0.2', fontsize='8')


        # anchor_genes
        G.attr('node', fixedsize = 'true')
        for gene in anchor_interesting_geneset:
            color='lightgray'
            if gene == anchor_gene:
                color='lightgreen'
            if 0: #make special anchor genes colored
                special=gene_special[gene]
                if 'omim' in special:
                    color='pink'
                elif ('exome' in special or 'lof' in special):
                    color='orange'
                elif 'coloc' in special:
                    color='lightyellow'
                if gene == anchor_gene:
                    color='lightgreen'

            G.node(gene, style = 'filled', fillcolor = color, shape = 'rectangle',
                   width = '0.6',  height = '0.2', fontsize = '8')

        network_genes = anchor_interesting_geneset.union(neighbors)

        for gene in network_genes:
            # PPI edges
            if gene in network_ppi:
                for pair in network_ppi[gene].intersection(network_genes):
                    if pair < gene:
                        continue
                    G.edge(gene, pair, color='dimgray', arrowhead='none')

            # TF edges
            if gene in TFsource_target:
                for target in TFsource_target[gene].intersection(network_genes):
                    G.edge(gene, target, color='darkgreen',
                           arrowhead='normal', arrowsize='0.5')

            if gene in TFtarget_source:
                for source in TFtarget_source[gene].intersection(network_genes):
                    G.edge(source, gene, color='darkgreen',
                           arrowhead='normal', arrowsize='0.5')

        G.render()



    return None


def plot_scorebreakdown(datafile, results_dir, base_filename):
    if SKIPPLOT:
        return None
    df = pd.read_csv(datafile, delimiter='\t')
    # m = ~df['gene_pattern'].isnull()
    # m2 = df[m]['gene_pattern'].values

    df2 = df.groupby([df.gene_pattern]).count()['locus_name']
    fig1 = plt.gcf()
    df2.plot(kind='bar')
    fig1.suptitle('categories: dist-special-ppi-TF')
    file = os.path.join(results_dir, base_filename +
                        '_bestscore_breakdown.png')
    fig1.savefig(file)
    plt.close()


def plot_nchange(pass_locus_change, results_dir, base_filename):
    if SKIPPLOT:
        return None
    # for each pass, number of loci that changed active genes
    nchange_series = []
    for p in pass_locus_change:
        nchange_series.append(len(pass_locus_change[p]))

    plt.plot(range(len(nchange_series)), nchange_series)
    plt.title('number of loci which changed as pass number')
    file = os.path.join(results_dir, base_filename + '_nchange_series.png')
    plt.savefig(file)
    plt.close()


def plot_changed_loci(pass_locus_change, results_dir, base_filename):
    if SKIPPLOT:
        return None
    # after burn-in, which loci changed genes and how often
    locus_nchange_series = dict()
    for p in range(len(pass_locus_change)):
        for locus, nchanges in pass_locus_change[p]:
            locus_nchange_series[locus] = nchanges

    loci_changed = []
    loci_nchange = []
    for locus in sorted(locus_nchange_series):
        loci_changed.append(locus)
        loci_nchange.append(locus_nchange_series[locus])

    plt.bar(loci_changed, loci_nchange)
    labels = [x[0:25] for x in loci_changed]
    plt.xticks(np.arange(len(labels)), labels,
               rotation=45, fontsize=8, ha='right')
    plt.title('nchanges of each locus, post burn-in')
    plt.tight_layout()
    file = os.path.join(results_dir, base_filename + '_loci_changed.png')
    plt.savefig(file)
    plt.close()  # clf()


def plot_convergence(maxdiff_list, results_dir, base_filename):
    if SKIPPLOT:
        return None
    # # EVALUATION: convergence
    plt.plot(range(len(maxdiff_list)), maxdiff_list)
    plt.title('network covergence, post burn-in')
    file = os.path.join(results_dir, base_filename + '_convergence.png')
    plt.savefig(file)
    plt.close()


def plot_networkscore(networkscore_list, locus_activeset, results_dir, base_filename):
    if SKIPPLOT:
        return None
    # # EVALUATION: network total score
    nloci = len(locus_activeset)
    network_avg_score = (1 / nloci) * np.array(networkscore_list)
    plt.plot(network_avg_score)
    plt.title('average score of network')
    file = os.path.join(results_dir, base_filename + '_networkscore.png')
    plt.savefig(file)
    plt.close()


def plot_scoregap(scores_logfile, results_dir, base_filename):
    if SKIPPLOT:
        return None
    # # EVALUATION: which score (distance, special, ppi, TF) lead to selecting final gene
    # # coded in that order, i.e. (1000) means selected gene had best distance score in locus
    gap_logfile = os.path.join(
        results_dir, base_filename + '_finalpass_gap.txt')

    fp = open(gap_logfile, 'w')
    fp.write('\t'.join(['locus_name', 'score_gap', 'candidate_gene',
                        'candidate_score', 'penulatimate_gene', 'penultimate_score']) + '\n')

    df = pd.read_csv(scores_logfile, delimiter='\t')
    locus_scoregap = []
    for locus in df['locus_name'].unique():
        tmp = df[df['locus_name'] == locus]
        locus_df_nrows = len(tmp) - 2
        locus_df = tmp.iloc[0:locus_df_nrows, :].copy()

        locus_df.sort_values(by=['total_score'], ascending=False, inplace=True)
        if len(locus_df) < 2:
            continue  # locus_gap = 99
        locus_gap = locus_df.iloc[0, :]['total_score'] - \
            locus_df.iloc[1, :]['total_score']
        locus_scoregap.append(locus_gap)

        locus_str = '\t'.join([locus, str(locus_gap),  locus_df.iloc[0, :]['gene_name'], str(
            locus_df.iloc[0, :]['total_score']), locus_df.iloc[1, :]['gene_name'], str(locus_df.iloc[1, :]['total_score'])])
        fp.write(locus_str + '\n')

    fp.close()

    plt.hist(locus_scoregap, bins=range(0, 30))
    plt.xticks(range(0, 30, 2))
    plt.title(
        'Total-score difference between highest- \n and second highest ranked genes of each locus')
    file = os.path.join(results_dir, base_filename + '_gap.png')
    plt.savefig(file)
    plt.close()


def plot_feature_weights(feature_score, results_dir, base_filename):
    if SKIPPLOT:
        return None
    distfac_active_list = feature_score['distfac_active']
    distfac_inactive_list = feature_score['distfac_inactive']
    plt.plot(distfac_active_list, 'r')
    plt.plot(distfac_inactive_list, 'b')
    plt.title('Convergence of Distance Scale Parameters')
    plt.xlabel('Number of passes through network')
    plt.legend(["active", "inactive"])
    plt.grid()
    file = os.path.join(results_dir, base_filename +
                        '_distance_scale_conv.pdf')
    plt.savefig(file)
    plt.clf()

    feature_list = []
    for feature in feature_score:
        if feature not in ['distfac_active', 'distfac_inactive']:
            feature_list.append(feature)
            plt.plot(feature_score[feature])

    plt.title('Convergence of Functional Score Weights')
    plt.xlabel('Number of passes through network')
    plt.legend(feature_list)  # (["omim", "exome", "coloc"])
    plt.grid()
    file = os.path.join(results_dir, base_filename +
                        '_functional_weight_conv.pdf')
    plt.savefig(file)
    plt.clf()


def plot_maxdistances(snp_gene_distance):
    if SKIPPLOT:
        return None
    max_distances = []
    for snp in snp_gene_distance:
        # k_genes = list(snp_gene_distance[snp].keys())[0:k]
        dist = []
        for gene in snp_gene_distance[snp]:
            dist.append(snp_gene_distance[snp][gene])
        max_dist = np.max(dist)
        max_distances.append(max_dist)

    plt.hist(max_distances)
    plt.savefig(results_dir + 'max_distanes.png')
    plt.close()
    return None


def plot_locus_ngenes_vs_width(locus_geneset, gene_bp_dict, gene_special, results_dir):
    if SKIPPLOT:
        return None
    locus_width_list = []
    locus_ngene_list = []
    locus_color_list = []
    for locus in locus_geneset:
        if locus.startswith('loc'):
            continue
        locus_ngene = len(locus_geneset[locus])
        # set color
        color = 'b'
        locus_genefeatures = set()
        for gene in locus_geneset[locus]:
            locus_genefeatures = locus_genefeatures.union(gene_special[gene])

        if 'coloc' in locus_genefeatures:
            color = 'g'

        if 'exome' in locus_genefeatures or 'lof' in locus_genefeatures:
            color = 'orange'

        if 'omim' in locus_genefeatures:
            color = 'r'

        locus_minbp = 1e24
        locus_maxbp = 0
        for gene in locus_geneset[locus]:
            locus_minbp = min(locus_minbp, int(gene_bp_dict[gene]['gene_tss']))
            locus_maxbp = max(locus_minbp, int(gene_bp_dict[gene]['gene_end']))
        locus_width = (locus_maxbp - locus_minbp)
        locus_width_list.append(locus_width)
        locus_ngene_list.append(locus_ngene)
        locus_color_list.append(color)

    width_ngene_color = tuple(
        zip(locus_width_list, locus_ngene_list, locus_color_list))
    width_ngene_color = sorted(width_ngene_color, key=lambda x: x[0])
    locus_width_list, locus_ngene_list, locus_color_list = list(
        zip(*width_ngene_color))

    colors = {'omim': 'red', 'exome': 'orange',
              'coloc': 'green', 'none': 'blue'}
    labels = list(colors.keys())
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[label])
               for label in labels]

    plt.tick_params(bottom=False)
    plt.xticks([])
    plt.scatter([str(x) for x in locus_width_list],
                locus_ngene_list, color=locus_color_list)
    plt.legend(handles, labels)

    plt.title('number of genes in each locus vs locus width')
    plt.xlabel('locus width')

    filename = os.path.join(results_dir, 'loci_gene_vs_width.pdf')
    logger.warning(f'skipping plotting {filename}')
    plt.savefig(filename)
    plt.clf()

    return None
