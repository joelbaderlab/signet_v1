#!/usr/bin/env python
#
# Copyright (c) 2023 Joel Bader and Zeinab Mousavi, Johns Hopkins University
#

import argparse
import yaml
import logging
import sys
import numpy as np
import random
import utils_data
import utils_ppi
import utils_plot
import pandas as pd
import os as os
import copy
import math
import pickle
from collections import defaultdict
import matplotlib.pyplot as plt
import graphviz as gv
from scipy.stats import chi2_contingency
from scipy.special import gammaln
from scipy import optimize
from scipy.stats import cauchy
from scipy.stats import iqr
#

import pdb

import warnings
warnings.simplefilter("error")

logging.basicConfig(format='%(levelname)s %(name)s.%(funcName)s: %(message)s')
logger = logging.getLogger('signet')
logger.setLevel(logging.INFO)

FEATURE_LIST = ['omim', 'exome', 'coloc' ] #Must be in order of importance for special_initialization to work

def get_args():
    parser = argparse.ArgumentParser(description='Select genes')
    parser.add_argument('-c', '--config', help='config file',required=False, type=str, default='./signet_config.yaml')

    parser.add_argument('-n', '--project_dirname', help='project name', required=False, type=str, default='TopMed') # default='TopMed' default='Autism'

    parser.add_argument('-p', '--phenotype_list', help='phenotype list', required=False, nargs='*', default=['QT', 'HR', 'QRS', 'PR', 'JT']) # default=['QT', 'HR', 'QRS', 'PR', 'JT'] default=['ASD']
    parser.add_argument('-k', '--k_closestgenes', help='number of closest genes',
                        required=False, type=int, default=1)
    parser.add_argument('-d', '--flankdistance', help='distance range from which to pick genes',required=False, type=int, default=250000)
    parser.add_argument('-s0', '--seed_start', help='seed0',required=False, type=int, default=0)
    parser.add_argument('-s1', '--seed_end', help='seed1',required=False, type=int, default=1)
    parser.add_argument('--init_rand', help='randomly initialize selected network', action='store_true')  # default is False
    parser.add_argument('--nodegcor', help='degree corrected null', action='store_true')  # default is False
    parser.add_argument(
        '--plot', help='plot network', action='store_true')  # default is False
    parser.add_argument('--nooptim', help='degree corrected null', action='store_true')  # default is False
    parser.add_argument('--notTF', help='transcription factor', action='store_false')  # default is True
    parser.add_argument('--unannotated', help='include unannotated genes', action='store_true')  # default is False
    parser.add_argument('--pseudogene', help='include pseudogene', action='store_true')  # default is False
    parser.add_argument('--no_antisense', help='exclude antisense genes', action='store_true')  # default is False
    parser.add_argument('--notPPI', help='not_run_ppi',action='store_false')


    args = parser.parse_args()
    logger.info(f'args: {args}')
    
    return(args)


class Config():
    # generate run configuration object from user-provided config file in yaml format
    def __init__(self, path):
        f_h = open(path, 'r')
        config = yaml.safe_load(f_h)
        f_h.close()

        self.__dict__ = config
        logger.info('config %s', str(self.__dict__))
        return


def write_pheno_locus_genes(filename, phenotype, pheno_locus_to_genes):
    with open(filename, 'w') as f:
        f.write('\t'.join(['phenotype', 'locus_name', 'gene_name', 'gene_locusdist_kbp']) + '\n')
        for locus, gene_distance in pheno_locus_to_genes.items():
            for (g, d) in gene_distance:
                f.write('\t'.join([phenotype, locus, g, str(d)]) + '\n')
    return None



def write_locus_genes(filename, locus_gene_distance, locus_stats, gene_bp_dict, gene_phenoset, gene_type):
    logger.info('writing %s', filename)
    with open(filename, 'w') as f:
        f.write('\t'.join([ 'locus_chr', 'locus_start', 'locus_end', 'locus_width', 'locus_name', 'locus_ngene', 'locus_pheno',  'gene_locusdist_kbp', 'gene_name', 'gene_type', 'gene_start', 'gene_end', 'gene_tss', 'gene_strand']) + '\n')
        # sort genes by distance and then by name
        loci = sorted(locus_gene_distance.keys())
        for locus in loci:
            tuples = locus_gene_distance[locus].items()
            mysort = sorted(tuples, key=lambda x: (x[1], x[0]))


            locus_chr = str(int(locus_stats[locus]['locus_chromosome'])).zfill(2)


            locus_start = str(int(locus_stats[locus]['locus_start']))
            locus_end = str(int(locus_stats[locus]['locus_end']))
            locus_width = str(locus_stats[locus]['locus_width'])
            locus_ngene = str(int(locus_stats[locus]['locus_ngene']))

            for (gene, distance) in mysort:
                g_type = gene_type[gene]
                phenostr = '+'.join(sorted(gene_phenoset[gene]))
                dist_str = str(float(distance) / 1000.0)
                gene_tss = str(int(gene_bp_dict[gene]['gene_tss']))
                gene_end = str(int(gene_bp_dict[gene]['gene_end']))
                gene_start = str(int(gene_bp_dict[gene]['gene_start']))
                gene_strand = gene_bp_dict[gene]['gene_strand']

                f.write('\t'.join([locus_chr, locus_start, locus_end, locus_width, locus, locus_ngene, phenostr, dist_str, gene, g_type, gene_start, gene_end, gene_tss, gene_strand]) + '\n')

    return None



def write_snp_genes(filename, snp_gene_distance, gene_phenoset, snp_bp_dict):
    logger.info('writing %s', filename)

    with open(filename, 'w') as f:
        f.write('\t'.join(['snp_name', 'snp_chr', 'snp_loc', 'gene_name', 'gene_snpdist_kbp', 'snp_ngenes', 'snp_mindistgene']) + '\n')
        # sort genes by distance and then by name
        snp_list = sorted(snp_gene_distance.keys())
        for snp in snp_list:
            if snp.startswith('loc'):
                continue

            snp_chr = str(int(snp_bp_dict[snp]['chromosome'])).zfill(2)
            snp_bp = str(int(snp_bp_dict[snp]['rsid_start']))

            tuples = snp_gene_distance[snp].items()
            mysort = sorted(tuples, key=lambda x: (x[1], x[0]))
            ngenes = len(mysort)

            min_dist = min(mysort)[1]
            snp_mindistgene_list = [x[0] for x in mysort if x[1]==min_dist] #incase there are more than one gene at min dist
            snp_mindistgenes = '+'.join(snp_mindistgene_list )


            for (g, d) in mysort:
                dist_str = str(float(d) / 1000.0)
                f.write('\t'.join([snp, snp_chr, snp_bp, g,  dist_str, str(ngenes) , snp_mindistgenes]) + '\n')
    return None


def write_locus_activeset(filename, locus_activeset, gene_special):
    logger.info('writing %s', filename)
    with open(filename, 'w') as f:
        f.write('\t'.join(['locus_name', 'selected_gene', 'gene_special']) + '\n')
        for locus in sorted(locus_activeset.keys()):
            for gene in sorted(locus_activeset[locus]):
                specialstr = ''
                if gene_special[gene]:
                    specialstr = '+'.join(sorted(list(gene_special[gene])))
                f.write('\t'.join([locus, gene, specialstr]) + '\n')
    return None


def get_gene_distance(snp_gene_distance):
    gene_distance = dict()
    # first find the closest distance
    for snp in snp_gene_distance:
        for gene in snp_gene_distance[snp]:
            d = snp_gene_distance[snp][gene]
            if not gene in gene_distance:
                gene_distance[gene] = d
            elif (d < gene_distance[gene]):
                gene_distance[gene] = d
    return (gene_distance)


def get_gene_signeddistance(snp_gene_signeddist):
    gene_distance = dict()
    # first find the closest distance
    for snp in snp_gene_signeddist:
        for gene in snp_gene_signeddist[snp]:
            d = snp_gene_signeddist[snp][gene]
            if not gene in gene_distance:
                gene_distance[gene] = d
            elif (abs(d) < abs(gene_distance[gene])):
                gene_distance[gene] = d
    return (gene_distance)


def get_gene_closestsnps(snp_gene_distance, gene_distance):
    gene_closestsnps = dict()
    for gene in gene_distance:
        gene_closestsnps[gene] = set()
        d0 = gene_distance[gene]
        for snp in snp_gene_distance:
            if gene in snp_gene_distance[snp]:
                if snp_gene_distance[snp][gene] == d0:
                    gene_closestsnps[gene].add(snp)

        assert len(gene_closestsnps[gene]) > 0, ('error, did not find a closest snp for gene %s' % gene)
    return(gene_closestsnps)



def write_gene_closest(filename, gene_distance, gene_closestsnps):
    logger.info('writing gene_closestto %s', filename)
    fp = open(filename, 'w')
    fp.write('\t'.join(['gene_name', 'gene_closestsnpdist_kbp', 'gene_closestsnp']) + '\n')
    for gene in gene_distance:
        d = str(float(gene_distance[gene])/1000.0)
        snpstr = ';'.join(sorted(gene_closestsnps[gene]))
        fp.write('\t'.join([gene, str(d), snpstr]) + '\n')
    fp.close()
    return None


def get_locus_stats(snp_bp_dict, locus_geneset, gene_bp_dict, results_dir):

    locus_stats = defaultdict(dict)
    df = pd.DataFrame()
    for locus in locus_geneset:
        ngene = 0
        bp_list = []
        chr_set = set()


        if not locus.startswith('loc'):
            locus_snps = locus.split('+')
            for snp in locus_snps:
                bp_list.append(int(snp_bp_dict[snp]['rsid_start']))

        for gene in locus_geneset[locus]:
            ngene += 1
            gene_tss = int(gene_bp_dict[gene]['gene_tss'])
            gene_end = int(gene_bp_dict[gene]['gene_end'])
            bp_list.append(gene_tss)
            bp_list.append(gene_end)
            locus_chromosome = gene_bp_dict[gene]['chromosome']
            chr_set.add(locus_chromosome)
        assert len(chr_set) == 1 , ("uh-oh locus spans over more than one chromosome")

        locus_start = min(bp_list)
        locus_end = max(bp_list)

        locus_width = str(float(locus_end - locus_start)/1000.0)
        locus_ngene = str(ngene)

        locus_stats[locus] = {'locus_chromosome':locus_chromosome , 'locus_start':str(locus_start), 'locus_end':str(locus_end), 'locus_width': locus_width, 'locus_ngene': locus_ngene}

    return locus_stats

def write_loci_aggregatedstats(results_dir, locus_stats):
    ngene_list = []
    width_list = []
    n_omimloci = 0


    for locus in locus_stats:
        bp_list = []
        if locus.startswith('loc'):
            n_omimloci += 1
        width_list.append(float(locus_stats[locus]['locus_width']))
        ngene_list.append(int(locus_stats[locus]['locus_ngene']))

    ngene = sum(ngene_list)

    loci_stats = {'nloci': len(width_list), 'n_omimloci': n_omimloci, 'width_min': min(width_list), 'width_med': np.median(
        width_list), 'width_max': max(width_list), 'ngene': ngene, 'ngene_min': min(ngene_list), 'ngene_med': np.median(
        ngene_list), 'ngene_max': max(ngene_list)}

    df = pd.DataFrame([loci_stats])
    filename = os.path.join(results_dir, 'loci_aggregatedstats.csv')
    df.to_csv(filename, index=True, header=True)
    return None



def initialize_activeset_random(locus_geneset, gene_distance, results_dir):
    print('initialize activeset randomly')
    locus_activeset = defaultdict()
    # otherwise there is randomness in set, no matter the seed!
    locuslist_sorted = sorted(locus_geneset.keys())
    for locus in locuslist_sorted:
        locus_genelist = sorted(locus_geneset[locus])
        recs = [(g, gene_distance[g]) for g in locus_genelist]
        locus_activeset[locus] = {random.choice(recs)[0]}

    init_activeset_df = pd.DataFrame([locus_activeset]).T
    filename = os.path.join(results_dir, 'init_activeset.csv')
    init_activeset_df.to_csv(filename)
    return(locus_activeset)

def initialize_activeset_special(locus_geneset, gene_special,  gene_distance, results_dir, feature_list):
    print('initialize activeset')

    locus_activeset = defaultdict(set)
    for locus in sorted(locus_geneset):
        generecs_list = [ ]
        for gene in locus_geneset[locus]:
            mylist = [1 if f in gene_special[gene] else 0 for f in feature_list ]
            mylist.append(-1 * gene_distance[gene])
            mylist.append(gene)
            generecs_list.append(mylist)

        generecs_list.sort()
        bestgene = generecs_list[-1][-1]
        locus_activeset[locus] = {bestgene}

    filename = os.path.join(results_dir, 'init_activeset.csv')
    logger.info(f'writing activeset to {filename}')
    init_activeset_df = pd.DataFrame([locus_activeset]).T
    init_activeset_df.to_csv(filename)
    logger.info(f'... done writing activeset')

    return(locus_activeset)




def feature_scorelist_update(feature_scorelist, locus_geneset, locus_activeset, gene_special, gene_distance,  universe_specialcounts):
    active_specialcounts = get_specialcounts(locus_activeset, gene_special )
    ngenes = len(get_networkgenes(locus_geneset))
    for feature in active_specialcounts:
        n_activegenes = len(locus_activeset) #FOR NOW ONE GENE PER LOCUS
        n_11 = active_specialcounts[feature]
        n_01 = n_activegenes - n_11 #active, but not feature
        n_10 = universe_specialcounts[feature] - n_11#not active, but has feature
        n_00 = (ngenes - len(locus_activeset)) - (n_10)
        num = (n_11+1)/(n_01+1)
        den = (n_10+1)/(n_00+1)
        score = math.log(num/den)
        feature_scorelist[feature].append(score)

    active_distances = get_active_distance(locus_activeset, gene_distance)
    active_distances = [ abs(x) for x in active_distances ]


    res_active = sum(active_distances) / len(active_distances)
    feature_scorelist['distfac_active'].append(res_active)

    return feature_scorelist

def get_distancescore_cauchy(distance, feature_scorelist):
    a = feature_scorelist['distfac_active'][-1]
    b = feature_scorelist['distfac_inactive'][-1]

    ret =  ((distance * distance) + (b * b)) / ((distance * distance) + (a * a))
    ret = np.log(ret)
    return(ret)


def get_distancescore_noinactive_old(distance, feature_scorelist):
    a = feature_scorelist['distfac_active'][-1]

    ret =  1 / ((distance * distance) + (a * a))
    ret = np.log(ret)
    return(ret)


def get_distancescore_noinactive(distance, feature_scorelist):
    a = feature_scorelist['distfac_active'][-1]
    ret = - np.log(2*a) - abs(distance) * (1./a)

    return (ret)


def get_specialscore(gene_functions, feature_scorelist):
    score = 0
    for function in gene_functions:
        score += feature_scorelist[function]
    ret = score
    return(ret)


def bayesscore_poisson(e, total, lambdazero):
    ret = gammaln(e + 1.0) - (e * math.log(lambdazero)) + lambdazero
    # leave out the constant baseline
    # ret -= math.log(total)
    return (ret)

def get_poisson_score(e, total, lambdazero):
    val = bayesscore_poisson(e, total, lambdazero)
    ret = val
    if (e < lambdazero):
        val0 = bayesscore_poisson(lambdazero, total, lambdazero)
        ret = -1* (val0 - abs(val - val0))
    return(ret)

def get_pois_fb_score(etot, e0):
    # score_pois is poisson distribution approximation
    # P(E|E0) = (E0^E / E!) exp(-E0)
    # score = ln[ P(E|E0) / P(E0/E0) ]
    # score_fb is fully bayes, integrating E0 to get P(E|alt) = const, often 1/(# of pairs plus 1)
    # and then
    # score_fb = ln[ const / P(E|E0) ] = - ln P(E|E0) + const (and we ignore the constant)

    # protect against a zero-valued etot
    e1 = etot
    if (e1 == 0):
        # score for half an edge or half of e0, whichever is smaller
        e1 = min(0.5 * e0, 0.5)
    
    score_pois = 0.0
    score_fb = 0.0
    if (e0 > 0):
        score_pois = e1 * math.log(e1 / e0) - (e1 - e0)
        score_fb = gammaln(e1 + 1.0) - (e1 * math.log(e0)) + e0
        if (e1 < e0):
            score_pois = -1.0 * score_pois
            baseline = gammaln(e0 + 1.0) - (e0 * math.log(e0)) + e0
            score_fb = baseline + abs(score_fb - baseline)
    elif (e0 == 0):
        assert(etot == 0), 'if e0 is zero, should not have edges'
    return(score_pois, score_fb)

PPISCORE_MEMORY = dict()
def get_ppi_score_within(activeset, network_ppi, degfac_ppi):
    #logger.info(f'getting ppi score ...')
    vset = activeset.intersection(set(network_ppi.keys()))
    #logger.info(f'geneset {len(activeset)}, network {len(network_ppi)}, intersection {len(vset)}')
    
    mykey = tuple(sorted(vset))
    if mykey in PPISCORE_MEMORY:
        # logger.info(f'... returning cached ppiscore, {len(PPISCORE_MEMORY)} configurations seen')
        return(PPISCORE_MEMORY[mykey])

    etot = 0
    for v1 in vset:
        etot += len(network_ppi[v1].intersection(vset))
    etot = 0.5 * etot
    e0 = 0
    degfac_list = [ degfac_ppi[v1] for v1 in vset ]
    degfacsq_list = [ x*x for x in degfac_list ]
    degfac_sum = sum(degfac_list)
    # remove the self terms
    e0 = (degfac_sum * degfac_sum) - sum(degfacsq_list)
    # account for undirected network
    e0 = 0.5 * e0
    assert(e0 >= 0), 'expect non-negative e0'

    # protect against a zero-valued etot
    e1 = etot
    if (e1 == 0):
        # score for half an edge or half of e0, whichever is smaller
        e1 = min(0.5 * e0, 0.5)

    # myscore is poisson distribution approximation
    # myscorefb is fully bayesian using 1 for prob(E|alt)
    myscore = 0.0
    myscorefb = 0.0
    if (e0 > 0):
        myscore = e1 * math.log(e1 / e0) - (e1 - e0)
        myscorefb = gammaln(e1 + 1.0) - (e1 * math.log(e0)) + e0
        if (e1 < e0):
            assert(myscore > 0.0), 'expect scores to be positive'
            myscore = -1.0 * myscore
            baseline = gammaln(e0 + 1.0) - (e0 * math.log(e0)) + e0
            myscorefb = baseline + abs(myscorefb - baseline)
    elif (e0 == 0):
        assert(etot == 0), 'if e0 is zero, should not have edges'

    (score_pois, score_fb) = get_pois_fb_score(etot, e0)
    assert(score_pois == myscore)
    assert(score_fb == myscorefb)

    ret = (myscore, myscorefb, etot, e0)
    PPISCORE_MEMORY[mykey] = ret
    return(PPISCORE_MEMORY[mykey])

TFSCORE_MEMORY = dict()
def get_tf_score_within(activeset, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev):
    #logger.info(f'getting tf score ...')
    tfset = set(network_gri.keys())
    targetset = set(network_gri_rev.keys())
    griset = tfset.union(targetset)
    vset = activeset.intersection(griset)
    
    mykey = tuple(sorted(vset))
    if mykey in TFSCORE_MEMORY:
        return(TFSCORE_MEMORY[mykey])

    # each edge has a direction
    # to count each edge once,
    # we count in the direction from each active TF to its active targets
    # have to be careful because network_gri and degfac may be undefined for vertices with 0 edges
    activetfset = tfset.intersection(vset)
    activetargetset = targetset.intersection(vset)
    etot = 0
    for v1 in activetfset:
        etot += len(network_gri[v1].intersection(vset))
    e0 = 0
    degfac_list = [ degfac_gri[v1] for v1 in activetfset ]
    degfac_rev_list = [ degfac_gri_rev[v1] for v1 in activetargetset ]
    selfset = activetfset.intersection(activetargetset)
    degfacsq_list = [ degfac_gri[v1] * degfac_gri_rev[v1] for v1 in selfset ]
    degfac_sum = sum(degfac_list)
    degfac_rev_sum = sum(degfac_rev_list)
    # remove the self terms
    e0 = (degfac_sum * degfac_rev_sum) - sum(degfacsq_list)
    assert(e0 >= 0), 'expect non-negative e0'

    (score_pois, score_fb) = get_pois_fb_score(etot, e0)
    
    ret = (score_pois, score_fb, etot, e0)
    TFSCORE_MEMORY[mykey] = ret
    # logger.info(f'*** new tf configuration number {len(TFSCORE_MEMORY)} {score_fb} {etot} {e0} ***')
    return(TFSCORE_MEMORY[mykey])


def get_locus_genescores(locus, locus_activeset, locus_geneset,  locus_score, locus_snptype,
                         gene_distance, gene_special, gene_type, feature_scorelist,
                         network_ppi, xnull, scores_logfile, gene_PPIdegree, degree_correction, ppi_Ksquared,
                         TFsource_target, reverse_TF, TF_Ksquared, gene_TFindegree, gene_TFoutdegree, ppi_run, TF_run,
                         degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev):

    # get the score for each gene
    fp = open(scores_logfile, 'a')

    # active genes at all the other loci
    other_loci = set(locus_activeset.keys()) - { locus }
    other_activegenes = set()
    for mylocus in other_loci:
        other_activegenes.update(locus_activeset[mylocus])

    #logger.info(f'{len(other_loci)} other loci, {len(other_activegenes)} other active genes')

    # pairs of universe network
    totalpair = 0.5 * len(locus_activeset) * (len(locus_activeset) - 1)

    recs = [ ]
    
    for gene in sorted(locus_geneset[locus]):

        currentset = other_activegenes.union({ gene })

        score_distance = get_distancescore_noinactive(gene_distance[gene], feature_scorelist)

        dist_str = str(float(gene_distance[gene]) / 1000.0)
        # the way to initialize gene_special is as an empty set
        score_special = 0
        for feature in gene_special[gene]:
            score_special += feature_scorelist[feature][-1]

        (ppi_score, ppi_score_fb, ppi_etot, ppi_e0) = get_ppi_score_within(currentset, network_ppi, degfac_ppi)
        score_ppi = ppi_score_fb

        edges_gene_to_others = int(0)  # for printing in output file
        other_edges = int(0)  # for printing in output file
        ppi_connected_genes_str = '.'  # for printing in output file

        # total number of ppi between other_activegenes
        # skip because takes too long to calculate
        # other_edges = utils_ppi.count_undirected_within_set(other_activegenes, network_ppi)


        ppi_connected_genes = { }
        if gene in network_ppi:
            ppi_connected_genes = network_ppi[gene].intersection(currentset)
        edges_gene_to_others = len(ppi_connected_genes)

        ppi_toks = []
        for g in sorted(ppi_connected_genes):
            specialstr = ''
            if gene_special[g]:
                specialstr = ':' + ','.join(sorted(list(gene_special[g])))
            ppi_toks.append(g + specialstr)
        ppi_connected_genes_str = ';'.join(ppi_toks)


        # Transcription Factor
        score_TF = 0.0
        TF_out_active_names = '.'  # for printing in output file
        TF_in_active_names = '.'  # for printing in output file
        TFedges_gene_str = '.'

        other_TFedges = utils_ppi.count_directed_within_set(other_activegenes, TFsource_target)
        outconnected_genes = set()
        inconnected_genes = set()
        if gene in TFsource_target:
            outconnected_genes = TFsource_target[gene].intersection(other_activegenes)
        if gene in reverse_TF:
            inconnected_genes = reverse_TF[gene].intersection(other_activegenes)
        TFedges_gene_to_others = len(outconnected_genes.union(inconnected_genes))

        Din, Dout, DinDout = utils_ppi.get_directedD_other(other_activegenes, gene_TFoutdegree, gene_TFindegree)

        TFx_gene = utils_ppi.get_directedlambda_gene(gene, gene_TFoutdegree, gene_TFindegree, Din, Dout, DinDout, TF_Ksquared)
        score_TF = get_poisson_score(other_TFedges + TFedges_gene_to_others, totalpair, TFx_gene)

        (tf_score_pois, tf_score_fb, tf_etot, tf_e0) = get_tf_score_within(currentset,
                                                                           network_gri, network_gri_rev,
                                                                           degfac_gri, degfac_gri_rev)
        score_TF = tf_score_fb

        TF_outdegree_str = '0'
        TF_out_active_cnt = '0'
        if (outconnected_genes):
            TF_outdegree_str = str(int(len(TFsource_target[gene])))
            TF_out_active_cnt = str(int(len(outconnected_genes)))

            tmp = []
            for g in sorted([x for x in outconnected_genes]):
                specialstr = ''
                if gene_special[g]:
                    specialstr = ':' +  ','.join(sorted(list(gene_special[g])))
                tmp.append(g + specialstr)

            TF_out_active_names = ';'.join(tmp)

        TF_indegree_str = '0'
        TF_in_active_cnt = '0'
        if inconnected_genes: 
            TF_indegree_str = str(int(len(reverse_TF[gene])))
            TF_in_active_cnt = str(int(len(inconnected_genes)))
            tmp = []
            for g in sorted([x for x in inconnected_genes]):
                specialstr = ''
                if gene_special[g]:
                    specialstr = ':' +  ','.join(sorted(list(gene_special[g])))
                tmp.append(g + specialstr)
            TF_in_active_names = ';'.join(tmp)

        score = score_distance + score_special + score_ppi + score_TF

        recs.append((score, gene))

        gene_special_str = ','.join([x for x in gene_special[gene]])
        if not gene_special_str:
            gene_special_str = '.'

        mystr = '\t'.join([locus, gene, gene_special_str,
                           str(score_special), dist_str, str(score_distance),
                           str(score_ppi), str(int(ppi_etot)), str(int(edges_gene_to_others)),
                           ppi_connected_genes_str, str(int(gene_PPIdegree[gene])),
                           str(score_TF),  TF_indegree_str, TF_in_active_cnt, TF_in_active_names,
                           TF_outdegree_str, TF_out_active_cnt, TF_out_active_names,
                           str(score), gene_type[gene], locus_snptype[locus], '.'])
        #logger.info('%s', mystr)
        fp.write(mystr + '\n')

    fp.close()

    return recs


def onelocus(locus, locus_activeset, locus_nchange, locus_activelist_series, locus_geneset,  locus_score, locus_mindsistgeneset, locus_snptype, gene_distance, gene_special, gene_type,  feature_scorelist,  network_ppi,  xnull, scores_logfile, gene_PPIdegree, degree_correction, ppi_Ksquared, TFsource_target, reverse_TF, TF_Ksquared, gene_TFindegree, gene_TFoutdegree, ppi_run, TF_run,
             degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev):

    recs = get_locus_genescores(locus, locus_activeset, locus_geneset,  locus_score, locus_snptype, gene_distance, gene_special, gene_type, feature_scorelist,  network_ppi,  xnull, scores_logfile, gene_PPIdegree, degree_correction, ppi_Ksquared, TFsource_target, reverse_TF, TF_Ksquared, gene_TFindegree, gene_TFoutdegree, ppi_run, TF_run,
                                degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev)

    recs = sorted(recs)
    (bestscore, bestgene) = recs[-1]

    bestgenes = set()
    for score, gene in recs:
        ###the AND condition is to avoid selecting a gene that is in neither category.
        if score == bestscore:
            bestgenes.add(gene)


    candidate_bestgene = random.choice(sorted(bestgenes)) # sorted so randomness is reproducable!


    write_scorebreakdown(locus, candidate_bestgene, scores_logfile, network_ppi, TFsource_target, reverse_TF)
    new_geneset = {candidate_bestgene}

    locus_nchange[locus]['changed'] = False
    if locus_activeset[locus] != new_geneset:
        locus_nchange[locus]['nchange'] += 1
        locus_nchange[locus]['changed'] = True

    locus_activeset[locus] = new_geneset
    locus_score[locus].append(bestscore)
    locus_activelist_series[locus].append(list(new_geneset)[0])  # ONLY 1


    return None


def write_scorebreakdown(locus, candidate_bestgene, scores_logfile, network_ppi, TFsource_target, reverse_TF):
    # # EVALUATION: which score (distance, special, ppi, TF) lead to selecting final gene
    # # coded in that order, i.e. (1000) means selected gene had best distance score in locus
    # changed this so that we have a "1" for special or ppi or TF as long as it is valid for the gene, even if its not the highest in that locus
    locus_scorebreakdown = []
    df = pd.read_csv(scores_logfile, delimiter='\t')



    locus_df = df[df['locus_name'] == locus].copy()



    max_distance_score = np.max(locus_df['distance_score'])
    max_special_score = np.max(locus_df['special_score'])
    max_ppi_score = np.max(locus_df['ppi_score'])
    max_TF_score = np.max(locus_df['TF_score'])
    max_total_score = np.max(locus_df['total_score'])


    candidate_bestgene_distance_score = locus_df[locus_df['gene_name'] == candidate_bestgene]['distance_score'].values[0]
    candidate_bestgene_special_score = locus_df[locus_df['gene_name']== candidate_bestgene]['special_score'].values[0]
    candidate_bestgene_ppi = locus_df[locus_df['gene_name']== candidate_bestgene]['ppi_deg_active_cnt'].values[0]
    candidate_bestgene_TF_score = locus_df[locus_df['gene_name'] == candidate_bestgene]['TF_score'].values[0]

    candidate_bestgene_TF = locus_df[locus_df['gene_name']== candidate_bestgene]['TF_in_active_cnt'].values[0] + locus_df[locus_df['gene_name']== candidate_bestgene]['TF_out_active_cnt'].values[0]

    special_score_locus ='x'
    if (len(np.unique(locus_df['special_score'])) > 1) and (candidate_bestgene_special_score == 0):
        special_score_locus = '0'
    if (candidate_bestgene_special_score > 0):
        special_score_locus = '1'
    if (len(np.unique(locus_df['special_score'])) == 1) and (candidate_bestgene_special_score == 0):
        special_score_locus = '_'

    ppi_score_locus = '_'
    if candidate_bestgene_ppi > 0:
        ppi_score_locus = 1

    TF_score_locus = '_'
    if candidate_bestgene_TF != '.':
        TF_score_locus = 1



    score_breakdown = ['s', int(candidate_bestgene_distance_score == max_distance_score), special_score_locus, ppi_score_locus, TF_score_locus]

    score_breakdown_str = ''.join([str(x) for x in score_breakdown])



    df.loc[(df['locus_name'] == locus) & (df['gene_name'] == candidate_bestgene) , 'gene_pattern'] = score_breakdown_str
    df.to_csv(scores_logfile, sep='\t', index=False)


def get_active_distance(locus_activeset, gene_distance):
        dist = []
        for locus in locus_activeset:
            if locus.startswith('loc'):
                continue
            for gene in locus_activeset[locus]: #only one gene but this is a set
                break
            dist.append(gene_distance[gene])
        distances = np.array(dist)
        return distances

def get_inactive_distance(locus_geneset, locus_activeset, gene_distance):
        dist = []
        for locus in locus_geneset: #assume two final genes with same scores have same distance
            if locus.startswith('loc'):
                continue
            for active_gene in locus_activeset[locus]: #only one gene but this is a set
                break
            for gene in locus_geneset[locus]:
                if gene.startswith('locus_name'):
                    break
                if gene == active_gene:
                    continue
                dist.append(gene_distance[gene])

        distances = np.array(dist)
        return distances

# def distance_func(x, distances):
#     n = len(distances)
#     return (1/(x**2)) - (2/n)*np.sum(1/(x**2 + distances**2))  # only one real root at x = 1

def onepass(locus_activeset, locus_activelist_series, locus_geneset, locus_score,
            locus_nchange, pass_locus_change, locus_mindsistgeneset, gene_distance,
            gene_special, gene_type, locus_snptype, feature_scorelist,
            network_ppi,  xnull, scores_logfile, gene_PPIdegree, degree_correction, ppi_Ksquared,
            TFsource_target, reverse_TF, TF_Ksquared, gene_TFindegree, gene_TFoutdegree, ppi_run, TF_run,
            degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev):

    locuslist_sorted = sorted(locus_geneset.keys())
    locuslist_randomorder = random.sample(
        locuslist_sorted, len(locuslist_sorted))

    fp = open(scores_logfile, 'w')
    fp.write('\t'.join(['locus_name', 'gene_name', 'gene_special',
                        'special_score',  'gene_locusdist_kbp', 'distance_score', 'ppi_score', 'edges_in_otheractive', 'ppi_deg_active_cnt', 'ppi_deg_active_names', 'ppi_deg', 'TF_score', 'TF_indeg', 'TF_in_active_cnt', 'TF_in_active_names', 'TF_outdeg', 'TF_out_active_cnt', 'TF_out_active_names','total_score', 'gene_type', 'snp_type', 'gene_pattern']) + '\n')
    fp.close()

    for locus in locuslist_randomorder:

        onelocus(locus, locus_activeset, locus_nchange, locus_activelist_series, locus_geneset,  locus_score,
                 locus_mindsistgeneset, locus_snptype, gene_distance, gene_special, gene_type,  feature_scorelist,
                 network_ppi,  xnull, scores_logfile, gene_PPIdegree, degree_correction, ppi_Ksquared,
                 TFsource_target, reverse_TF, TF_Ksquared, gene_TFindegree, gene_TFoutdegree, ppi_run, TF_run,
                 degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev)

    npass = len(pass_locus_change)
    changed_loci = []
    for locus in locus_nchange:
        if locus_nchange[locus]['changed']:
            changed_loci.append((locus, locus_nchange[locus]['nchange']))

    pass_locus_change[npass] = changed_loci

    return None


def get_networkgenes(network):
    ret = set()
    for locus in network:
        ret.update(network[locus])
    return (ret)

def get_specialcounts(network_dict, gene_special ):
    ##count how many genes in each special category:
    special_counts = dict()
    for f in FEATURE_LIST:
        special_counts[f] = 0
    for locus in network_dict:
        if locus.startswith('loc'):
            continue
        for gene in network_dict[locus]:
            for special in gene_special[gene]:
                assert special in special_counts, 'unknown special feature: ' + special
                special_counts[special] += 1
    return special_counts


def write_extended_finaldf(scores_logfile, seed, degree_correction, optim_run, results_dir, base_filename):
    df = pd.read_csv(scores_logfile, delimiter='\t')
    final_df = pd.DataFrame(columns = df.columns)

    final_df['exp_score'] = 0

    for locus in df['locus_name'].unique():

        locus_df = df[df['locus_name'] == locus].copy()
        locus_df['gene_mindist']= '.'
        locus_df['gene_nonzero']= '.'


        max_total_score = np.max(locus_df['total_score'])
        mindist = min(locus_df['gene_locusdist_kbp'])
        mindist_genelist = locus_df[locus_df['gene_locusdist_kbp']==mindist]['gene_name'].values.tolist()
        mindist_genestr = '+'.join(sorted(mindist_genelist))

        maxppi = max(locus_df['ppi_score'])
        ppi_diff = len(np.unique(locus_df['ppi_score']))



        scores = np.array(locus_df['total_score'])


        exp_scores = np.exp(scores-max_total_score)
        sum_scores = np.sum(exp_scores)
        exp_scores = exp_scores/sum_scores
        locus_df['exp_score'] = exp_scores

        for idx, row in locus_df.iterrows():


            if row['gene_name'] in mindist_genelist:
                locus_df.at[idx, 'gene_mindist'] = 'mindist'

            if not row['gene_pattern'] == '.':
                locus_df.at[idx, 'gene_nonzero'] = 'nonzero'


            if final_df.empty:
                final_df = locus_df.loc[[idx]]
            else:
                final_df = pd.concat([final_df, locus_df.loc[[idx]]])


    # write the new df:
    file = os.path.join(results_dir, base_filename + '_extended_finalpass_scores.txt')
    logger.info(f'writing {file}')
    final_df.to_csv(file, sep='\t')


    return  final_df


def write_feature_log(feature_scorelist, results_dir, base_filename):
    for feature in feature_scorelist:
        feature_logfile = os.path.join(results_dir, base_filename + '_'  + feature + '_log.txt')
        df = pd.DataFrame(feature_scorelist[feature], columns=[feature])
        df.to_csv(feature_logfile, sep='\t', index=False)


def write_networkscore_log(networkscore_list, seed, results_dir, base_filename):
    logfile = os.path.join(results_dir, base_filename + '_networkscore_log.txt')
    df = pd.DataFrame(networkscore_list, columns=['network_score'])
    df['seed'] = seed
    df.to_csv(logfile, sep='\t', index=False)


def run_experiment(phenotype_list, k, FLANKDIST, degree_correction, optim_run, network_plot,
                   gene_PPIdegree, ppi_Ksquared, network_ppi, xnull, ppi_run,
                   TF_run, TFsource_target, reverse_TF, TF_Ksquared,
                   gene_TFindegree, gene_TFoutdegree, locus_geneset, networkgenes,
                   gene_special, locus_mindsistgeneset, gene_distance, gene_signeddistance,
                   gene_type, locus_snptype, external_genes, seed, results_dir, seedless_base_filename,
                   counts_log, unannotated, pseudogene, antisense,
                   degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev, init_rand):
    random.seed(seed)
    logger.info(f'*** starting with seed {seed} ***')

    #base_filename = 'seed' + str(seed) + '_ppi' + str(ppi_run) + '_TF' + str(TF_run) + '_degcor' + str(degree_correction) + '_optim' + str(optim_run) +'_unannotated' + str(unannotated) + '_pseudogene' + str(pseudogene) + '_no-antisense' + str(antisense)
    base_filename = f'seed{seed:02d}'

    scores_logfile = os.path.join(results_dir, base_filename + '_finalpass_scores.txt')

    pseleted_logfile = os.path.join(results_dir, base_filename + '_pselected.txt')

    initscores_logfile = os.path.join(results_dir, seedless_base_filename + '_initialconfig_scores.txt')
    fp = open(initscores_logfile, 'w')
    fp.write('\t'.join(['locus_name', 'gene_name', 'gene_special',
                        'special_score',  'gene_locusdist_kbp', 'distance_score', 'ppi_score',
                        'edges_in_otheractive', 'ppi_deg_active_cnt', 'ppi_deg_active_names', 'ppi_deg',
                        'TF_score', 'TF_indeg', 'TF_in_active_cnt', 'TF_in_active_names',
                        'TF_outdeg', 'TF_out_active_cnt', 'TF_out_active_names',
                        'total_score', 'gene_type', 'snp_type', 'gene_pattern']) + '\n')
    fp.close()
    
    ##########################
    # Network Optimization
    ##########################


    logger.info('initializing ...')

    universe_specialcounts = get_specialcounts(locus_geneset, gene_special )
    logger.info(f'universe specialcounts {universe_specialcounts}')

    if init_rand:
        locus_activeset = initialize_activeset_random(locus_geneset, gene_distance, results_dir)
    if not init_rand:
        locus_activeset = initialize_activeset_special(locus_geneset, gene_special,  gene_distance, results_dir, FEATURE_LIST)

    init_activeset = copy.deepcopy(locus_activeset)


    logger.info('... done initializing active set')

    pass_locus_change = dict()
    locus_score = defaultdict(list)
    networkscore_list = []

    locus_nchange = defaultdict(dict)
    for locus in locus_activeset:
        locus_nchange[locus]['nchange'] = 0
        locus_nchange[locus]['changed'] = False

    locus_activelist_series = defaultdict(list)
    for locus in locus_activeset:
        locus_activelist_series[locus].append(list(locus_activeset[locus])[0])

    gene_nselected = dict()  # number of times selected
    gene_pselected = dict()  # probability selected, equal to ntimes / npass
    # initialize

    for locus in locus_geneset:
        for gene in locus_geneset[locus]:
            gene_nselected[gene] = 0  # integer number of times selected
            # number of times divided by number of passes but as floats!!!!!
            gene_pselected[gene] = 0.0

    npass = 10
    npass_completed = 0
    CONVERGENCE = 0.01

    maxdiff_list = []
    maxdiff_gene_list = []
    converged = False

    feature_scorelist = defaultdict(list)
    feature_scorelist = feature_scorelist_update(feature_scorelist, locus_geneset, locus_activeset, gene_special, gene_distance, universe_specialcounts)

    # calculate network of the initial configuration
    networkscore = 0.0
    for locus in locus_activeset:
        for gene in locus_activeset[locus]:

            logger.info(f'... starting gene {gene} and locus {locus} getting score (of initalized activeset)')

            # usually the first arguments to get_locus_scores are locus, locus_activeset, locus_geneset
            # here the locus_geneset is replaced by locus_activeset
            #    to restrict the calculation to just the active gene at each locus
            initrecs = get_locus_genescores(locus, locus_activeset, locus_activeset,  locus_score, locus_snptype,
                                            gene_distance, gene_special, gene_type, feature_scorelist,  network_ppi,  xnull,
                                            initscores_logfile, gene_PPIdegree, degree_correction, ppi_Ksquared,
                                            TFsource_target, reverse_TF, TF_Ksquared, gene_TFindegree, gene_TFoutdegree,
                                            ppi_run, TF_run, degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev)

            # initrecs is a list of (score, gene) tuples
            # initrecs[0][0] is the scalar value of the score of the first gene
            # since there is only one active gene, this is the score of the active gene
            initscore = initrecs[0][0]
            locus_score[locus].append(initscore)
            # since the initscore was just appended, this score is added to the networkscore
            # this is probably not right because it adds the ppi and gri scores repeatedly
            networkscore += locus_score[locus][-1]

    networkscore_list.append(networkscore)


    for pass_number in range(npass):

        logger.info(f'... starting seed {seed} pass {pass_number + 1}')
        
        onepass(locus_activeset, locus_activelist_series, locus_geneset, locus_score, locus_nchange, pass_locus_change,
                locus_mindsistgeneset, gene_distance, gene_special, gene_type, locus_snptype, feature_scorelist,
                network_ppi,  xnull, scores_logfile, gene_PPIdegree, degree_correction, ppi_Ksquared,
                TFsource_target, reverse_TF, TF_Ksquared, gene_TFindegree, gene_TFoutdegree,
                ppi_run, TF_run, degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev)

        feature_scorelist = feature_scorelist_update(feature_scorelist, locus_geneset, locus_activeset,
                                                     gene_special, gene_distance, universe_specialcounts)

        # if (network_plot and (not pass_number % 10)):
        #     utils_plot.plot_network(len(pass_locus_change), results_dir, locus_activeset,
        #                  network_ppi, phenotype_list, gene_phenoset, gene_special, degree_correction, TFsource_target)

        assert len(locus_score) == len(locus_geneset), 'uh-oh locus_score has more/less loci than locus_geneset'

        # calculate network score of this pass
        networkscore_list_old = []
        networkscore = 0.0
        for locus in locus_score:
            networkscore += locus_score[locus][-1]
        networkscore_list_old.append(networkscore)

        # calculate network of the current configuration
        final_locus_score = defaultdict()
        networkscore = 0.0
        for locus in locus_activeset:
            for gene in locus_activeset[locus]:

                finalrecs = get_locus_genescores(locus, locus_activeset, locus_activeset,  final_locus_score,
                                                 locus_snptype, gene_distance, gene_special, gene_type, feature_scorelist,
                                                 network_ppi,  xnull, initscores_logfile, gene_PPIdegree, degree_correction, ppi_Ksquared,
                                                 TFsource_target, reverse_TF, TF_Ksquared, gene_TFindegree, gene_TFoutdegree,
                                                 ppi_run, TF_run,
                                                 degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev)
                finalscore = finalrecs[0][0]
                networkscore += finalscore

        networkscore_list.append(networkscore)


        for locus in locus_activeset:
            for g in locus_activeset[locus]:
                gene_nselected[g] = gene_nselected[g] + 1


        # maximum difference in pselected
        maxdiff = 0.0
        for g in gene_nselected:
            oldval = gene_pselected[g]
            newval = float(gene_nselected[g]) / float(pass_number + 1)
            gene_pselected[g] = newval
            diff = abs(oldval - newval)
            if (diff > maxdiff):
                maxdiff = diff
        maxdiff_list.append(maxdiff)

        #parameters' maximum relative-difference
        max_paramdiff = 0.0
        for feature in feature_scorelist:
            oldval = float(feature_scorelist[feature][-2])
            newval = float(feature_scorelist[feature][-1])
            if oldval == 0:
                oldval =0.1 * CONVERGENCE
            diff = abs(newval - oldval) / oldval
            if (diff > max_paramdiff):
                max_paramdiff = diff

        current_pass = pass_number
        nchanges = len(pass_locus_change[current_pass])


        if ((maxdiff <= CONVERGENCE and max_paramdiff <= CONVERGENCE)  or (nchanges == 0)  ):
            logger.info(f'... seed {seed} finished after pass {pass_number + 1}')
            converged = True

            break

    if not converged:
        logger.info('uh-oh! did not converge after %d passes', pass_number)

    final_datafile =  os.path.join(results_dir, base_filename + '_finalpass_scores.txt')
    #utils_plot.plot_scorebreakdown(final_datafile, results_dir, base_filename)

    # # EVALUATION plots ###########

    # utils_plot.plot_nchange(pass_locus_change, results_dir, base_filename)
    # utils_plot.plot_changed_loci(pass_locus_change, results_dir, base_filename)
    # utils_plot.plot_convergence(maxdiff_list, results_dir, base_filename)
    # utils_plot.plot_networkscore(networkscore_list, locus_activeset, results_dir, base_filename)
    # utils_plot.plot_scoregap(scores_logfile, results_dir, base_filename)
    # utils_plot.plot_feature_weights(feature_scorelist, results_dir, base_filename)
    write_feature_log(feature_scorelist, results_dir, base_filename)
    write_networkscore_log(networkscore_list, seed, results_dir, base_filename)


    # # EVALUATION ###########
    list_gene_counts = [] #[networkscore]
    for locus in sorted(locus_geneset.keys()):
        for gene in sorted(locus_geneset[locus]):
            gene_count = 0
            if gene in locus_activeset[locus]:
                gene_count = 1
            list_gene_counts.append(gene_count)

    df = pd.read_csv(counts_log, sep='\t')
    if type(df. iloc[:, 0].values[0]) !=str :
        if np.isnan(df. iloc[:, 0].values[0]):
            df = df.iloc[1:, :]
    df[str(seed)] = list_gene_counts
    df.to_csv(counts_log, index=False, sep='\t')

    logger.warning('beware, write_extended_finaldf may cause an error')
    extendedfinal_df = write_extended_finaldf(scores_logfile, seed, degree_correction, optim_run, results_dir, base_filename)
    return None

def describe_network(vertex_neighbors):
    logger.info(f'{len(vertex_neighbors)} vertices')
    my_degree = dict()
    tot_edge = 0
    for v in vertex_neighbors:
        my_degree[v] = len(vertex_neighbors[v])
        tot_edge += my_degree[v]
    logger.info(f'{tot_edge} edges, or {tot_edge/2} if undirected')
    deghist = defaultdict(int)
    for (v, deg) in my_degree.items():
        deghist[deg] += 1
    logger.info(f'average degree {tot_edge/len(my_degree)}')
    logger.info(f'degree\tnvert')
    for d in sorted(deghist.keys()):
        logger.info(f'{d}\t{deghist[d]}')
    return None

def get_uniprot_subset(uniprot_to_hgnc, networkgenes):
    hgnc_to_uniprot = dict()
    for (u, mylist) in uniprot_to_hgnc.items():
        for h in set(mylist):
            if h not in hgnc_to_uniprot:
                hgnc_to_uniprot[h] = set()
            hgnc_to_uniprot[h].add(u)
    uniprot_subset = set()
    for h in networkgenes:
        if h in hgnc_to_uniprot:
            myset = hgnc_to_uniprot[h]
            if (len(myset) > 1):
                logger.warning(f'{h} maps to multiple uniprot {myset}')
                for u in myset:
                    logger.warning(f'... {u} -> {uniprot_to_hgnc[u]}')
            uniprot_subset = uniprot_subset.union(myset)
        else:
            logger.warning(f'{h} not in uniprot')
    logger.info(f'{len(networkgenes)} genes -> {len(uniprot_subset)} uniprot')
    return uniprot_subset

def get_network_ppi(uniprot_to_hgnc, uniprot_subset, hgnc_subset, ppi_file):
    logger.info(f'reading subset from {ppi_file}')
    uniprot_directed_dict = utils_ppi.read_ppi_subset(ppi_file, uniprot_subset)
    logger.info(f'converting uniprot to hgnc')
    protein_directed_dict = utils_ppi.convert_uniprot_to_hgnc(uniprot_directed_dict, uniprot_to_hgnc)
    logger.info(f'converting network to undirected')
    protein_proteins = utils_ppi.make_undirected(protein_directed_dict)

    network_ppi = dict()
    for p in protein_proteins:
        if p in hgnc_subset:
            network_ppi[p] = protein_proteins[p].intersection(hgnc_subset)
            network_ppi[p].discard(p)
            if len(network_ppi[p]) == 0:
                del network_ppi[p]
    return(network_ppi)

def write_networkfile(filename, network_dict):
    logger.info(f'writing network to {filename}')
    with open(filename, 'w') as fp:
        fp.write('\t'.join(['SIGNET_SOURCE','SIGNET_DEST']) + '\n')
        for v1 in sorted(network_dict):
            myset = network_dict[v1]
            v2 = ' '.join(sorted(list(myset)))
            fp.write('\t'.join([v1,v2]) + '\n')
    return None

def read_networkfile(filename):
    logger.info(f'reading network from {filename}')
    ret = dict()
    with open(filename, 'r') as fp:
        headertoks = fp.readline().strip().split('\t')
        assert(len(headertoks) == 2)
        assert(headertoks[0] == 'SIGNET_SOURCE')
        assert(headertoks[1] == 'SIGNET_DEST')
        for line in fp:
            toks = line.strip().split('\t')
            assert(len(toks) == 2)
            (v1, v2str) = (toks[0], toks[1])
            v2set = set(v2str.split())
            assert v1 not in ret, 'repeated vertex: ' + v1
            ret[v1] = v2set
    return(ret)

def get_degfac_ppi(network_ppi):
    # probability of an edge in an undirected network is approximately d1 d2 / 2E
    # d1 and d2 are vertex degrees
    # E is the number of undirected edges
    # degfac for a vertex is d/sqrt(2E)
    degfac_ppi = dict()
    totedge = 0.0
    for (v1, myset) in network_ppi.items():
        mydegree = float(len(myset))
        degfac_ppi[v1] = mydegree
        totedge += mydegree
    totedge = 0.5 * totedge
    logger.info(f'totedge {totedge} totvert {len(network_ppi)}')
    prefac = 1.0 / math.sqrt(2.0 * totedge)
    for (v1, myval) in degfac_ppi.items():
        degfac_ppi[v1] = prefac * myval
    return(degfac_ppi)

def get_degfac_gri(network_gri, network_gri_rev):
    # probability of an edge in a directed network is approximately kin kout / E
    # kin and kout are in and out degrees
    # E is total number of directed edges
    # degfac is kout / sqrt(E)
    # degfac_rev is kin / sqrt(E)
    degfac_gri = dict()
    degfac_gri_rev = dict()
    totedge = 0
    totedge_rev = 0
    for (v1, myset) in network_gri.items():
        outdeg = len(myset)
        degfac_gri[v1] = outdeg
        totedge += outdeg
    prefac = 1.0/math.sqrt(float(totedge))
    for (v2, myset) in network_gri_rev.items():
        indeg = len(myset)
        degfac_gri_rev[v2] = indeg
        totedge_rev += indeg
    assert(totedge == totedge_rev), f'total edge count error {totedge} {totedge_rev}'

    for (v1, mydeg) in degfac_gri.items():
        degfac_gri[v1] = prefac * float(mydeg)
    for (v2, mydeg) in degfac_gri_rev.items():
        degfac_gri_rev[v2] = prefac * float(mydeg)

    logger.info(f'vout {len(degfac_gri)} vin {len(degfac_gri_rev)} nedge {totedge}')

    return(degfac_gri, degfac_gri_rev)
        
        

GENE_BP_DICT = None
def main():

    args = get_args()

    project_dirname = args.project_dirname
    phenotype_list = args.phenotype_list
    k = args.k_closestgenes
    FLANKDIST = args.flankdistance
    degree_correction = not args.nodegcor
    optim_run = not args.nooptim
    seed_start = args.seed_start
    seed_end = args.seed_end
    TF_run = args.notTF
    ppi_run = args.notPPI
    unannotated = args.unannotated
    pseudogene = args.pseudogene
    antisense = not args.no_antisense
    network_plot = args.plot
    init_rand = args.init_rand


    config = Config(args.config)

    path = config.path
    data_dir = config.data_dir
    results_dir = config.results_dir
    gene_filename  = config.gene_filename
    ppi_filename = config.ppi_filename
    protein_gene_filename  = config.protein_gene_filename
    trrust_filename = config.trrust_filename
    gene_file = os.path.join(data_dir, gene_filename)

    protein_gene_file = os.path.join(data_dir, protein_gene_filename)
    ppi_file = os.path.join(data_dir, ppi_filename)
    trrust_file = os.path.join(data_dir, trrust_filename)

    snploc_filename  = config.snploc_filename
    gwas_filename = config.gwas_filename

    phenostr = '+'.join(phenotype_list)
    flankstr = str(int( int(FLANKDIST) / 1000)) + 'kb'
    




    if init_rand == True:
        initstr = 'initrand' + str(init_rand)
        results_dir = os.path.join(results_dir, '_'.join([project_dirname, phenostr, flankstr, initstr, 'NewDist']))
    else:
        results_dir = os.path.join(results_dir, '_'.join([project_dirname, phenostr, flankstr, 'NewDist']))


    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    logger.info(f'results directory{results_dir}')

    pdb.set_trace()

    project_datadir = os.path.join(data_dir, project_dirname)
    snploc_file = os.path.join(project_datadir, snploc_filename)
    gwas_file = os.path.join(project_datadir, gwas_filename)

    logger.info(f'datadir {project_datadir}')
    logger.info(f'snploc {snploc_file}')
    logger.info(f'gwas {gwas_file}')
    
    maxdistances_plot = False

    # external genes:
    feature_geneset = utils_data.get_feature_geneset(project_datadir, phenotype_list)
    for myfeature in feature_geneset:
        logger.info(f'feature {myfeature} -> {len(feature_geneset[myfeature])} genes')
    external_genes = set().union(*list(feature_geneset.values()))
    logger.info(f'total number of genes with features: {len(external_genes)}')

    logger.info(f'reading genes from {gene_file} ...')
    GENE_BP_DICT = utils_data.read_gene_file(gene_file, unannotated, pseudogene, antisense, project_datadir, external_genes)
    with open(results_dir + '/gene_bp_dict.pickle', 'wb') as handle:
        pickle.dump(GENE_BP_DICT, handle, protocol=pickle.HIGHEST_PROTOCOL)

    ensembl_genes = set(GENE_BP_DICT.keys())
    logger.info(f'... {len(ensembl_genes)} genes')

    missing_genes = external_genes - ensembl_genes
    nmissing = len(missing_genes)
    logger.info(f'number external genes missing: {nmissing}')
    if missing_genes:
        print ('some genes missing, check for their other names', missing_genes)

    gene_type = defaultdict()
    for g in GENE_BP_DICT:
        gene_type[g] = GENE_BP_DICT[g]['gene_type']

    pheno_snp_geneset = defaultdict(dict)
    snp_gene_distance = defaultdict(dict)
    snp_gene_signeddist = defaultdict(dict)
    pheno_geneset = defaultdict(set)
    gene_phenoset = defaultdict(set)

    rep_snp = 0

    for ph in phenotype_list:

        pheno_snp_geneset[ph] = defaultdict(dict)
        pheno_snp_gene_distance = utils_data.get_network(k, FLANKDIST, GENE_BP_DICT,
                                                         ph, project_datadir, snploc_filename)

        for snp in pheno_snp_gene_distance:

            if snp in snp_gene_distance:
                rep_snp +=1

            # will rewrite repeat snps if snp is in >1 phen but ok
            snp_gene_distance[snp] = dict()
            genes, distances = list(zip(*pheno_snp_gene_distance[snp]))
            pheno_geneset[ph] = pheno_geneset[ph].union(set(genes))
            pheno_snp_geneset[ph][snp] = set(genes)

            for gene, distance in zip(genes, distances):
                snp_gene_distance[snp][gene] = abs(distance)
                snp_gene_signeddist[snp][gene] = distance
                gene_phenoset[gene] = gene_phenoset[gene].union([ph])


    logger.info(f'GWAS gene count {len(gene_phenoset)}')
                
    # for network genes, get distance to closest snp and the closest snp
    gene_distance = get_gene_distance(snp_gene_distance)
    gene_signeddistance = get_gene_signeddistance(snp_gene_signeddist)

    # pdb.set_trace()


    for gene in gene_distance:
        assert gene_distance[gene] == abs(gene_signeddistance[gene])

    gene_closestsnps = get_gene_closestsnps(snp_gene_distance, gene_distance)

    # if maxdistances_plot:
    #     utils_ppi.plot_maxdistances(snp_gene_distance)

    write_gene_closest(os.path.join(
        results_dir, 'gene_closest.txt'), gene_distance, gene_closestsnps)

    gene_special = utils_data.get_gene_featureset(gene_distance, feature_geneset)
    with open(results_dir + '/gene_special.pickle', 'wb') as handle:
        pickle.dump(gene_special, handle, protocol=pickle.HIGHEST_PROTOCOL)

    pheno_locus_gene_distance = defaultdict(dict)
    pheno_locus_geneset = defaultdict(dict)
    for ph in phenotype_list:
        pheno_snp_gene_distance = utils_data.get_network(k, FLANKDIST, GENE_BP_DICT,
                                                         ph, project_datadir, snploc_filename)


        locus_genes_distances = utils_data.merge_pheno_snps(pheno_snp_gene_distance)

        pheno_locus_gene_distance[ph] = defaultdict(dict)
        for locus in locus_genes_distances:
            genes, distances = list(zip(*locus_genes_distances[locus]))
            pheno_locus_geneset[ph][locus] = set(genes)
            pheno_locus_gene_distance[ph][locus] = defaultdict(dict)
            for gene, distance in zip(genes, distances):
                pheno_locus_gene_distance[ph][locus][gene] = distance

        filename = os.path.join(results_dir, ph + '_locus_gene_distance.txt')
        logger.info(f'... writing {filename}')
        write_pheno_locus_genes(filename, ph,locus_genes_distances)

    # keep track of gene-phenotype connections
    # concatenate the loci of all phenotypes
    union_locus_geneset = defaultdict(set)
    union_locus_gene_distance = defaultdict(dict)
    for ph in pheno_locus_geneset:
        union_locus_geneset.update(pheno_locus_geneset[ph])
        union_locus_gene_distance.update(pheno_locus_gene_distance[ph])
        logger.warning(f'{ph} check whether this overwrites previous distances')


    # check that all special genes are in network
    # can miss if they are located > FLANKDIST from all SNPs and are not in the first k-genes of any SNP
    nw_missedgenes = set(gene_special.keys()) - set(gene_distance.keys())
    nw_missedgenes_special = [gene_special[x] for x in nw_missedgenes]

    logger.info(f'{len(nw_missedgenes)} missed genes: {nw_missedgenes}')
    logger.info(f'missed special: {nw_missedgenes_special}')


    logger.info(f'before adding missing genes, {len(gene_distance)} genes and {len(union_locus_geneset)} loci')
    for gene in nw_missedgenes:

        logger.info(f'adding {gene} to network')
        if gene not in missing_genes:
            locus_name = 'locus_' + gene
            #locus_loc = int((int(GENE_BP_DICT[gene]['gene_tss'])+
            #                 int(GENE_BP_DICT[gene]['gene_end']))/2 - int(GENE_BP_DICT[gene]['gene_tss']))
            locus_loc = 0
            gene_distance[gene] = locus_loc
            gene_signeddistance[gene] = locus_loc
            gene_closestsnps[gene] = set([locus_name])
            snp_gene_distance[locus_name][gene] = locus_loc
            snp_gene_signeddist[locus_name][gene] = locus_loc
            union_locus_geneset[locus_name] = {gene}
            union_locus_gene_distance[locus_name][gene] = locus_loc
            gene_type[gene] = list(gene_special[gene])[0]

    logger.info(f'after adding missing genes, {len(gene_distance)} genes and {len(union_locus_geneset)} loci')


    locus_to_genes = utils_data.merge_loci(union_locus_geneset, union_locus_gene_distance)

    locus_gene_distance = dict()
    locus_geneset = dict()
    for locus in locus_to_genes:
        locus_gene_distance[locus] = dict()
        locus_geneset[locus] = set()
        for (g, d) in locus_to_genes[locus]:
            locus_gene_distance[locus][g] = d
            locus_geneset[locus].add(g)


    logger.info(f'after merging, {len(locus_geneset)} loci')


    locus_snptype = utils_data.locus_snptype(locus_geneset, project_datadir, gwas_filename, phenotype_list, delimiter='\t')

    locus_mindsistgeneset = defaultdict(set)
    for locus in locus_gene_distance:
        distance_gene = []
        for gene in locus_gene_distance[locus]:
            distance_gene.append((locus_gene_distance[locus][gene], gene))
        distance_gene = sorted(distance_gene)
        mindist = distance_gene[0][0]
        locus_mindsistgeneset[locus] = set(g for (d, g) in distance_gene if (d, g)[0]==mindist)

    with open(results_dir + '/locus_mindsistgeneset.pickle', 'wb') as handle:
        pickle.dump(locus_mindsistgeneset, handle, protocol=pickle.HIGHEST_PROTOCOL)

    mindist_genes = [*locus_mindsistgeneset.values()]
    mindist_genes = set.union(*mindist_genes)

    networkgenes = get_networkgenes(locus_geneset)
    logger.info('network has %d SNPs, %d loci, and %d genes', len(
        snp_gene_distance), len(locus_geneset), len(networkgenes))

    gene_locus = defaultdict()
    for locus in locus_geneset:
        for gene in locus_geneset[locus]:
            gene_locus[gene] = locus

    logger.info(f'{len(gene_locus)} genes in gene_locus')
    logger.info(f'{len(mindist_genes)} mindist genes')
            

    snp_bp_dict = defaultdict(dict)
    for phenotype in phenotype_list:
        SNPfile =  os.path.join(project_datadir, phenotype + snploc_filename)
        ph_snp_bp = utils_data.read_snp_file(SNPfile)
        snp_bp_dict.update(ph_snp_bp)

    snp_bp_df = pd.DataFrame.from_dict(snp_bp_dict).transpose()
    snp_bp_df.reset_index(inplace=True)
    snp_bp_df = snp_bp_df.rename(columns = {'index':'rsid'})
    filename = os.path.join(results_dir, 'snp_bp.csv')

    snp_bp_df.to_csv(filename, index=False, header=True)

    #utils_plot.plot_locus_ngenes_vs_width(locus_geneset, GENE_BP_DICT, gene_special, results_dir)
    locus_stats = get_locus_stats(snp_bp_dict, locus_geneset, GENE_BP_DICT, results_dir)

    snp_filename = os.path.join(results_dir, 'snp_to_genes.txt')
    write_snp_genes(snp_filename, snp_gene_distance, gene_phenoset, snp_bp_dict)

    write_loci_aggregatedstats(results_dir, locus_stats)

    filename = os.path.join(results_dir, 'locus_to_genes.txt')
    write_locus_genes(filename, locus_gene_distance, locus_stats, GENE_BP_DICT, gene_phenoset, gene_type)

    logger.info(f'finished contructing loci')
    logger.info(f'gathering network edges')

    network_ppi = None
    xnull = None
    gene_PPIdegree = None
    ppi_Ksquared = None

    # implement PPI interactions in the Lk network
    # fnull is the null hypothesis for the probability that a pair of proteins has an interaction
    # do we just count interactions between loci? or also within loci?
    # might as well count them all and divide by number of pairs

    network_ppi = None
    networkfile = os.path.join(results_dir, 'network_ppi.txt')
        
    if os.path.isfile(networkfile):
        logger.info(f'reading network_ppi from {networkfile}')
        network_ppi = read_networkfile(networkfile)

    else:
        logger.info(f'regenerating network_ppi, takes a minute')

        logger.info(f'reading {protein_gene_file}')
        uniprot_to_hgnc = utils_ppi.protein_gene(protein_gene_file)
        uniprot_subset = get_uniprot_subset(uniprot_to_hgnc, networkgenes)
        logger.info(f'{len(networkgenes)} map to {len(uniprot_subset)} uniprot')
        network_ppi = get_network_ppi(uniprot_to_hgnc, uniprot_subset, networkgenes, ppi_file)
        logger.info(f'writing network_ppi to {networkfile} for next time')
        write_networkfile(networkfile, network_ppi)

    with open(results_dir + '/network_ppi.pickle', 'wb') as handle:
        pickle.dump(network_ppi, handle, protocol=pickle.HIGHEST_PROTOCOL)

    max_networkedges = utils_ppi.count_undirected_within_set(networkgenes, network_ppi)
    logger.info('%d max possible PPI edges in network', max_networkedges)

    # no degree_correction:
    fnull = utils_ppi.get_fnull(networkgenes, network_ppi)
    v = float(len(locus_geneset))  # pick only one gene per locus
    xnull = 0.5 * v * (v - 1) * fnull
    
    m = utils_ppi.get_Enull(locus_geneset, network_ppi)
    
    nlocus = len(locus_geneset)
    logger.info('%d loci', nlocus)
    logger.info('%f degree (null)', fnull * (nlocus - 1.0))
    logger.info('%f edges (null)', 0.5 * fnull * nlocus * (nlocus - 1.0))
    
    # degree_correction:
    
    gene_PPIdegree, gene_connectedPPIgenes= utils_ppi.get_degree_connectedgenes(networkgenes, network_ppi, gene_locus)

    Dmax_squared = 0
    for gene in gene_PPIdegree:
        Dmax_squared += gene_PPIdegree[gene]**2
        
    ppi_Ksquared = 2 * max_networkedges
        
    D1max, D2max = utils_ppi.get_D_other(networkgenes, gene_PPIdegree)
    xnull_max = (1 / (2 * ppi_Ksquared)) * \
        ((D1max)**2 - (D2max))  # ==nedge_max]]

    # network_ppi should have the interactions, network_ppi[p] is the set of neighbors of p
    logger.info(f'network_ppi has ...')
    logger.info(f'... {len(network_ppi)} proteins compared to {len(networkgenes)} in network')

    logger.info(f'describing network_ppi')
    describe_network(network_ppi)

    logger.info(f'done with network_ppi')

    logger.info(f'getting network_gri')
    
    TF_Ksquared = None
    TFsource_target = None
    reverse_TF = None
    gene_TFoutdegree = None
    gene_TFindegree = None

    xnull = True
    UniverseTFsource_target, UniverseTFsource_target_sign = utils_ppi.Trrust_dict(trrust_file)

    TFsource_target = dict()
    for p in UniverseTFsource_target:
        if p in networkgenes:
            TFsource_target[p] = UniverseTFsource_target[p].intersection(networkgenes)
            TFsource_target[p].discard(p)
            if len(TFsource_target[p]) == 0:
                del TFsource_target[p]

    with open(results_dir + '/TFsource_target.pickle', 'wb') as handle:
        pickle.dump(TFsource_target, handle, protocol=pickle.HIGHEST_PROTOCOL)

    TFsource_target_sign = dict()
    for s in TFsource_target:
        TFsource_target_sign[s] = dict()
        for t in TFsource_target[s]:
            TFsource_target_sign[s][t] = UniverseTFsource_target_sign[s][t]

    reverse_TF = utils_ppi.reverse_network(TFsource_target)

    # for symmetry with ppi name as network_gri and network_gri_rev
    network_gri = TFsource_target
    network_gri_rev = reverse_TF

    logger.info(f'network_gri statistics:')
    describe_network(network_gri)
    logger.info(f'network_gri_rev statistics:')
    describe_network(network_gri_rev)

    logger.info(f'writing network_gri files')
    write_networkfile(os.path.join(results_dir, 'network_gri.txt'), network_gri)
    write_networkfile(os.path.join(results_dir, 'network_gri_rev.txt'), network_gri_rev)

    degfac_ppi = get_degfac_ppi(network_ppi)
    with open(os.path.join(results_dir, 'degfac_ppi.txt'), 'w') as fp:
        fp.write('VERTEX\tDEGFAC\n')
        for (v1, degfac) in sorted(degfac_ppi.items()):
            fp.write(f'{v1}\t{degfac}\n')

    (degfac_gri, degfac_gri_rev) = get_degfac_gri(network_gri, network_gri_rev)
    with open(os.path.join(results_dir, 'degfac_gri.txt'), 'w') as fp:
        fp.write('VERTEX\tDEGFAC\tDEGFAC_REV\n')
        myset = set(degfac_gri.keys()).union(set(degfac_gri_rev.keys()))
        for v in sorted(myset):
            fp.write(f'{v}\t{degfac_gri.get(v, 0)}\t{degfac_gri_rev.get(v, 0)}\n')
            
    gene_TFindegree, gene_connectedTFIngene = utils_ppi.get_degree_connectedgenes(networkgenes, reverse_TF, gene_locus)
    gene_TFoutdegree, gene_connectedTFOutgene = utils_ppi.get_degree_connectedgenes(
        networkgenes, TFsource_target, gene_locus)

    max_TFedges = utils_ppi.count_directed_within_set(networkgenes, TFsource_target)
    logger.info('%d max possible TF edges in network', max_TFedges)

    Din_max, Dout_max, DinDout_max = utils_ppi.get_directedD_other(
        networkgenes, gene_TFoutdegree, gene_TFindegree)

    TF_Ksquared = (1 / (max_TFedges)) * \
        (Din_max * Dout_max - DinDout_max)

    TFxnull_max = (1 / (TF_Ksquared)) * ((Din_max * Dout_max) -
                                         (DinDout_max))  # equal to max_TFedges


    seedless_base_filename = 'ppi' + str(ppi_run) + '_TF' + str(TF_run) + '_degcor' + str(degree_correction) + '_optim' + str(optim_run) +'_unannotated' + str(unannotated) + '_pseudogene' + str(pseudogene) + '_no-antisense' + str(antisense)
    seedless_base_filename = 'SIGNET'

    counts_log = os.path.join(results_dir, seedless_base_filename + '_counts.txt')
    fp = open(counts_log, 'w')
    fp.write('\t'.join(['', 'seed:'])+ '\n')
    fp.write('\t'.join(['', 'network_score:'])+ '\n')
    for locus in sorted(locus_geneset.keys()):
        for gene in sorted(locus_geneset[locus]):
            fp.write('\t'.join([locus, gene]) + '\n')
    fp.close()


    locus_nppigene = dict()
    for locus in locus_geneset:
        locus_nppigene[locus] = 0
        for gene in locus_geneset[locus]:
            if gene in network_ppi:
                locus_nppigene[locus]+=1

    with open(results_dir + '/locus_nppi.txt', 'w') as f:
        f.write('\t'.join(['locus_name', 'nppi']) + '\n')
        for locus in locus_nppigene:
            f.write('\t'.join([locus, str(locus_nppigene[locus])]) + '\n')

    with open(results_dir + '/gene_PPI_GRI.txt', 'w') as f:
        f.write('\t'.join(['gene_name', 'PPIdegree', 'TF_outdeg', 'TF_indeg']) + '\n')
        for locus in sorted(locus_geneset.keys()):
            for gene in sorted(locus_geneset[locus]):
                f.write('\t'.join([gene, str(gene_PPIdegree[gene]), str(gene_TFoutdegree[gene]), str(gene_TFindegree[gene])]) + '\n')

    logger.info('*** FINISHED NETWORK CONSTRUCTION ***')

    init_countsdf = pd.DataFrame()
    network_locus_list = []
    network_gene_list = []
    for locus in sorted(locus_geneset.keys()):
        for gene in sorted(locus_geneset[locus]):
            network_locus_list.append(locus)
            network_gene_list.append(gene)
    init_countsdf['locus_name'] = network_locus_list
    init_countsdf['gene_name'] = network_gene_list

    allruns_countlog = os.path.join(results_dir, seedless_base_filename + '_counts_allruns.txt')
    if not os.path.exists(allruns_countlog):
        init_countsdf.to_csv(allruns_countlog, index=False, sep='\t')

    for seed in range(seed_start, seed_end):

        logger.info(f'*** starting seed {seed} ***')
        #base_filename = 'seed' + str(seed) + '_ppi' + str(ppi_run) + '_TF' + str(TF_run) + '_degcor' + str(degree_correction) + '_optim' + str(optim_run) +'_unannotated' + str(unannotated) + '_pseudogene' + str(pseudogene) + '_no-antisense' + str(antisense)
        base_filename = f'seed{seed:02d}'

        run_experiment(phenotype_list, k, FLANKDIST, degree_correction, optim_run, network_plot,
                       gene_PPIdegree, ppi_Ksquared, network_ppi, xnull, ppi_run,
                       TF_run, TFsource_target, reverse_TF, TF_Ksquared,
                       gene_TFindegree, gene_TFoutdegree, locus_geneset, networkgenes,
                       gene_special, locus_mindsistgeneset, gene_distance, gene_signeddistance,
                       gene_type, locus_snptype, external_genes, seed,
                       results_dir, seedless_base_filename, counts_log,
                       unannotated, pseudogene, antisense,
                       degfac_ppi, network_gri, network_gri_rev, degfac_gri, degfac_gri_rev, init_rand)

        logger.info(f'... finished seed {seed}')

        
        counts_df = pd.read_csv(allruns_countlog, sep='\t')
        file = os.path.join(results_dir, base_filename + '_extended_finalpass_scores.txt')




        pass_df = pd.read_csv(file, sep='\t')
        loci = counts_df['locus_name']



        gene_selection_list = []
        for locus in np.unique(loci):
            locus_df = pass_df[pass_df['locus_name'] == locus]
            selectedgene_row = locus_df['gene_pattern'] !='.'
            gene_selection = selectedgene_row.astype(int).values.tolist()
            gene_selection_list.extend(gene_selection)



        counts_df[str(seed)] = gene_selection_list
        counts_df.to_csv(allruns_countlog, index=False, sep='\t')


    #finished seed_end iterations
    #calculate prob_selected
    file = os.path.join(results_dir, seedless_base_filename + '_counts_allruns.txt')
    counts_df = pd.read_csv(file, sep='\t')
    nselected = counts_df.iloc[:, 2:].sum(axis=1) # number times the gene was selected
    nseed = len(counts_df.columns) - 2  # number of 'other' columns
    counts_df['gene_prob_selected'] = (nselected / nseed)

    counts_filename = os.path.join(results_dir, seedless_base_filename + '_counts_allruns_prob.txt')
    counts_df.to_csv(counts_filename, index=False, sep='\t')


    logger.info(f'*** concluded successfully ***')

    return None

if __name__ == "__main__":
    main()
