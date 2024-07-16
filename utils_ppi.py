#!/usr/bin/env python
# coding: utf-8


import numpy as np
import collections
import copy
import csv
import pdb

from collections import Counter


import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def count_edges(network):
    # counts directed edges
    n = 0
    for a in network.keys():
        for b in network[a]:
            n += 1
    return(n)


def make_undirected(network):
    ret = copy.deepcopy(network)
    m = 0
    for a in network.keys():
        for b in network[a]:
            if (b == a):
                ret[a].remove(b)
                m += 1
            else:
                ret[b].update([a])

    n1 = count_edges(network)
    n2 = count_edges(ret)
    logger.info(
        'directed network with %d edges -> undirected network with %d edges, %d', n1, n2, m)
    return(ret)


def read_ppi(filename):
    ret = collections.defaultdict(set)
    m = 0
    with open(filename, newline='') as fp:
        reader = csv.DictReader(fp, delimiter='\t')
        for row in reader:
            ppi_type = row['evidence_type']
            if ('exp' in ppi_type) or ('ortho' in ppi_type):
                (a, b) = (row['uniprot1'], row['uniprot2'])
                ret[a].update([b])
            else:
                m += 1

    fp.close()
    n = count_edges(ret)
    logger.info(
        '%s -> %d directed PPI interactions and %d removed_pred_edges', filename, n, m)
    return(ret)

def read_ppi_subset(filename, uniprot_subset):
    ret = collections.defaultdict(set)
    m = 0
    with open(filename, newline='') as fp:
        reader = csv.DictReader(fp, delimiter='\t')
        for row in reader:
            (a, b, evidence) = (row['uniprot1'], row['uniprot2'], row['evidence_type'])
            if (a in uniprot_subset) and (b in uniprot_subset):
                if ('exp' in evidence) or ('ortho' in evidence):
                    ret[a].update([b])
            else:
                m += 1

    fp.close()
    n = count_edges(ret)
    logger.info(
        '%s -> %d directed PPI interactions and %d skipped', filename, n, m)
    return(ret)

# network is a dict of sets, vertex to its neighbors
# u_to_h is a dict of sets, uniprot to hgncs


def convert_uniprot_to_hgnc(network, u_to_h):
    ret = collections.defaultdict(set)
    for u1 in network:
        if u1 not in u_to_h:
            # logger.info('no uniprot for %s', u1)
            continue
        set1 = u_to_h[u1]
        for u2 in network[u1]:
            if u2 not in u_to_h:
                # logger.info('no uniprot for %s', u2)
                continue
            set2 = u_to_h[u2]
            for a in set1:
                if a not in ret:
                    ret[a] = set()
                for b in set2:
                    ret[a].add(b)
    return(ret)


def gene_protein(filename):
    # Find mapping of protein name to Gene names
    fp = open(filename, 'r')
    ret = collections.defaultdict(list)
    next(fp)  # skip header line
    for line in fp:
        toks = line.strip().split()

        (a, b) = (toks[0], toks[2:])
        for gene in b:
            ret[gene].append(a)
    fp.close()
    return(ret)


def protein_gene(filename):
    fp = open(filename, 'r')
    ret = collections.defaultdict(list)
    next(fp)
    for line in fp:
        toks = line.strip().split()
        (a, b) = (toks[0], toks[2:])
        for gene in b:
            ret[a].append(gene)
    fp.close()
    return (ret)


def gene_gene(gene_list, protein_gene_dict, gene_protein_dict, protein_protein_dict):
    gene_gene_dict = collections.defaultdict(list)
    for gene in gene_list:
        protein1_list = []
        protein1_list.extend(gene_protein_dict[gene])
        protein2_list = []
        for x in protein1_list:
            protein2_list.extend(protein_protein_dict[x])
        gene2_set = set()
        for y in protein2_list:
            gene2_set = gene2_set.union(set(protein_gene_dict[y]))
        gene2_nw = set(gene2_set).intersection(set(gene_list))
        gene2_nw.discard(gene)
        if len(gene2_nw) > 0:
            gene_gene_dict[gene].extend(list(gene2_nw))

    return (gene_gene_dict)


# how many undericted edges within a set of vertices?
def count_undirected_within_set(myset, network_interactions):
    n = 0
    for a in myset:

        if a not in network_interactions:
            continue
        for b in network_interactions[a]:
            if (b < a):
                continue
            if b not in myset:
                continue
            n += 1
    return(int(n))

# how many dericted edges within a set of vertices?


def count_directed_within_set(myset, network_interactions):
    n = 0
    for a in myset:

        if a not in network_interactions:
            continue
        for b in network_interactions[a]:
            if b not in myset:
                continue
            n += 1
    return(int(n))


def count_vertex_to_set(v, myset, network_interactions):
    n = 0
    if v in network_interactions:
        for a in network_interactions[v]:
            if a in myset:
                n += 1
    return(n)


def vertex_to_set(v, myset, network_interactions):
    connected = set()
    if v in network_interactions:
        for a in network_interactions[v]:
            if a in myset:
                connected.add(a)
    return(connected)


# calculate fraction of occupied edges between genes for an undirected network
# don't include self-edges
def get_fnull(genes, network_interactions):
    logger.info('calculating fnull')
    nedge = 0
    mylist = sorted(genes)
    for v in mylist:
        if v not in network_interactions:
            continue
        for w in mylist:
            # skip self-edges
            if (w >= v):
                continue
            if w in network_interactions[v]:
                nedge += 1

    nvert = len(mylist)
    npair = float(0.5 * nvert * (nvert - 1.0))
    fnull = float(nedge) / float(npair)
    # the average degree should be (fnull)(nvert - 1)
    avgdegree = fnull * (nvert - 1.0)
    logger.info('ngenes %d nedge %d npair %d fnull %f avgdeg %f',
                nvert, nedge, npair, fnull, avgdegree)
    return(fnull)


def get_Enull(locus_geneset, network_interactions):
    logger.info('calculating Enull')

    Enull = 0
    loci = sorted(locus_geneset.keys())
    for locusX in loci:
        if locusX.startswith('loc'):  # ONLY FOR DEBUGGING
            continue

        ngenes_locusX = len(locus_geneset[locusX])

        for v in sorted(locus_geneset[locusX]):
            if v not in network_interactions:
                continue

            for locusY in loci:  # find v-w connections in network, weighted by loci ngene)
                if locusY >= locusX:
                    continue
                if locusY.startswith('loc'):  # ONLY FOR DEBUGGING
                    continue

                v_locusY_nedge = 0
                ngenes_locusY = len(locus_geneset[locusY])

                for w in sorted(locus_geneset[locusY]):
                    if w in network_interactions[v]:
                        v_locusY_nedge += 1

                Enull = Enull + (v_locusY_nedge /
                                 (ngenes_locusX * ngenes_locusY))

    return(Enull)


def get_degree_connectedgenes(geneset, network_interactions, gene_locus):
    logger.info('computing degree of genes')

    gene_connectedgenes = collections.defaultdict(set)
    gene_degree = dict()
    mylist = sorted(geneset)
    for v in mylist:
        locus_set = set()
        ndegree = 0
        if v not in network_interactions:
            gene_degree[v] = 0
            continue
        v_network_interactions = network_interactions[v].intersection(mylist)
        for w in v_network_interactions:
            gene_connectedgenes[v].update([w])
            w_locus = gene_locus[w]

            if not (w_locus in locus_set):
                locus_set = locus_set.union(set([w_locus]))
                ndegree += 1

        gene_degree[v] = int(ndegree)
    return gene_degree, gene_connectedgenes


def get_D_other(otherset, gene_degree):
    D1 = 0
    D2 = 0
    for g in otherset:
        D1 += gene_degree[g]
        D2 += gene_degree[g]**2
    return D1, D2


def get_lambda_gene(gene, gene_degree, D1, D2, c):
    f_g = gene_degree[gene]
    lambda_gene = (1 / (2 * c)) * ((D1 + f_g)**2 - (D2 + f_g**2))
    return lambda_gene


def Trrust_dict(filename):
    logger.info(f'reading Trrust data from {filename}')
    fp = open(filename, 'r')
    TFsource_target = dict()
    TFsource_target_sign = dict()
    next(fp)  # skip header line
    for line in fp:
        toks = line.strip().split()
        assert(len(toks) == 4), 'bad line: ' + line
        (a, b, c) = (toks[0], toks[1], toks[2])
        if not a in TFsource_target:
            TFsource_target[a] = set()
            TFsource_target_sign[a] = dict()
        TFsource_target[a].update([b])
        TFsource_target_sign[a][b] = c
    fp.close()
    return (TFsource_target, TFsource_target_sign)


def reverse_network(network_interactions):
    logger.info('reversing network direction')
    reverse_interactions = collections.defaultdict(set)
    for v in network_interactions:
        for w in network_interactions[v]:
            reverse_interactions[w].update([v])
    return reverse_interactions


def get_directedD_other(otherset, gene_outdegree, gene_indegree):
    Din = 0
    Dout = 0
    DinDout = 0

    for g in otherset:
        Din += int(gene_indegree[g])
        Dout += int(gene_outdegree[g])
        DinDout += int(gene_indegree[g] * gene_outdegree[g])
    return Din, Dout, DinDout


def get_directedlambda_gene(gene, gene_outdegree, gene_indegree, Din, Dout, DinDout, c):
    gene_din = gene_indegree[gene]
    gene_dout = gene_outdegree[gene]
    lambda_gene = (1 / c) * (gene_din * Dout +
                             gene_dout * Din + Din * Dout - DinDout)

    return lambda_gene
