#!/usr/bin/env python
# coding: utf-8


import numpy as np
import csv
from collections import defaultdict
import warnings
from collections import Counter
from operator import itemgetter
import pandas as pd
import os
import re
import pdb
import yaml
import glob

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def isInt(mystr):
    try:
        int(mystr)
    except ValueError:
        return False
    else:
        return True

def read_gene_file(gene_file, unannotated, pseudogene, antisense, project_datadir, external_genes):
    antisense_skipped = []
    unannotated_skipped = []
    pseudo_skipped = []

    ret = defaultdict(dict)
    synonym_ret = defaultdict()
    notIntSet = set()
    notIntCnt = 0
    with open(gene_file) as fp:

        reader = csv.DictReader(fp, delimiter='\t')

        nrow = 0
        for row in reader:
            nrow += 1

            if (row['Gene stable ID'] == '[success]'):
                # this is a success indicator from the wget command downloading the file
                # should be the last line of the file if it's here
                # we'll continue rather than check
                logger.info(f'... skipping [success], row number {nrow}')
                continue

            (gene, chromosome, gene_tss, gene_end, gene_start, gene_type, gene_description, gene_synonym) = (
                row['Gene name'], row['Chromosome/scaffold name'], row['Transcription start site (TSS)'], row['Gene end (bp)'], row['Gene start (bp)'], row['Gene type'], row['Gene description'], row['Gene Synonym'])

            # if chromosome.startswith('CHR'):
            #   continue

            # require the chromosome to be an integer
            # skip X, Y, mitochondrial, patches, ...
            if not isInt(chromosome):
                notIntSet.add(chromosome)
                notIntCnt += 1
                continue

            synonym_ret[gene_synonym] = gene

            if unannotated == False:
                if re.search(r'\A[A][a-zA-Z]\d{6}', gene):
                    unannotated_skipped.append(gene)
                    continue

            accepted_types = set({'omim', 'coloc', 'protein_coding'})
            if pseudogene == False:
                if (((gene.startswith("RNU")) or gene_type not in accepted_types) and gene not in external_genes):
                    pseudo_skipped.append(gene)
                    continue
                    #{'SSXP10', 'GBAP1', 'LINC00881'} gene_type is IndRNA but they are in coloc geneset, so we keep them

            if antisense == True:
                if "-AS" in gene:
                    antisense_skipped.append(gene)
                    continue

            if gene not in ret:
                ret[gene]['chromosome'] = chromosome

                # distances might be stored as strings
                # convert to integers
                igs = int(gene_start)
                ige = int(gene_end)
                itss = int(gene_tss)
                
                # gene_start and gene_end are the position on the chromosome
                # start should always be smaller than end
                
                assert(igs <= ige), f'bad gene endpoints: {row}'

                ret[gene]['gene_start'] = gene_start
                ret[gene]['gene_end'] = gene_end

                # change to using the 5' end of the gene rather than the tss of the first transcript
                # also keep track of strand
                mystrand = '+' if abs(igs-itss) <= abs(ige-itss) else '-'
                ret[gene]['gene_strand'] = mystrand
                ret[gene]['gene_tss'] = gene_start if (mystrand == '+') else gene_end
                
                
                #ret[gene]['gene_tss'] = gene_tss
                #if gene_tss == gene_end:
                #    ret[gene]['gene_tss'] = gene_start
                #    if gene_start == gene_end:
                #        logger.info('gene %s has width 0. check data', gene)
                
                ret[gene]['gene_type'] = gene_type
                ret[gene]['gene_description'] = gene_description
                ret[gene]['gene_synonym'] = gene_synonym

    fp.close()

    logger.info(f'skipped {notIntCnt} rows mapping to {len(notIntSet)} non-integer chromosomes')

    ensembl_genes = set(ret.keys())
    missing_genes = external_genes - ensembl_genes
    for gene in missing_genes:
        if gene in synonym_ret:
            if synonym_ret[gene] in ret:
                ret[gene] = ret.pop(synonym_ret[gene])

            else:
                logger.info(
                    'gene(s) %s not found in ensemble_gene or ensemble_gene_synonym', missing_genes)

    n_unannotated_skipped = len(np.unique(unannotated_skipped))
    n_pseudo_skipped = len(np.unique(pseudo_skipped))
    n_antisense_skipped = len(np.unique(antisense_skipped))

    logger.info('found and skipped %d unannotated genes',
                n_unannotated_skipped)
    logger.info('skipped %d pseudo genes', n_pseudo_skipped)
    logger.info('skipped %d antisense genes', n_antisense_skipped)

    return(ret)


READ_SNP_FILE_CACHE = dict()


def read_snp_file(file, delimiter='\t'):

    if file in READ_SNP_FILE_CACHE:
        logger.info(f'returning cached version of {file}')
        return(READ_SNP_FILE_CACHE[file])

    fp = open(file, 'r')
    ret = defaultdict(dict)

    next(fp)  # this skips the first line (header line)

    for line in fp:
        toks = line.strip().split(delimiter)
        (rsid, chromosome, rsid_start, rsid_end) = (
            toks[0], toks[2], toks[3], toks[4])

        if chromosome.startswith('CHR'):
            continue

        if rsid in ret:
            continue  # for speed, since GWAS files can have repeat snp's

        ret[rsid]['chromosome'] = chromosome
        ret[rsid]['rsid_start'] = rsid_start
        ret[rsid]['rsid_end'] = rsid_end

    fp.close()
    READ_SNP_FILE_CACHE[file] = ret
    return(ret)


def read_gwas_file(data_dir, gwas_filename, phenotype_list, delimiter='\t'):

    for phenotype in phenotype_list:
        file = os.path.join(data_dir, phenotype + gwas_filename)
        ret = dict()

        logger.info(f'reading GWAS file {file}')
        with open(file) as ifile:
            reader = csv.DictReader(ifile, delimiter=delimiter)
            for row in reader:
                rsid = row['SNPS']  # * 10 ^ (-8)
                if ((rsid not in ret and float(row['P-VALUE']) <= 5 * 10**(-8)) or ((rsid in ret) and (float(ret[rsid]['P']) > float(row['P-VALUE'])))):
                    ret[rsid] = defaultdict(list)
                    ret[rsid]['pubmed'] = row['PUBMEDID']
                    ret[rsid]['type'] = row['CONTEXT']
                    # might have to revisit this parsing to ensure the numbers are correct
                    # not used currently
                    sample_description = row['INITIAL SAMPLE SIZE'].replace(",", "")
                    # ret[rsid]['Neff'] = re.search('\d+|$', sample_description).group()
                    ret[rsid]['Neff'] = sample_description
                    ret[rsid]['P'] = "{:.50f}".format(float(row['P-VALUE']))

    return(ret)


def locus_snptype(locus_geneset, data_dir, gwas_filename, phenotype_list, delimiter='\t'):

    snp_type = read_gwas_file(
        data_dir, gwas_filename, phenotype_list, delimiter)

    ret = defaultdict()
    for locus in locus_geneset:

        locus_snptype = []
        locus_snps = locus.split('+')
        locus_snptype = [snp_type[snp]['type']
                         for snp in locus_snps if snp in snp_type]
        locus_snptype_string = '+'.join(sorted(locus_snptype))
        ret[locus] = locus_snptype_string

    return(ret)


def get_snp_gene_distance(gene_bp_dict, snp_bp_dict):
   # calculate distance to all genes on same chromosome
    ret = defaultdict(list)

    for variant in snp_bp_dict.keys():

        variant_list = []

        # find position of variant
        variant_chromosome = snp_bp_dict[variant]['chromosome']
        variant_bp = int(snp_bp_dict[variant]['rsid_start'])
        # find distance to all genes in that chromosome:
        for gene in gene_bp_dict.keys():
            if not gene_bp_dict[gene]['chromosome'] == variant_chromosome:
                continue
            else:
                gene_tss = int(gene_bp_dict[gene]['gene_tss'])

                # this is incorrect, want the signed distance based on strand
                # if plus strand, then distance = variant - tss
                # if minus strand, then distance = tss - variant
                # negative means snp is upstream, positive means snp is downstream
                distance = (gene_tss - variant_bp)

            variant_list.append((gene, distance))

        ret[variant] = sorted(variant_list, key=lambda x: abs(x[1]))

    return (ret)


def get_network(k, FLANKDIST, gene_bp_dict, phenotype, data_dir, snploc_filename):

    # SNP file is the SNP information from ensembl
    snploc_file = os.path.join(data_dir, phenotype + snploc_filename)
    snp_bp_dict = read_snp_file(snploc_file)

    # for each variant, a list of all genes in the chromosome sorted by distance
    phen_dict_raw = get_snp_gene_distance(gene_bp_dict, snp_bp_dict)

    # if want to limit to the closest k genes:
    # phen_dict[variant] = phen_dict[variant][0:k]

    phen_dict = dict()
    for variant in phen_dict_raw.keys():
        # closest k genes
        variant_genedist = phen_dict_raw[variant][0:k]
        # add more genes based on distance of kth gene
        for gene, distance in phen_dict_raw[variant][k:]:

            if abs(distance) < FLANKDIST:
                variant_genedist.append((gene, distance))
            else:
                break
        phen_dict[variant] = variant_genedist

    # get unique genes
    uniq_genes = dict()
    for v in phen_dict:
        for g in phen_dict[v]:
            uniq_genes[g] = uniq_genes.get(g, 0) + 1
    gene_hist = dict()
    for g in uniq_genes:
        cnt = uniq_genes[g]
        gene_hist[cnt] = gene_hist.get(cnt, 0) + 1

    logger.info(
        f'phenotype {phenotype}: {len(phen_dict)} variants, {len(uniq_genes)} genes, {str(gene_hist)}')

    return phen_dict


def get_feature_geneset(datadir, phenotype_list, delimiter='\t', encoding=None):
    feature_geneset = defaultdict(set)
    for file in glob.glob(datadir + '/*_genes.txt'):
        filename = os.path.basename(file)
        splits = filename.split('_')

        if len(splits) == 3:
            if not any(phenotype in splits[0] for phenotype in phenotype_list):
                continue

        feature = filename.split('_')[-2]

        df = pd.read_csv(file, delimiter=delimiter, encoding=encoding)
        geneset = set(df['gene_name'].values)
        feature_geneset[feature] = feature_geneset[feature].union(geneset)
    return (feature_geneset)


def get_geneset(data_dir, filename, delimiter='\t', encoding=None):
    geneset = set()
    file = os.path.join(data_dir, filename)
    df = pd.read_csv(file, delimiter=delimiter, encoding=encoding)
    geneset = set(df['gene_name'].values)
    return (geneset)


def omim_suptable(omim_file, omim_suptable_file):
    df = pd.read_csv(omim_file, sep='\t')
    df = df.replace(np.nan, '')
    dup = df[df.duplicated(subset=['gene_name'], keep=False)]
    for g in dup['gene_name']:
        syndromes = ','.join(
            sorted(df[df['gene_name'] == g]['Syndrome'].values))
        subtypes = ','.join(
            filter(None, sorted(df[df['gene_name'] == g]['Syndrome Subtype'].values)))
        df = df.append({'gene_name': g, 'Syndrome': syndromes,
                        'Syndrome Subtype': subtypes}, ignore_index=True)
        df.drop_duplicates(subset=['gene_name'], keep='last')

    df = df.sort_values(by=['gene_name'])
    df.to_csv(omim_suptable_file, sep='\t', index=False)


def get_omim_pheno(omim_file):
    fp = open(omim_file, 'r')
    next(fp)  # this skips the first line (header line)
    ret = defaultdict(set)
    for line in fp:
        toks = line.strip().split('\t')
        gene, pheno = toks[0], toks[1]
        ret[gene].add(pheno)
    return (ret)


def get_gene_featureset(gene_distance, feature_geneset):
    gene_featureset = defaultdict(set)

    for feature in feature_geneset:
        for gene in feature_geneset[feature]:
            gene_featureset[gene].add(feature)

    for gene in gene_distance:
        if gene not in gene_featureset:
            gene_featureset[gene] = set()

    return (gene_featureset)


def get_gene_featureset_onespecialblocked(gene_distance, feature_geneset, blockedspecialgene):
    gene_featureset = defaultdict(set)

    for feature in feature_geneset:
        for gene in feature_geneset[feature]:
            gene_featureset[gene].add(feature)

    for gene in gene_distance:
        if gene not in gene_featureset:
            gene_featureset[gene] = set()

    gene_featureset[blockedspecialgene] = set()
    return (gene_featureset)


def independent_snps_exclude_genes(independent_snps, snp_gene_distance):
    #independent_snps = {'rs10040989', 'rs1468333'}
    # for snp in sorted(independent_snps):
    #     print(snp)
    #     print(sorted(snp_gene_distance[snp].items()))

    gene_indsnp_distance = defaultdict(dict)
    for snp in sorted(independent_snps):

        for gene in sorted(snp_gene_distance[snp]):
            if gene not in gene_indsnp_distance:
                gene_indsnp_distance[gene]['snp'] = snp
                gene_indsnp_distance[gene]['dist'] = snp_gene_distance[snp][gene]

            else:

                registered_distance = gene_indsnp_distance[gene]['dist']
                current_snp_distance = snp_gene_distance[snp][gene]

                if current_snp_distance < registered_distance:
                    registered_snp = gene_indsnp_distance[gene]['snp']
                    del(snp_gene_distance[registered_snp][gene])
                    gene_indsnp_distance[gene]['snp'] = snp
                    gene_indsnp_distance[gene]['dist'] = snp_gene_distance[snp][gene]
                else:
                    del(snp_gene_distance[snp][gene])

    # # print(gene_indsnp_distance)
    # for snp in sorted(independent_snps):
    #     print(snp)
    #     print(sorted(snp_gene_distance[snp].items()))
    # pdb.set_trace()

    return snp_gene_distance


def merge_pheno_snps(snp_to_genes):
    # snp_to_genes is a dictionary, key = snp, value = list of genes
    # output should be similar, but merging SNPs whose genes overlap
    # the important thing is that if v is a value, then we can build a dictionary seen(v)
    # and i don't know if that will work if v is a tuple
    # but also if a gene is assigned to multiple variants, the distances could be different
    # so we really just want the gene names
    locus_to_genes = defaultdict(list)  # genes as a list
    gene_to_locus = dict()  # gene points to a single locus
    for snp in snp_to_genes:
        # i think this is the same as genes, distances = zip(*snp_to_genes) but i'm never sure
        genes = [x[0] for x in snp_to_genes[snp]]
        distances = [x[1] for x in snp_to_genes[snp]]
        # figure out if any of the genes is already in a locus
        merge = dict()
        for g in genes:
            if g in gene_to_locus:
                merge[gene_to_locus[g][0]] = True

        if not merge:
            for i in range(len(genes)):
                g = genes[i]
                d = abs(distances[i])
                gene_to_locus[g] = (snp, d)
                locus_to_genes[snp].append((g, d))
        else:
            all_genes_dist = dict()
            for i in range(len(genes)):
                g = genes[i]
                d = abs(distances[i])
                all_genes_dist[g] = d

            toks = [snp]
            for m in merge:
                toks = toks + m.split('+')
                other_snp_genes = [x[0] for x in locus_to_genes[m]]
                for g in other_snp_genes:
                    d = gene_to_locus[g][1]
                    if g in all_genes_dist:
                        d = min(all_genes_dist[g], d)
                    all_genes_dist[g] = d
                del locus_to_genes[m]
            locusname = '+'.join(sorted(toks))
            locus_gene_distance = [(g, d) for (g, d) in all_genes_dist.items()]
            locus_to_genes[locusname] = sorted(
                locus_gene_distance, key=lambda x: x[1])

            for g, d in all_genes_dist.items():
                gene_to_locus[g] = (locusname, d)

    return locus_to_genes


def merge_loci(union_locus_geneset, union_locus_gene_distance):

    locus_to_genes = defaultdict(list)  # genes as a list
    gene_to_locus = dict()  # gene points to a single locus

    for locus in union_locus_geneset:
        genes = union_locus_geneset[locus]
        # figure out if any of the genes is already in a locus
        merge = dict()
        for g in genes:
            if g in gene_to_locus:
                merge[gene_to_locus[g][0]] = True  # merge: locus: True

        if not merge:
            for g in genes:
                d = abs(union_locus_gene_distance[locus][g])
                gene_to_locus[g] = (locus, d)
                locus_to_genes[locus].append((g, d))

        else:
            all_genes_dist = dict()
            for g in genes:

                d = abs(union_locus_gene_distance[locus][g])
                all_genes_dist[g] = d

            toks = [locus]

            for m in merge:
                toks = toks + m.split('+')
                other_snp_genes = [x[0] for x in locus_to_genes[m]]
                for g in other_snp_genes:
                    d = gene_to_locus[g][1]
                    if g in all_genes_dist:
                        d = min(all_genes_dist[g], d)
                    all_genes_dist[g] = d
                del locus_to_genes[m]
            snpset = set()

            for t in toks:
                for s in t.split('+'):
                    snpset.add(s)
            locusname = '+'.join(sorted(snpset))
            locus_gene_distance = [(g, d) for (g, d) in all_genes_dist.items()]
            locus_to_genes[locusname] = sorted(
                locus_gene_distance, key=lambda x: x[1])

            for g, d in all_genes_dist.items():
                gene_to_locus[g] = (locusname, d)

    return locus_to_genes


# def merge_loci_exclude(union_locus_geneset, union_locus_gene_distance):

#     locus_to_genes = defaultdict(list)  # genes as a list
#     gene_to_locus = dict()  # gene points to a single locus
#     independent_snps = {'rs10040989', 'rs1468333', 'rs1321311', 'rs236349'}

#     pdb.set_trace()
#     for locus in {'rs10040989', 'rs1468333'}:  # union_locus_geneset:
#         genes = union_locus_geneset[locus]
#         # figure out if any of the genes is already in a locus
#         merge = dict()
#         for g in genes:
#             if g in gene_to_locus:
#                 merge[gene_to_locus[g][0]] = True  # merge: locus: True

#         pdb.set_trace()

#         if not merge:
#             for g in genes:
#                 d = abs(union_locus_gene_distance[locus][g])
#                 gene_to_locus[g] = (locus, d)
#                 locus_to_genes[locus].append((g, d))

#             pdb.set_trace()

#         else:
#             all_genes_dist = dict()
#             for g in genes:

#                 d = abs(union_locus_gene_distance[locus][g])
#                 all_genes_dist[g] = d

#             toks = [locus]

#             pdb.set_trace()

#             if merge.split('+').intersection(independent_snps) and locus in independent_snps: #these are independent snps/signals so do not merge but need to find which gene goes to which independent locus

    #         else:
    #             #not indepenent SNPs so merge
    #             for m in merge:
    #                 toks = toks + m.split('+')
    #                 other_snp_genes = [x[0] for x in locus_to_genes[m]]
    #                 for g in other_snp_genes:
    #                     d = gene_to_locus[g][1]
    #                     if g in all_genes_dist:
    #                         d = min(all_genes_dist[g], d)
    #                     all_genes_dist[g] = d
    #                 del locus_to_genes[m]
    #             snpset = set()

    #             pdb.set_trace()
    #             for t in toks:
    #                 for s in t.split('+'):
    #                     snpset.add(s)
    #             locusname = '+'.join(sorted(snpset))
    #             locus_gene_distance = [(g, d) for (g, d) in all_genes_dist.items()]
    #             locus_to_genes[locusname] = sorted(
    #                 locus_gene_distance, key=lambda x: x[1])

    #             pdb.set_trace()

    #             for g, d in all_genes_dist.items():
    #                 gene_to_locus[g] = (locusname, d)

    #             pdb.set_trace()

    # return locus_to_genes


def GWAScat_rows(phenotype_list, GWAS_dir, project_datadir):
    for phenotype in phenotype_list:
        filename = phenotype + '_GWAS_rows.tsv'
        inputfile = os.path.join(GWAS_dir, filename)

        df = pd.read_csv(inputfile, delimiter='\t')

        rsid_df = df[['SNPS']]
        rsid_list = df['SNPS'].values

        # extract variant rsid in a file to be uploaded on Ensembl
        outfile = GWAS_dir + phenotype + '_GWAScat_rsid.txt'
        with open(outfile, 'w', newline='') as myfile:
            wr = csv.writer(myfile)
            wr.writerow(rsid_list)

        # extract SNP data
        GWAScat_df = df[['SNPS', 'PUBMEDID', 'MAPPED_TRAIT', 'STRONGEST SNP-RISK ALLELE',
                         'P-VALUE', 'OR or BETA', '95% CI (TEXT)', 'CONTEXT']]

        GWAScat_df['Risk_Allele'] = [
            str.split(x, '-')[1] for x in GWAScat_df['STRONGEST SNP-RISK ALLELE'].values]
        logger.info('writing %s', filename)
        GWAScat_df.to_csv(path_or_buf=project_datadir +
                          phenotype + '_GWAScat.txt', sep='\t')
    return None


def get_emagmagenes(file_list, results_dir, project_dirname, delimiter=' '):

    ret = set()

    for file in file_list:

        fp = open(file, 'r')

        next(fp)  # this skips the first line (header line)

        for line in fp:
            gene = line.strip().split(delimiter)[0]
            ret.add(gene)

    fp.close()

    filename = os.path.join(
        results_dir, project_dirname + '_emagma_NCBIgenes.txt')

    fp = open(filename, 'w+')
    for gene in ret:
        fp.write(gene + '\n')
    fp.close()

    return(ret)
