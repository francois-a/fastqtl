#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import numpy as np
import pandas as pd
import os
import gzip
from datetime import datetime

parser = argparse.ArgumentParser(description='Filter significant SNP-gene pairs from FastQTL results using FDR cutoff')
parser.add_argument('permutation_results', help='FastQTL output')
parser.add_argument('fdr', type=np.double, help='False discovery rate (e.g., 0.05)')
parser.add_argument('annotation_gtf', help='Annotation in GTF format')
parser.add_argument('snp_lookup', help='Tab-delimited file with columns: Chr, Pos, VariantID, Ref_b37, Alt, RS_ID_dbSNP135_original_VCF, RS_ID_dbSNP142_CHG37p13, Num_alt_per_site')
parser.add_argument('--nominal_results', default='', help='FastQTL output from nominal pass')
parser.add_argument('--nominal_results_unnormalized', nargs=2, default='', help='FastQTL output file (nominal pass), and units')
args = parser.parse_args()

#------------------------------------------------------------------------------
# 1. eGenes (permutation output): add gene and variant information
#------------------------------------------------------------------------------
gene_dict = {}
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Parsing GTF ... ', end='', flush=True)
# add: gene_name, gene_chr, gene_start, gene_end, strand
with open(args.annotation_gtf) as gtf:
    for row in gtf:
        row = row.strip().split('\t')
        if row[0][0]=='#' or row[2]!='gene': continue
        attributes = row[8].split('; ',5)
        gene_id = attributes[0].split()[1].replace('"','')
        gene_name = attributes[4].split()[1].replace('"','')
        gene_dict[gene_id] = [gene_name, row[0], row[3], row[4], row[6]]
print('done.', flush=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Loading variant ID lookup table ... ', end='', flush=True)
snp_lookup_df = pd.read_csv(args.snp_lookup, index_col=2, sep='\t',
    dtype={'chr':str,'variant_pos':np.int64,  'variant_id':str, 'ref':str, 'alt':str, 'num_alt_per_site':np.int32})
print('done.', flush=True)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Annotating permutation results (eGenes) ... ', end='', flush=True)
gene_df = pd.read_csv(args.permutation_results, sep='\t', index_col=0)
gene_info = pd.DataFrame(data=[gene_dict[i] for i in gene_df.index], columns=['gene_name', 'gene_chr', 'gene_start', 'gene_end', 'strand'], index=gene_df.index)
gene_df = pd.concat([gene_info, gene_df], axis=1)
assert np.all(gene_df.index==gene_info.index)
gene_df = gene_df.join(snp_lookup_df, on='variant_id')  # add variant information

col_order = ['gene_name', 'gene_chr', 'gene_start', 'gene_end', 'strand',
    'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'variant_id', 'tss_distance']\
    +list(snp_lookup_df.columns)\
    +['minor_allele_samples', 'minor_allele_count', 'maf', 'ref_factor',
    'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta', 'qval', 'pval_nominal_threshold']
gene_df = gene_df[col_order]
snp_lookup_df = None

outname = os.path.split(args.permutation_results)[1].split('.txt.gz')[0]+'.annotated.txt.gz'
with gzip.open(outname, 'wt') as f:
    gene_df.to_csv(f, sep='\t', float_format='%.6g')
print('done.', flush=True)
    
#------------------------------------------------------------------------------
# 2. variant-gene pairs: output new file with all significant pairs
#------------------------------------------------------------------------------
if args.nominal_results:
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Loading variant-gene pairs ... ', end='', flush=True)
    pairs_df = pd.read_csv(args.nominal_results, sep='\t', index_col=1,
        dtype={'gene_id':str, 'variant_id':str, 'tss_distance':np.int32, 'pval_nominal':np.float32, 'slope':np.float32, 'slope_se':np.float32})
    if args.nominal_results_unnormalized:
        # append slope and slope_se in unnormalized units
        raw_pairs_df = pd.read_csv(args.nominal_results_unnormalized[0], sep='\t', usecols=['variant_id', 'slope', 'slope_se'], index_col=0,
            dtype={'variant_id':str, 'slope':np.float32, 'slope_se':np.float32})
        raw_pairs_df.rename(columns={'slope':'slope_'+args.nominal_results_unnormalized[1], 'slope_se':'slope_'+args.nominal_results_unnormalized[1]+'_se'}, inplace=True)
        assert(np.all(pairs_df.index==raw_pairs_df.index))
        pairs_df = pd.concat([pairs_df, raw_pairs_df], axis=1)
        raw_pairs_df = None
    gene_groups = pairs_df.groupby('gene_id', sort=False)
    print('done.', flush=True)

    # apply FDR threshold:
    egene_df = gene_df.loc[gene_df['qval']<=args.fdr, ['pval_nominal_threshold', 'pval_nominal', 'pval_beta']].copy()
    egene_df.rename(columns={'pval_nominal': 'min_pval_nominal'}, inplace=True)

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Parsing significant variant-gene pairs ... ', end='', flush=True)
    signif_pairs_df = []
    for gene_id in egene_df.index:
        igene_df = gene_groups.get_group(gene_id)
        igene_df = igene_df[igene_df['pval_nominal']<=egene_df.loc[gene_id, 'pval_nominal_threshold']]    
        igene_df = igene_df.merge(egene_df.loc[[gene_id]], left_on='gene_id', right_index=True)
        signif_pairs_df.append(igene_df)
    signif_pairs_df = pd.concat(signif_pairs_df, axis=0)    
    outname = os.path.split(args.nominal_results)[1].split('.allpairs.txt.gz')[0]+'.signifpairs.txt.gz'
    with gzip.open(outname, 'wt') as f:
        signif_pairs_df.to_csv(f, sep='\t', float_format='%.6g')    
    print('done.', flush=True)
