#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import numpy as np
import pandas as pd
import os
import gzip
from datetime import datetime
import subprocess
import io

parser = argparse.ArgumentParser(description='Filter significant SNP-gene pairs from FastQTL results using FDR cutoff')
parser.add_argument('permutation_results', help='FastQTL output')
parser.add_argument('fdr', type=np.double, help='False discovery rate (e.g., 0.05)')
parser.add_argument('annotation_gtf', help='Annotation in GTF format')
parser.add_argument('--snp_lookup', default='', help='Tab-delimited file with columns: chr, variant_pos, variant_id, ref, alt, num_alt_per_site, rs_id_dbSNP...')
parser.add_argument('--nominal_results', default='', help='FastQTL output from nominal pass')
parser.add_argument('--nominal_results_unnormalized', nargs=2, default='', help='FastQTL output file (nominal pass), and units')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

#------------------------------------------------------------------------------
# 1. eGenes (permutation output): add gene and variant information
#------------------------------------------------------------------------------
gene_dict = {}
print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Parsing GTF', flush=True)
# add: gene_name, gene_chr, gene_start, gene_end, strand
with open(args.annotation_gtf) as gtf:
    for row in gtf:
        row = row.strip().split('\t')
        if row[0][0]=='#' or row[2]!='gene': continue
        # get gene_id and gene_name from attributes
        attr = dict([i.split() for i in row[8].replace('"','').split(';') if i!=''])
        gene_dict[attr['gene_id']] = [attr['gene_name'], row[0], row[3], row[4], row[6]]

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Annotating permutation results (eGenes)', flush=True)
gene_df = pd.read_csv(args.permutation_results, sep='\t', index_col=0)
gene_info = pd.DataFrame(data=[gene_dict[i] for i in gene_df.index], columns=['gene_name', 'gene_chr', 'gene_start', 'gene_end', 'strand'], index=gene_df.index)
gene_df = pd.concat([gene_info, gene_df], axis=1)
assert np.all(gene_df.index==gene_info.index)

col_order = ['gene_name', 'gene_chr', 'gene_start', 'gene_end', 'strand',
    'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df', 'variant_id', 'tss_distance']
if args.snp_lookup:
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Adding variant annotations from lookup table', flush=True)
    # intersect lookup table with variant_ids (col 7 in permutation_results; col 3 in snp_lookup)
    cmd = "awk 'NR==FNR {v[$7]; next} $3 in v' <(zcat "+args.permutation_results+") <(zcat "+args.snp_lookup+")"
    s = subprocess.check_output(cmd, shell=True, executable='/bin/bash')
    snp_lookup_df = pd.read_csv(io.StringIO(s.decode()), index_col=2, sep='\t',
        dtype={'chr':str, 'variant_pos':np.int64, 'variant_id':str, 'ref':str, 'alt':str, 'num_alt_per_site':np.int32})
    gene_df = gene_df.join(snp_lookup_df, on='variant_id')  # add variant information
    col_order += list(snp_lookup_df.columns)
col_order += ['minor_allele_samples', 'minor_allele_count', 'maf', 'ref_factor',
    'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta', 'qval', 'pval_nominal_threshold']
gene_df = gene_df[col_order]

outname = os.path.join(args.output_dir, os.path.split(args.permutation_results)[1].split('.txt.gz')[0]+'.annotated.txt.gz')
with gzip.open(outname, 'wt') as f:
    gene_df.to_csv(f, sep='\t', float_format='%.6g')

#------------------------------------------------------------------------------
# 2. variant-gene pairs: output new file with all significant pairs
#------------------------------------------------------------------------------
if args.nominal_results:
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Filtering significant variant-gene pairs', flush=True)

    # eGenes (apply FDR threshold)
    egene_df = gene_df.loc[gene_df['qval']<=args.fdr, ['pval_nominal_threshold', 'pval_nominal', 'pval_beta']].copy()
    egene_df.rename(columns={'pval_nominal': 'min_pval_nominal'}, inplace=True)
    egene_ids = set(egene_df.index)
    threshold_dict = egene_df['pval_nominal_threshold'].to_dict()

    # process by chunks to reduce memory usage
    signif_df = []
    mask = []
    for i,chunk in enumerate(pd.read_csv(args.nominal_results, sep='\t', iterator=True, chunksize=1000000, index_col=1,
        dtype={'gene_id':str, 'variant_id':str, 'tss_distance':np.int32,
            'ma_samples':np.int32, 'ma_count':np.int32, 'maf':np.float32,
            'pval_nominal':np.float64, 'slope':np.float32, 'slope_se':np.float32})):
        chunk = chunk[chunk['gene_id'].isin(egene_ids)]
        m = chunk['pval_nominal']<chunk['gene_id'].apply(lambda x: threshold_dict[x])
        signif_df.append(chunk[m])
        mask.append(m)
        print('Chunks processed: {0:d}'.format(i+1), end='\r', flush=True)
    signif_df = pd.concat(signif_df, axis=0)

    if args.nominal_results_unnormalized:
        signif_raw_df = []
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Filtering significant variant-gene pairs (unnormalized)', flush=True)
        for i,chunk in enumerate(pd.read_csv(args.nominal_results_unnormalized[0], sep='\t', iterator=True, chunksize=1000000,
            usecols=['gene_id', 'variant_id', 'slope', 'slope_se'], index_col=1, dtype={'gene_id':str, 'variant_id':str, 'slope':np.float32, 'slope_se':np.float32})):
            chunk = chunk[chunk['gene_id'].isin(egene_ids)]
            signif_raw_df.append(chunk[mask[i]])
            print('Chunks processed: {0:d}'.format(i+1), end='\r', flush=True)
        signif_raw_df = pd.concat(signif_raw_df, axis=0)
        signif_raw_df.rename(columns={'slope':'slope_'+args.nominal_results_unnormalized[1], 'slope_se':'slope_'+args.nominal_results_unnormalized[1]+'_se'}, inplace=True)
        assert np.all(signif_df.index==signif_raw_df.index)
        signif_df = pd.concat([signif_df, signif_raw_df[signif_raw_df.columns[1:]]], axis=1)

    signif_df = signif_df.merge(egene_df, left_on='gene_id', right_index=True)

    outname = os.path.join(args.output_dir, os.path.split(args.nominal_results)[1].split('.allpairs.txt.gz')[0]+'.signifpairs.txt.gz')
    with gzip.open(outname, 'wt') as f:
        signif_df.to_csv(f, sep='\t', float_format='%.6g')

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Completed annotation', flush=True)
