#!/usr/bin/env python3
# Author: Francois Aguet
import argparse
import os
import subprocess
import contextlib
from datetime import datetime

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

def get_cmd(args):
    cmd = os.path.join(fastqtl_dir, 'bin', 'fastQTL')+' --vcf '+args.vcf+' --bed '+args.bed+' --window '+args.window\
        +' --maf-threshold '+args.maf_threshold+' --ma-sample-threshold '+args.ma_sample_threshold
    if args.covariates:
        cmd += ' --cov '+args.covariates
    if args.phenotype_groups:
        cmd += ' --grp '+args.phenotype_groups
    if args.threshold:
        cmd += ' --threshold '+args.threshold
    if args.permute:
        cmd += ' --permute '+' '.join([str(p) for p in args.permute])
    if args.best_variant_only:
        cmd += ' --report-best-only'
    if args.seed:
        cmd += ' --seed '+args.seed
    if args.exclude_samples:
        cmd += ' --exclude-samples '+args.exclude_samples
    cmd += ' --chunk {} {}'.format(*args.chunk)\
        + ' --out '+args.prefix+'_chunk{0:03d}.txt.gz'.format(args.chunk[0])\
        + ' --log '+args.prefix+'_chunk{0:03d}.log'.format(args.chunk[0])
    return cmd


parser = argparse.ArgumentParser(description='Run FastQTL')
parser.add_argument('vcf', help='Genotypes in VCF 4.1 format')
parser.add_argument('bed', help='Phenotypes in UCSC BED extended format')
parser.add_argument('prefix', help='Prefix for output file name')
parser.add_argument('chunk', type=int, nargs=2, help='Chunk and total number of chunks, e.g., 1 100')
parser.add_argument('--covariates', default='', help='Covariates')
parser.add_argument('--phenotype_groups', default='', help='File with mapping of phenotype_id to group_id (gene_id)')
parser.add_argument('--permute', default=None, type=str, nargs='+', help='Number of permutations, e.g. [1000, 10000] (adaptive). Default: None (run nominal pass)')
parser.add_argument('--best_variant_only', action='store_true')
parser.add_argument('--window', default='1e6', help='Cis-window size. Default values is 1Mb (1e6).')
parser.add_argument('--threshold', default='', help='Output only significant phenotype-variant pairs with a p-value below threshold (default 1)')
parser.add_argument('--maf_threshold', default='0.0', help='Include only genotypes with minor allele frequency >=maf_threshold (default 0)')
parser.add_argument('--ma_sample_threshold', default='0', help='Include only genotypes with >=ma_sample_threshold samples carrying the minor allele (default 0)')
parser.add_argument('--seed', default=None, help='Random number generator seed')
parser.add_argument('--exclude_samples', default=None, help='')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()
fastqtl_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

with cd(args.output_dir):
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Processing chunk {}/{}'.format(*args.chunk), flush=True)
    s = subprocess.check_call(get_cmd(args), shell=True, executable='/bin/bash')
    assert s==0
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Done.', flush=True)
