#!/usr/bin/env python3
# Author: Francois Aguet
import argparse
import os
import numpy as np
import subprocess
import gzip
import multiprocessing as mp
import contextlib
from datetime import datetime
import tempfile
import glob

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

def get_cmd(args, chunk):
    cmd = os.path.join(fastqtl_dir, 'bin', 'fastQTL')+' --vcf '+args.vcf+' --bed '+args.bed+' --window '+args.window \
        +' --maf-threshold '+args.maf_threshold \
        +' --ma-sample-threshold '+args.ma_sample_threshold \
        +' --interaction-maf-threshold '+args.interaction_maf_threshold
    if args.covariates:
        cmd += ' --cov '+args.covariates
    if args.phenotype_groups:
        cmd += ' --grp '+args.phenotype_groups
    if args.threshold:
        cmd += ' --threshold '+args.threshold
    if args.permute:
        cmd += ' --permute '+' '.join([str(p) for p in args.permute])
    if args.interaction:
        cmd += ' --interaction '+args.interaction
    if args.best_variant_only:
        cmd += ' --report-best-only'
    if args.seed:
        cmd += ' --seed '+args.seed
    if args.exclude_samples:
        cmd += ' --exclude-samples '+args.exclude_samples
    if args.exclude_sites:
        cmd += ' --exclude-sites '+args.exclude_sites
    cmd += ' --chunk '+str(chunk)+' '+args.chunks\
        + ' --out '+args.prefix+'_chunk{0:03d}.txt.gz'.format(chunk)\
        + ' --log '+args.prefix+'_chunk{0:03d}.log'.format(chunk)
    return cmd

def perm_worker(inputs):
    args = inputs[0]
    chunk = inputs[1]
    cmd = get_cmd(args, chunk)
    print('Processing chunk '+str(chunk), flush=True)
    s = subprocess.check_call(cmd, shell=True, executable='/bin/bash', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print('Finished chunk '+str(chunk), flush=True)
    return s


parser = argparse.ArgumentParser(description='Run FastQTL')
parser.add_argument('vcf', help='Genotypes in VCF 4.1 format')
parser.add_argument('bed', help='Phenotypes in UCSC BED extended format')
parser.add_argument('prefix', help='Prefix for output file name')
parser.add_argument('--covariates', default='', help='Covariates')
parser.add_argument('--phenotype_groups', default='', help='File with mapping of phenotype_id to group_id (gene_id)')
parser.add_argument('--chunks', default='100', help='Number of chunks, minimum: #chromosomes')
parser.add_argument('--permute', default=None, type=str, nargs='+', help='Number of permutations, e.g. [1000, 10000] (adaptive). Default: None (run nominal pass)')
parser.add_argument('--interaction', default=None, type=str, help='Interaction term')
parser.add_argument('--best_variant_only', action='store_true')
parser.add_argument('--window', default='1e6', help='Cis-window size. Default values is 1Mb (1e6).')
parser.add_argument('--threshold', default='', help='Output only significant phenotype-variant pairs with a p-value below threshold (default 1)')
parser.add_argument('--maf_threshold', default='0.0', help='Include only genotypes with minor allele frequency >=maf_threshold (default 0)')
parser.add_argument('--ma_sample_threshold', default='0', help='Include only genotypes with >=ma_sample_threshold samples carrying the minor allele (default 0)')
parser.add_argument('--interaction_maf_threshold', default='0', help='MAF threshold for interactions, applied to lower and upper half of samples')
parser.add_argument('--fdr', default=0.05, type=np.double)
parser.add_argument('--seed', default=None, help='Random number generator seed')
parser.add_argument('--exclude_samples', default=None, help='')
parser.add_argument('--exclude_sites', default=None, help='')
parser.add_argument('--qvalue_lambda', default=None, help='lambda parameter for pi0est in qvalue.')
parser.add_argument('-t', '--threads', default=8, type=int, help='Number of threads')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()
fastqtl_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running FastQTL on {0:d} threads.'.format(args.threads), flush=True)

with cd(args.output_dir):
    with mp.Pool(processes=args.threads) as pool:
        pdata_res = [pool.map_async(perm_worker, ((args,k),)) for k in np.arange(1,int(args.chunks)+1)]
        pool.close()
        pool.join()
    for res in pdata_res:  # check exit status
        assert res.get()[0]==0

    with tempfile.NamedTemporaryFile(mode='w+') as chunk_list_file, \
         tempfile.NamedTemporaryFile(mode='w+') as log_list_file:

        # write chunk and log paths to file
        chunk_files = sorted(glob.glob(args.prefix+'_chunk*.txt.gz'))
        chunk_list_file.write('\n'.join(chunk_files)+'\n')
        chunk_list_file.flush()
        log_files = sorted(glob.glob(args.prefix+'_chunk*.log'))
        log_list_file.write('\n'.join(log_files)+'\n')
        log_list_file.flush()

        # merge chunks
        cmd = 'python3 '+os.path.join(fastqtl_dir, 'python', 'merge_chunks.py') \
            +' {} {} {} --fdr {} -o .'.format(chunk_list_file.name, log_list_file.name, args.prefix, args.fdr)
        if args.qvalue_lambda:
            cmd += ' --qvalue_lambda {}'.format(args.qvalue_lambda)
        if args.permute:
            cmd += ' --permute'
        subprocess.check_call(cmd, shell=True)

        # remove chunk files
        for f in chunk_files + log_files:
            os.remove(f)
