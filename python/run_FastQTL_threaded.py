#!/usr/bin/env python3
# Author: Francois Aguet
# Last update: 03/03/2016

import argparse
import os
import numpy as np
import subprocess
import gzip
import multiprocessing as mp


def get_cmd(args, chunk):
    cmd = os.path.join(fastqtl_dir, 'bin', 'fastQTL')+' --vcf '+args.vcf+' --bed '+args.bed+' --window '+args.window\
        +' --maf-threshold '+args.maf_threshold+' --ma-sample-threshold '+args.ma_sample_threshold
    if args.covariates:
        cmd += ' --cov '+args.covariates
    if args.threshold:
        cmd += ' --threshold'+args.threshold
    if args.permute:
        cmd += ' --permute '+' '.join([str(p) for p in args.permute])
    if args.seed:
        cmd += ' --seed '+args.seed
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
parser.add_argument('--chunks', default='100', help='Number of chunks, minimum: #chromosomes')
parser.add_argument('--permute', default=None, type=str, nargs='+', help='Number of permutations, e.g. [1000, 10000] (adaptive). Default: None')
parser.add_argument('--window', default='1e6', help='Cis-window size. Default values is 1Mb (1e6).')
parser.add_argument('--threshold', default='', help='Output only significant phenotype-variant pairs with a p-value below threshold (default 1)')
parser.add_argument('--maf_threshold', default='0.0', help='Include only genotypes with minor allele frequency >=maf_threshold (default 0)')
parser.add_argument('--ma_sample_threshold', default='0', help='Include only genotypes with >=ma_sample_threshold samples carrying the minor allele (default 0)')
parser.add_argument('--fdr', default=0.05, type=np.double)
parser.add_argument('--seed', default=None, help='Random number generator seed')
parser.add_argument('-t', '--threads', default='8', help='Number of threads')
args = parser.parse_args()
fastqtl_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

print('Starting pool', flush=True)
with mp.Pool(processes=int(args.threads)) as pool:
    pdata_res = [pool.map_async(perm_worker, ((args,k),)) for k in np.arange(1,int(args.chunks)+1)]
    pool.close()
    pool.join()
for res in pdata_res:  # check exit status
    assert res.get()[0]==0
print('done', flush=True)

if args.permute:
    # merge chunks
    print('Merging chunks ... ', end='', flush=True)
    cmd = 'zcat '+args.prefix+'_chunk*.txt.gz | gzip -c -1 > '+args.prefix+'.txt.gz && rm '+args.prefix+'_chunk*.txt.gz'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    cmd = 'cat '+args.prefix+'_chunk*.log > '+args.prefix+'.egenes.log && rm '+args.prefix+'_chunk*.log'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    print('done.', flush=True)
    
    # calculate q-values (R script also adds header)
    print('Calculating q-values', flush=True)
    cmd = 'Rscript '+os.path.join(fastqtl_dir, 'R', 'calculateSignificanceFastQTL.R')\
        +' '+args.prefix+'.txt.gz '+str(args.fdr)+' '+args.prefix+'.egenes.txt.gz'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    os.remove(args.prefix+'.txt.gz')
else:
    # merge chunks
    print('Merging chunks ... ', end='', flush=True)
    with gzip.open('header_chunk.txt.gz', 'wb') as f:
        f.write(b'gene_id\tvariant_id\ttss_distance\tpval_nominal\tslope\tslope_se\n')
    cmd = 'zcat header_chunk.txt.gz '+args.prefix+'_chunk*.txt.gz | gzip -c -1 > '+args.prefix+'.txt.gz && rm '+args.prefix+'_chunk*.txt.gz'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    os.remove('header_chunk.txt.gz')
    cmd = 'cat '+args.prefix+'_chunk*.log > '+args.prefix+'.allpairs.log && rm '+args.prefix+'_chunk*.log'
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')
    print('done.', flush=True)
    
    os.rename(args.prefix+'.txt.gz', args.prefix+'.allpairs.txt.gz')
