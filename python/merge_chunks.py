#!/usr/bin/env python3
# Author: Francois Aguet
import argparse
import os
import numpy as np
import subprocess
import gzip
import contextlib
from datetime import datetime

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)

parser = argparse.ArgumentParser(description='Run FastQTL')
parser.add_argument('chunk_list', help='List of chunks')
parser.add_argument('log_list', help='List of chunks')
parser.add_argument('prefix', help='Prefix for output file name')
parser.add_argument('--permute', action='store_true')
parser.add_argument('--fdr', default=0.05, type=np.double)
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()
fastqtl_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

with cd(args.output_dir):
    if args.permute:
        # merge chunks
        print('Merging chunks ... ', end='', flush=True)
        cmd = 'xargs zcat < '+args.chunk_list+' | gzip -c -1 > '+args.prefix+'.txt.gz'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        cmd = 'xargs cat < '+args.log_list+' > '+args.prefix+'.egenes.log'
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
        with gzip.open('header_chunk.txt.gz', 'wb') as f:  # order from analysisNominal.cpp
            f.write(('\t'.join([
                'gene_id',
                'variant_id',
                'tss_distance',
                'ma_samples',
                'ma_count',
                'maf',
                'pval_nominal',
                'slope',
                'slope_se',
            ])+'\n').encode('utf-8'))
        cmd = 'zcat header_chunk.txt.gz <(xargs cat < '+args.chunk_list+') | gzip -c -1 > '+args.prefix+'.txt.gz'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        os.remove('header_chunk.txt.gz')
        cmd = 'xargs cat < '+args.log_list+' > '+args.prefix+'.allpairs.log'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        print('done.', flush=True)

        os.rename(args.prefix+'.txt.gz', args.prefix+'.allpairs.txt.gz')
