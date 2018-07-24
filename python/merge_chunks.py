#!/usr/bin/env python3
# Author: Francois Aguet
import argparse
import os
import numpy as np
import subprocess
import gzip
import contextlib
from datetime import datetime
import tempfile


@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


def merge_chunks(chunk_list_file, header, output_dir, prefix):
    """Merge FastQTL output chunks and add header"""
    with tempfile.NamedTemporaryFile(mode='w+b', dir=output_dir) as header_file:
        with gzip.open(header_file, mode='wt') as gz:
            gz.write('\t'.join(header)+'\n')
        header_file.flush()
        cmd = 'xargs zcat {} < {} | gzip -c -1 > {}'.format(
            header_file.name, chunk_list_file, os.path.join(output_dir, prefix+'.txt.gz'))
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')


def merge_logs(log_list_file, output_dir, prefix):
    """Merge FastQTL logs"""
    cmd = 'xargs cat < {} > {}'.format(
        log_list_file, os.path.join(output_dir, prefix+'.log'))
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run FastQTL')
    parser.add_argument('chunk_list', help='List of chunks')
    parser.add_argument('log_list', help='List of chunks')
    parser.add_argument('prefix', help='Prefix for output file name')
    parser.add_argument('--permute', action='store_true')
    parser.add_argument('--fdr', default=0.05, type=np.double)
    parser.add_argument('--qvalue_lambda', default=None, help='lambda parameter for pi0est in qvalue.')
    parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
    args = parser.parse_args()
    fastqtl_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Merging chunks', flush=True)
    if args.permute:
        header = [
            'gene_id', 'num_var', 'beta_shape1', 'beta_shape2', 'true_df', 'pval_true_df',
            'variant_id', 'tss_distance', 'ma_samples', 'ma_count', 'maf', 'ref_factor',
            'pval_nominal', 'slope', 'slope_se', 'pval_perm', 'pval_beta'
        ]
        with open(args.chunk_list) as f:
            with gzip.open(f.readline().strip(), 'rt') as f2:
                ncol = len(f2.readline().strip().split('\t'))
        if ncol==19:  # grp_permute
            header += ['group_id', 'group_size']
            print('  * group permutation output detected')
        prefix = args.prefix
    else:  # nominal
        header = [
            'gene_id',
            'variant_id', 'tss_distance', 'ma_samples', 'ma_count', 'maf',
            'pval_nominal', 'slope', 'slope_se',
        ]
        prefix = args.prefix+'.allpairs'

    with cd(args.output_dir):
        merge_chunks(args.chunk_list, header, args.output_dir, prefix)
        merge_logs(args.log_list, args.output_dir, prefix)

        if args.permute:
            print('Calculating q-values', flush=True)
            cmd = 'Rscript '+os.path.join(fastqtl_dir, 'R', 'calculateSignificanceFastQTL.R')\
                +' '+args.prefix+'.txt.gz '+str(args.fdr)+' '+args.prefix+'.genes.txt.gz'
            if args.qvalue_lambda is not None:
                cmd += ' --lambda '+args.qvalue_lambda
            subprocess.check_call(cmd, shell=True, executable='/bin/bash')
            os.remove(args.prefix+'.txt.gz')
            os.rename(args.prefix+'.log', args.prefix+'.genes.log')

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Done', flush=True)
