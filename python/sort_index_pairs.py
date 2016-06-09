#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='')
parser.add_argument('fastqtl_pairs_file', help='FastQTL output from nominal pass')
parser.add_argument('--tmp_dir', default='', help='tmp directory for sorting')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

print('Sorting and indexing '+args.fastqtl_pairs_file+' (tabix) ... ', end='', flush=True)
file_base = os.path.split(args.fastqtl_pairs_file)[1].split('.txt.gz')[0]
file_name = os.path.join(args.output_dir, file_base+'.sorted.txt.gz')

cmd = 'cat <(echo -e "#chr\tpos\tvariant_id\tgene_id\ttss_distance\tpval_nominal\tslope\tslope_se") \
    <(zcat testNominal.allpairs.txt.gz | tail -n+2 | awk \'{print $2,$1,$3,$4,$5,$6}\' OFS="\t" | awk -F"_" \'{print $1"\t"$2"\t"$0}\')'
if args.tmp_dir:
    cmd += ' | sort -T '+args.tmp_dir+' -n -k1,1 -k2,2 | bgzip -c > '+file_name
else:
    cmd += ' | sort -n -k1,1 -k2,2 | bgzip -c > '+file_name

subprocess.check_call(cmd, shell=True, executable='/bin/bash')
subprocess.call('tabix -f -b 2 -e 2 '+file_name, shell=True, executable='/bin/bash')
print('done.')
