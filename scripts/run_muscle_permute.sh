#!/bin/bash

DIR=/home/jjzhu/source_code/fastqtl-1
PROG=$DIR/bin/fastQTL

# input files
VCF=/scratch/PI/sabatti/controlled_access_data/GTEx_eQTL_analysis/data_copy/phg000520.v2.GTEx_MidPoint_Imputation.genotype-calls-vcf.c1/GTEx_Analysis_20150112_OMNI_2.5M_5M_450Indiv_chr1to22_genot_imput_info04_maf01_HWEp1E6_ConstrVarIDs.vep.vcf.gz
BED=/scratch/PI/sabatti/controlled_access_data/our_analysis/muscle_data/Muscle_Skeletal_Analysis.v6p.FOR_QC_ONLY.normalized.expr.bed.gz
COV=/scratch/PI/sabatti/controlled_access_data/our_analysis/muscle_data/Muscle_Skeletal_Analysis.covariates.txt

DDIR=/scratch/PI/sabatti/controlled_access_data/our_analysis/muscle_fastqtl_out
mkdir -p $DDIR

$PROG  --vcf $VCF  --bed $BED --cov $COV \
       --region 22:17000000-18000000 \
       --window 1e6 \
       --ma-sample-threshold 10 \
       --maf-threshold 0.01 \
       --permute 1000 \
       --out $DDIR/fastqtl.permuted.txt.gz 
