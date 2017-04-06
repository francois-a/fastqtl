#!/bin/bash
TDIR=/home/jjzhu/src
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TDIR/boost_1_58_0/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TDIR/cnpy/lib/


DIR=/home/jjzhu/source_code/fastqtl-1
ODIR=$DIR/example/cis_eqtl_mtx
PROG=$DIR/bin/fastQTL
mkdir -p $ODIR

$PROG  --vcf $DIR/example/genotypes.vcf.gz \
               --bed $DIR/example/phenotypes.bed.gz \
               --region 22:17000000-18000000 \
               --window 1e6 \
               --ma-sample-threshold 10 \
               --maf-threshold 0.01 \
               --mtx $ODIR \
               --out $ODIR/fastqtl.nonminal.txt.gz 

# echo ""
# echo ""
# echo "CHECK OUTPUT MATRICES:"
# echo ""
# for I in 0 1 2 3 4
# do
#     FNAME=$DIR/example/cis_eqtl_mtx/$I.txt
#     wc -l $FNAME 
#     head -n3 $FNAME | cut -f1-3
#     awk '{print NF}' $FNAME | sort -nu | tail -n 1
# done
