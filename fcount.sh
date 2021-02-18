#!/bin/bash
GTF=/data/sequence/enzo_cardiomyocyte_maturation/refgenome/Mus_musculus.GRCm38.85.gtf
OUT=kevinwatt_YAP.mx
featureCounts -Q 10 -T 20 -a $GTF -o $OUT *bam

MXIN=$OUT
MXOUT=`echo $MXIN | sed 's/.mx/_gnames.mx/'`
GNAMES=/data/sequence/enzo_cardiomyocyte_maturation/refgenome/Mus_musculus.GRCm38.85.gnames.txt

cut -f1,7- $MXIN | sed 1d | head -1 > header.txt

cut -f1,7- $MXIN | sed 1,2d | sort -k 1b,1 \
| join -t $'\t' -1 1 -2 1 $GNAMES - \
| sed 's/\t/_/' | cat header.txt - > $MXOUT
