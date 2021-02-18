#!/bin/bash
MXIN=$1
MXOUT=`echo $MXIN | sed 's/.mx/_gnames.mx/'`
GNAMES=../refgenome/Mus_musculus.GRCm38.85.gnames.txt

cut -f1,7- $MXIN | sed 1d | head -1 > header.txt

cut -f1,7- $MXIN | sed 1,2d | sort -k 1b,1 \
| join -t $'\t' -1 1 -2 1 $GNAMES - \
| sed 's/\t/_/' | cat header.txt - > $MXOUT
