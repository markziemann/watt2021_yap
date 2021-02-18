#!/bin/bash
REF=/data/projects/mziemann/rna-seq-qc/mouse_rRNA.fa
for FQZ in *fastq.gz ; do
echo $FQZ
FQ=`echo $FQZ | sed 's/.gz//'`
zcat $FQZ | head -4000000 | fastq_quality_trimmer -t 30 -l 20 -Q33 \
| tee $FQ | bwa aln -t 30 $REF - | bwa samse $REF - $FQ \
| samtools view -uSh - \
| samtools sort - ${FQ}.sort &
done
wait
parallel samtools index ::: *bam
for i in *bam ; do samtools flagstat $i > ${i}.stats & done ; wait
