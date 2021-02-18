#!/bin/bash
DIR=/data/sequence/enzo_cardiomyocyte_maturation/refgenome
GTF=/data/sequence/enzo_cardiomyocyte_maturation/refgenome/Mus_musculus.GRCm38.85.gtf

for FQZ in `ls *fastq.gz` ; do

FQ=`echo $FQZ | sed 's/.gz//'`

pigz -dc $FQZ | fastq_quality_trimmer -t 20 -l 20 -Q33 > $FQ

STAR genomeLoad LoadAndKeep --genomeDir $DIR --readFilesIn $FQ --runThreadN 30 \
--sjdbGTFfile $GTF --outSAMattributes NH HI NM MD

rm $FQ
mv Aligned.out.sam ${FQ}.STAR.sam
mv Log.final.out ${FQ}_starlog.txt

( samtools view -uSh ${FQ}.STAR.sam | samtools sort - ${FQ}.STAR
rm ${FQ}.STAR.sam
samtools index ${FQ}.STAR.bam
samtools flagstat ${FQ}.STAR.bam > ${FQ}.STAR.bam.stats ) &

done
STAR genomeLoad Remove --genomeDir $DIR
wait

