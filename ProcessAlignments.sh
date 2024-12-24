#!/bin/bash

WD="/home/oscar/Desktop/BioinformaticProjects/UpDownRegulated/Project/data"
BAMS=$WD"/Alignment"
THREAD="4"
SORTED_BAMS=$WD"/""SortedBAMS"
OUT=$WD"/html/AlignmentQualities"
AlQualRep=$WD"/data/html/Alignment_Quality_Data"


mkdir $SORTED_BAMS
mkdir $OUT

# Sort and indexing
for i in $(ls $BAMS | grep '.bam' | grep -v '.bai')
do
    O=${i%.bam}
    samtools sort -@ $THREAD -o $SORTED_BAMS"/"$O"_sorted.bam" $BAMS"/"$i 
    samtools index -@ $THREAD $SORTED_BAMS"/"$O"_sorted.bam" $SORTED_BAMS"/"$O"_sorted.bam.bai"
    samtools stat --threads $THREAD $SORTED_BAMS"/"$O"_sorted.bam" > $OUT"/"$O"_quality.txt"
done

multiqc $OUT -o $AlQualRep
