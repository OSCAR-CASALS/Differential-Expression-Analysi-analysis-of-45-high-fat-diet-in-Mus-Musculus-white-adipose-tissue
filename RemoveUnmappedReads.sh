#!/bin/bash

WD="/home/oscar/Desktop/BioinformaticProjects/UpDownRegulated/Project"
BAMS=$WD"/data/SortedBAMS"
OUT=$WD"/data/FilteredReads"

THREAD="4"

STATS=$WD"/data/html/Filtered_Alignment"

REPORT=$WD"/data/html/Filtered_Alignment_Report"

mkdir $OUT

mkdir $STATS

for i in $(ls $BAMS | grep '.bam' | grep -v '.bai')
do
    O=${i%.bam}
    echo "Removing unmaped reads from: " $i

    samtools view -b -F 4 $BAMS"/"$i > $OUT"/"$O"_unmapped_reads_filtered.bam"

    echo "Sorting: " $i

    samtools sort -@ $THREAD -o $OUT"/Temp_"$O"_unmapped_reads_filtered.bam" $OUT"/"$O"_unmapped_reads_filtered.bam"

    mv $OUT"/Temp_"$O"_unmapped_reads_filtered.bam" $OUT"/"$O"_unmapped_reads_filtered.bam"

    echo "Indexing: " $i

    samtools index -@ $THREAD $OUT"/"$O"_unmapped_reads_filtered.bam" $OUT"/"$O"_unmapped_reads_filtered.bam.bai"

    echo "Computing Statistics: " $i

    samtools stat --threads $THREAD $OUT"/"$O"_unmapped_reads_filtered.bam" > $STATS"/"$O"_quality.txt"
    
done

multiqc $STATS -o $REPORT