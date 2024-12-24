GENOME="/home/oscar/Desktop/BioinformaticProjects/UpDownRegulated/Project/data/genome"
FASTA=$GENOME"/GCF_000001635.27_GRCm39_genomic.fna"
GTF=$GENOME"/genomic.gtf"


INDEX="/home/oscar/Desktop/BioinformaticProjects/UpDownRegulated/Project/data/genome/grcm38/genome"
SPLICE_SITES=$GENOME"/splice_sites.ss"
EXONS=$GENOME"/exons.exon"

DO_INDEX="FALSE"

WD="/home/oscar/Desktop/BioinformaticProjects/UpDownRegulated/Project"
FASTQ=$WD"/data/trimmed_data"
OUT=$WD"/data/Alignment"

THREADS="4"

HISAT2_PATH="/home/oscar/Desktop/BioinformaticProjects/UpDownRegulated/Project/HISAT2/hisat2-2.2.1"

mkdir $OUT

if [ "$DO_INDEX" = "TRUE" ]; then

    # Creating Splice Site and Exon files.

    python $HISAT2_PATH"/"extract_splice_sites.py $GTF > $SPLICE_SITES
    python $HISAT2_PATH"/"extract_exons.py $GTF > $EXONS

    # Index
    
    $HISAT2_PATH"/"hisat2-build -p $THREADS --exon $EXONS $FASTA $INDEX
fi

for i in $(ls $FASTQ)
do
    O=${i%.fastq.gz}
    $HISAT2_PATH"/"hisat2 --phred33 --dta -x $INDEX -U $FASTQ"/"$i -S $OUT"/"$O".sam" -p $THREADS --summary-file $OUT"/"$O.txt #--known-splicesite-infile $SPLICE_SITES
    samtools view -bS $OUT"/"$O".sam" > $OUT"/"$O".bam"
    rm $OUT"/"$O".sam"
done