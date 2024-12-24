# trimming

WD="/home/oscar/Desktop/BioinformaticProjects/UpDownRegulated/Project"

FASTQ=$WD"/data/fastq"
OUT=$WD"/data/trimmed_data"
QUALITY=$WD"/data/html/fastqc_trimmed"
THREADS=8

mkdir $OUT
mkdir $QUALITY

for i in $(ls $FASTQ)
do
    O=${i%.fastq}
    trimmomatic SE -phred33 -threads $THREADS \
    $FASTQ"/"$i \
    $OUT"/"$O"_trimmed.fastq" \
    ILLUMINACLIP:data/custom_adapters.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    HEADCROP:12

    fastqc $OUT"/"$O"_trimmed.fastq" -o $QUALITY
done

multiqc $QUALITY -o "data/html/Trimmed_Quality"