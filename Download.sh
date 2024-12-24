#!/bin/bash

outDir=$1

threads=$2

sra_id=$3

orDir=$(pwd)

if [[ -z $sra_id ]]; then
    echo "ERROR: SRA ID missing."
    exit 1
fi

if [[ -z $outDir ]]; then
    echo "WARNING: If not specified the output directory will be created in: "$(pwd)
    outDir=$orDir
fi

if [[ -d $outDir"/"$sra_id ]]; then
    echo "ERROR: "$outDir"/"$sra_id "already exists."
    exit 1
fi

if [[ -z $threads ]]; then
    echo "WARNING: If not specified the number of threads will be set to 1."
    threads=1
fi

if [[ -d $outDir ]]; then
    echo "WARNING: Output directory already exists."
else 
    #Create Output directory
    mkdir $outDir
fi

echo "Output directory: " $outDir
#Prefetch
#echo "Changing working directory to: $OutputDirectory"

echo "prefetch -O $outDir $sra_id"
prefetch -O $outDir $sra_id

#fasterqdump
echo "fasterq-dump $sra_id --threads $threads --progress"
cd $outDir
fasterq-dump $sra_id --threads $threads --progress

#fastqc

if [ -f "${sra_id}.fastq" ]; then
    echo "${sra_id}"
    fastqc ${sra_id}.fastq
else
    fastqc ${sra_id}_1.fastq
    fastqc ${sra_id}_2.fastq
fi

#Removing temporal directory created by prefetch
rm -r $sra_id

#gzip ${sra_id}_1.fastq
#gzip ${sra_id}_2.fastq

cd $orDir


