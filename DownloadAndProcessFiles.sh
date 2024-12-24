#!/bin/bash

METADATA="Metadata.csv"
THREADS=8
WD="/home/oscar/Desktop/BioinformaticProjects/UpDownRegulated/Project"
mkdir $WD"/data"
mkdir $WD"/data/fastq"

ID=$(awk -F ',' '{print $1}' $METADATA | tail -n+2)

for i in $ID
do
	bash Download.sh $WD"/data/fastq" $THREADS $i
done
