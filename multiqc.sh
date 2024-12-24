#!/bin/bash
WD="/home/oscar/Desktop/BioinformaticProjects/UpDownRegulated/Project"
OUT=$WD"/data/html/"

multiqc $WD"/data/html/fastqc/" -o $OUT"/Raw_Quality"