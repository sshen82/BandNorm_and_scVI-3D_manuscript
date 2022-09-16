#!/bin/bash


START=$(date +%s)
echo $START

python3.7  scVI-3D.py -b "${bandMax}" -c "${chrom}" -r "${resolution}" -i "${inPath}" -o "${outPath}" -cs "${cellSummary}" -g "${genomeSize}" -pool "${strategy}" -br -n 100 -gpu -p "${cpuN}" -pca 15 -up -tp -v


END=$(date +%s)
DIFF=$(( $END - $START ))
echo "It took $DIFF seconds"
