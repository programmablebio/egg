#!/bin/bash
#SBATCH -c 20                              # Request 1 core (max is 20)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=250G                          # 40GB RAM
#SBATCH -o logs/DSS_2024-01-21_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/DSS_2024-01-21_%j.err                 #
#SBATCH --mail-type=ALL

#Arguments:
#$1: input samplesheet (.csv)

source activate EMseq-R-DSS
echo "Starting: ${1}"
Rscript DSS-batch.R $1 ${1%.csv}_diff_methyl_loci.csv ${1%.csv}_diff_methyl_regions.csv
echo "Done calling regions. Zipping output."
gzip ${1%.csv}_diff_methyl_loci.csv
