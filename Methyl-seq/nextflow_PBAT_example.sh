#!/bin/bash
#SBATCH -c 20                              # Request 20 cores (max is 20)
#SBATCH -t 2-00:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=128G                          # 40GB RAM
#SBATCH -o EMseq-NF-logs_2024-01-11/EMseq_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e EMseq-NF-logs_2024-01-11/EMseq_%j.err                 #
#SBATCH --mail-type=ALL

#Parameters:
#$1: input csv file

source activate EMseq
echo "Processing ${1}"
nextflow run nf-core/methylseq \
--input ${1} \
--outdir EMseq-PBAT-out_${1%.csv} \
--fasta /n/groups/church/oogenesis/EMseq/EMseq-data/genome/grch38_core+bs_controls.fa \
--pbat \
-w /n/scratch/users/m/mdp817/NEXTFLOW/${1%.csv} \
-ansi-log false
