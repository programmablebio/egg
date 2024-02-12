#!/bin/bash
#SBATCH -c 20                              # Request 20 cores (max is 20)
#SBATCH -t 1-12:00                         # Runtime in D-HH:MM format
#SBATCH -p medium                           # Partition to run in
#SBATCH --mem=128G                          # 40GB RAM
#SBATCH -o /n/groups/church/oogenesis/EMseq-NF-logs/2024-01-02-meth_sample01_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e /n/groups/church/oogenesis/EMseq-NF-logs/2024-01-02-meth_sample01_%j.err                 #
#SBATCH --mail-type=ALL

source activate EMseq
nextflow run nf-core/methylseq \
--input /n/groups/church/oogenesis/EMseq-NF-samplesheets/sample01.csv \
--outdir 2024-01-02_EMseq-out_01 \
--fasta /n/groups/church/oogenesis/EMseq/EMseq-data/genome/grch38_core+bs_controls.fa \
-w /n/scratch/users/m/mdp817/NEXTFLOW01 \
-ansi-log false
