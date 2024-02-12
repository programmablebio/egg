#!/bin/bash
#SBATCH -c 5                              # Request 1 core (max is 20)
#SBATCH -t 0-04:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=120G                          # 40GB RAM
#SBATCH -o logs/CpG_multiaggregate_cov_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e logs/CpG_multiaggregate_cov_%j.err                 #
#SBATCH --mail-type=END

#Arguments:
#    $1: .cov.gz file with methylation data
source activate EMseq-postprocess
echo "Processing: aggregating ${1} by all .bed files specified"

#Make a temporary uncompressed .bed file
zcat $1 | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$5","$5+$6}' > "${1%.cov.gz}.bed"

echo "Finished making temporary .bed file"

#Run the chromosome aggregation
#Check if output file exists:
if ! [ -f "${1%.cov.gz}.chrom" ]; then
    echo "Aggregating chromosomes"
    CpG_aggregation.py -i "${1%.cov.gz}.bed" -b /n/groups/church/oogenesis/EMseq/EMseq-data/annotations/hg38+controls.chrom.bed -o "${1%.cov.gz}.chrom" -t count &
fi

#Run the CpG island aggregation
if ! [ -f "${1%.cov.gz}.CpG_islands_aggregated.bed" ]; then
    echo "Aggregating CpG islands"
    CpG_aggregation.py -i "${1%.cov.gz}.bed" -b /n/groups/church/oogenesis/EMseq/EMseq-data/annotations/2023-12-29_CpGislands_export.bed -o "${1%.cov.gz}.CpG_islands_aggregated.bed" -t count &
fi

#Run the promoter aggregation
if ! [ -f "${1%.cov.gz}.prom900up400down_aggregated.bed" ]; then
    echo "Aggregating promoters"
    CpG_aggregation.py -i "${1%.cov.gz}.bed" -b /n/groups/church/oogenesis/EMseq/EMseq-data/annotations/Hs_EPDnew_006_hg38_900up400down.bed -o "${1%.cov.gz}.prom900up400down_aggregated.bed" -t count &
fi

#Run the imprinting aggregation
if ! [ -f "${1%.cov.gz}.imprints_aggregated.bed" ]; then
    echo "Aggregating imprints"
    CpG_aggregation.py -i "${1%.cov.gz}.bed" -b /n/groups/church/oogenesis/EMseq/EMseq-data/annotations/human_imprintome_hg38_ICRs_coordinates.bed -o "${1%.cov.gz}.imprints_aggregated.bed" -t count &
fi

#Run the transposable element aggregation
if ! [ -f "${1%.cov.gz}.TEs_aggregated.bed" ]; then
    echo "Aggregating TEs"
    CpG_aggregation.py -i "${1%.cov.gz}.bed" -b /n/groups/church/oogenesis/EMseq/EMseq-data/annotations/2024-01-06_RepeatMasker_UCSC_Export.bed -o "${1%.cov.gz}.TEs_aggregated.bed" -t count
fi

wait

echo "Done aggregating. Removing temporary .bed file"
rm "${1%.cov.gz}.bed"
