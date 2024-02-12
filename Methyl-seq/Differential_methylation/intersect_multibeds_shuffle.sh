#!/bin/bash


#Run bedtools intersect for each set of features
mkdir -p shuffled
    
for seed in {1..10} ; do
    echo "Running permutation ${seed} on file ${1}"
    bedtools shuffle -noOverlapping -seed $seed -i <(sed '1d' ${1}) -g /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/DMR/Output/hg38.autosomes | \
    bedtools intersect -a stdin -b /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/2024-01-06_RepeatMasker_UCSC_Export.bed \
    /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/2023-12-29_CpGislands_export.bed \
    /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/Hs_EPDnew_006_hg38_900up400down.bed \
    /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/human_imprintome_hg38_ICRs_coordinates.bed \
    -header -C -names TEs CpG_islands Promoters ICRs > shuffled/${1}_multi-intersected-counts_shuffled_${seed}.bed
done