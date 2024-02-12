#!/bin/bash


#Run bedtools intersect for each set of features

#TEs:

bedtools intersect -a <(sed '1d' ${1}) \
    -b /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/2024-01-06_RepeatMasker_UCSC_Export.bed \
    /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/2023-12-29_CpGislands_export.bed \
    /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/Hs_EPDnew_006_hg38_900up400down.bed \
    /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/human_imprintome_hg38_ICRs_coordinates.bed \
    -header -C -names TEs CpG_islands Promoters ICRs > intersected/${1}_multi-intersected-counts.bed