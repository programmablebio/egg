#!/bin/bash


#Run bedtools intersect for each set of features

mkdir -p intersected

#TEs:

bedtools intersect -a <(sed '1d' ${1}) \
    -b /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/2024-01-06_RepeatMasker_UCSC_Export.bed \
    -wao \
    -names TEs > intersected/${1%.bed}_intersect_TEs.bed
echo "Wrote: ${1%.bed}_intersect_TEs.bed" 
    
#CpG islands:

bedtools intersect -a <(sed '1d' ${1}) \
    -b /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/2023-12-29_CpGislands_export.bed \
    -wao \
    -names CpG_islands > intersected/${1%.bed}_intersect_CpG_islands.bed
echo "Wrote: ${1%.bed}_intersect_CpG_islands.bed" 


#Promoters:

bedtools intersect -a <(sed '1d' ${1}) \
    -b /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/Hs_EPDnew_006_hg38_900up400down.bed \
    -wao \
    -names promoters > intersected/${1%.bed}_intersect_promoters.bed
echo "Wrote: ${1%.bed}_intersect_promoters.bed" 

#ICRs

bedtools intersect -a <(sed '1d' ${1}) \
    -b /Users/merrickpiersonsmela/Documents/Church_Lab/oogenesis/Computational/EMseq/Bed_plots/base_beds/human_imprintome_hg38_ICRs_coordinates.bed \
    -wao \
    -names ICRs > intersected/${1%.bed}_intersect_ICRs.bed
echo "Wrote: ${1%.bed}_intersect_ICRs.bed"