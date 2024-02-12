# Load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(bsseq))
suppressPackageStartupMessages(library(dmrseq))
suppressPackageStartupMessages(library(DSS))
suppressPackageStartupMessages(library(argparse))

# Make sure metadata has the following structure:
#    Files   | Group
# -----------------------
# <file_name>| <int 1 or 2>

args = commandArgs(trailingOnly=TRUE)

# Read in data
metadata <- read.csv(args[1])

# get file names
cov_files <- metadata %>%
  pull(Files)

group1_files <- metadata %>%
  filter(Group == 1) %>%
  pull(Files)

group2_files <- metadata %>%
  filter(Group == 2) %>%
  pull(Files)

#Read the bismark files with dmrseq, then use DSS for the actual analysis

# make BSseq object
bismark_cov <- read.bismark(files = cov_files, rmZeroCov = TRUE)

# trt = vector of condition labels for each sample
trt <- metadata %>%
  pull(Group)
pData(bismark_cov)$Condition <- trt

# filter for loci that are not zero in more than 1 sample in a condition
sample_1 <- which(pData(bismark_cov)$Condition == 1)
cov_temp_1 <- bismark_cov[, sample_1]
loci_1_keep <- which(DelayedMatrixStats::rowCounts(getCoverage(cov_temp_1, type="Cov"),value=0)<2)

sample_2 <- which(pData(bismark_cov)$Condition == 2)
cov_temp_2 <- bismark_cov[, sample_2]
loci_2_keep <- which(DelayedMatrixStats::rowCounts(getCoverage(cov_temp_2, type="Cov"),value=0)<2)

loci.idx <- sort(unique(intersect(loci_1_keep, loci_2_keep)))
bismark_cov_filt <- bismark_cov[loci.idx,]


#Use DSS for the actual analysis
#It seems that the group1 and group2 names are the file names of the files in group1 and group2

dmlTest = DMLtest(bismark_cov_filt, group1=group1_files,
                  group2=group2_files,
                  smoothing=TRUE, ncores=20) #use 20 cores
write.csv(dmlTest, args[2])
message("Finished DML testing, starting DMR calling")
dmrs = callDMR(dmlTest) #use default threshold settings

write.csv(dmrs, args[3])


