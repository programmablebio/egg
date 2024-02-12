library(dplyr)
library(janitor)
library(DESeq2)

# read the raw counts data
o = read.delim("input_data/2024-02-12_raw_counts_combined.csv",sep=',')


#set the row names
o[is.na(o)] = 0
rownames(o) = o$Geneid
o$X = NULL
o <- select(o, -Geneid)

## Deseq for Differential Genes

# Make metadata dataframe
coldata_o = read.delim("input_metadata/2024-02-12_metadata_D5control_hPSCcombined.csv",sep=',')

#convert all non alphanumeric characters to underscores to avoid problems with R
#note, need to make sure the sample information doesn't get lost here
coldata_o$label <- gsub("[^[:alnum:]]", "_", coldata_o$label)

# Convert the 'label' column to a factor
coldata_o$label <- factor(coldata_o$label)
coldata_o$label <- relevel(coldata_o$label, ref = "Control")
ddso = DESeqDataSetFromMatrix(countData = o,colData = coldata_o, design = ~label + Study) #include Study in the design to deal with batch effects
ddso = DESeq(ddso)

# Get the unique levels of the 'label' column
label_levels <- levels(coldata_o$label)

# Iterate over each unique label level
for (i in label_levels) {
  # Check if the current level is 'Control'
  if (i == "Control") {
    # Get the default results
    res <- results(ddso)
  } else {
    # Get the results comparing the current level to 'Control'
    res <- results(ddso, contrast = c("label", i, "Control"))
  }
  # Create the file name
  file_name <- paste0("DEGs_vs_D5_batch_adjusted/",i, "_vs_D5_batch_adjusted.csv")

  # Write the results to a CSV file
  write.csv(as.data.frame(res), file = file_name)
  print("Wrote:")
  print(file_name)
}

#Export counts data for heatmap
countsData <- counts(ddso, normalized=TRUE)
write.csv(countsData, file="input_data/2024-02-12_normalized_counts_combined.csv")
