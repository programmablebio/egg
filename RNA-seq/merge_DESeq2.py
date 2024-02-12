#!usr/bin/env python

import pandas as pd
import os
import sys

#Merge the files in the folder provided as an argument
folder = sys.argv[1]

#Read gene symbols (this is a hardcoded .csv)
symbols = pd.read_csv("input_metadata/2024-02-06_ensembl_id_to_hgnc.csv")
symbols = symbols.drop_duplicates(subset='ensembl_gene_id_version', keep="first")

#Read and merge the data
data = {}

for file in os.listdir(folder):
    if file.endswith('.csv'):
        print("Reading: "+file)
        data[file[:-4]] = pd.read_csv(os.path.join(folder,file), index_col = 0)

merged_df = pd.DataFrame()

for key, value in data.items():
    merged_df["baseMean"] = value["baseMean"]
    merged_df[key+"_log2FoldChange"] = value["log2FoldChange"]
    merged_df[key+"_padj"] = value["padj"]

#Merge gene symbols and write output
merged_df = merged_df.merge(symbols, how = 'left', right_on = 'ensembl_gene_id_version', left_index = True)
merged_df.set_index("ensembl_gene_id_version", inplace = True, verify_integrity=True)
output_name = os.path.join(folder,"merged_DESeq2.csv")
merged_df.to_csv(output_name)
print("Wrote: "+output_name)

