#!usr/bin/env python

import sys
import os
import pandas as pd
import numpy as np
from bioinfokit import visuz
import matplotlib.pyplot as plt
import requests
from requests.structures import CaseInsensitiveDict

#Read the input data
fc_data = pd.read_csv(sys.argv[1])
fc_data.dropna(subset = ['hgnc_symbol'],inplace=True)

#get list of sample names
names = [col[:-5] for col in fc_data.columns if '_padj' in col]

#Set up the output directory for GO enrichment
plot_dir = "GO_enrich/"
os.makedirs(plot_dir,exist_ok=True)

##### RUN GO ENRICHMENT #####

#Thresholds for enrichment
fc_thresh = 1
pval_thresh = 0.05
go_thresh = 0.05

#PANTHERDB (version 17) parameters. Note, the current version is 10.5281/zenodo.6399963
url_base = "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList="
#9606 = human
#GO:0008150 = biological process
#alternatively, GO:0003674 = molecular function
url_tail = "&organism=9606&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR"
#url_tail = "&organism=9606&annotDataSet=GO%3A0003674&enrichmentTestType=FISHER&correction=FDR"

#Request headers
headers = CaseInsensitiveDict()
headers["connection"] = "keep-alive"
headers["keep-alive"] = "timeout=3600, max=100"

for gene in names: #using "gene" as variable name from old code
    fc_col =  gene + "_log2FoldChange" 
    pval_col = gene + "_padj"
    
    degs_up = fc_data[(fc_data[pval_col] < pval_thresh) & (fc_data[fc_col] > fc_thresh)][['hgnc_symbol','ensembl_gene_id_version',fc_col,pval_col]]
    
    degs_up.to_csv(plot_dir+gene + '_up_' + str(fc_thresh) + '.csv')
    
    degs_down = fc_data[(fc_data[pval_col] < pval_thresh) & (fc_data[fc_col] < -fc_thresh)][['hgnc_symbol','ensembl_gene_id_version',fc_col,pval_col]]
    
    degs_down.to_csv(plot_dir+gene + '_down_' + str(fc_thresh) + '.csv')
    
    #Get enrichment data for upregulated genes
    gene_list = degs_up['hgnc_symbol'].tolist()
    
    if len(gene_list) > 1:
        full_url = url_base + ','.join(gene_list) + url_tail

        print("Requesting enrichment data for " + str(len(gene_list)) + " " + gene + " upregulated targets")
        #print(full_url)

        try:
            r = requests.post(full_url, headers = headers)
        except Exception as e:
            print("Failed: " + full_url)
            print(e)

        if r.status_code == 200:
            response_df = pd.DataFrame(r.json()['results']['result'])
            response_df[['id','label']] = pd.json_normalize(response_df['term'])
            response_df.drop(columns = 'term', inplace = True)
            signif = response_df[response_df['fdr'] < go_thresh]
            signif.to_csv(plot_dir+gene + '_' + str(fc_thresh) + "-fold_upregulated_GOs_biological_process.csv")
            #signif.to_csv(gene + '_' + str(fc_thresh) + "-fold_upregulated_GOs_molecular_function.csv")
        else: print("Failed: " + full_url)
        
    #Get enrichment data for downregulated genes
    gene_list = degs_down['hgnc_symbol'].tolist()
    full_url = url_base + ','.join(gene_list) + url_tail
    if len(gene_list) > 1:
        print("Requesting enrichment data for " + str(len(gene_list)) + " " + gene + " downregulated targets")

        try:
            r = requests.post(full_url)
        except:
            print("Failed: " + full_url)

        if r.status_code == 200:
            response_df = pd.DataFrame(r.json()['results']['result'])
            response_df[['id','label']] = pd.json_normalize(response_df['term'])
            response_df.drop(columns = 'term', inplace = True)
            signif = response_df[response_df['fdr'] < go_thresh]
            signif.to_csv(plot_dir+gene + '_' + str(fc_thresh) + "-fold_downregulated_GOs_biological_process.csv")
            #signif.to_csv(gene + '_' + str(fc_thresh) + "-fold_downregulated_GOs_molecular_function.csv")
        else: print("Request failed")

##### MAKE VOLCANO PLOTS #####

#Set up the output directory for volcano plots
plot_dir = "volcano_plots/"
os.makedirs(plot_dir,exist_ok=True)


#also, make volcano plots
for gene in names: #using "gene" as variable name from old code
    fc_col =  gene + "_log2FoldChange" 
    pval_col = gene + "_padj" 


    fc_thresh = 1
    pval_thresh = 10**-17

    degs0 = fc_data[(fc_data[pval_col] < pval_thresh) & (np.abs(fc_data[fc_col]) > fc_thresh)]


    fc_thresh = 4
    pval_thresh = 10**-10

    degs1 = fc_data[(fc_data[pval_col] < pval_thresh) & (np.abs(fc_data[fc_col]) > fc_thresh)]


    fc_thresh = 6
    pval_thresh = 10**-5

    degs2 = fc_data[(fc_data[pval_col] < pval_thresh) & (np.abs(fc_data[fc_col]) > fc_thresh)]


    fc_thresh = 8
    pval_thresh = 0.05

    degs3 = fc_data[(fc_data[pval_col] < pval_thresh) & (np.abs(fc_data[fc_col]) > fc_thresh)]

    gene_list = tuple(degs0['hgnc_symbol']) + tuple(degs1['hgnc_symbol']) + tuple(degs2['hgnc_symbol']) + tuple(degs3['hgnc_symbol'])

    plot_subset = fc_data.dropna(subset = [fc_col, pval_col])
    print(gene)
    try:
        visuz.GeneExpression.volcano(df=plot_subset, geneid = "hgnc_symbol", lfc=fc_col, pv=pval_col,
                                     dotsize = 3, valpha = 0.5,
                                     color = ("blue", "gray", "red"), xlm = (np.floor(plot_subset[fc_col].min()),np.ceil(plot_subset[fc_col].max()),1), ylm = (0,np.ceil(-np.log10(plot_subset[fc_col].min())),5), 
                                     sign_line = True, gfont = 5,
                                     genenames = gene_list, gstyle = 1, figname = plot_dir+'volcano_'+gene)#, show = True)
    except ValueError as e:
        print("overflow")
        visuz.GeneExpression.volcano(df=plot_subset, geneid = "hgnc_symbol", lfc=fc_col, pv=pval_col,
                                     dotsize = 3, valpha = 0.5,
                                     color = ("blue", "gray", "red"), xlm = (np.floor(plot_subset[fc_col].min()),np.ceil(plot_subset[fc_col].max()),1), ylm = (0,100,5), 
                                     sign_line = True, gfont = 5,
                                     genenames = gene_list, gstyle = 1, figname = plot_dir+'volcano_'+gene)#, show = True)

    except AssertionError as e:
        print(e)
        