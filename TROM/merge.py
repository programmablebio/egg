import pandas as pd
import numpy as np
atlas_labels = pd.read_csv('atlas/atlas_count_matrix.csv')
atlas_labels.set_index('Gene',inplace=True)
pgclc_atlas = pd.read_csv('atlas/pgclc_countmatrix_garyk.csv')
pgclc_atlas.set_index('Gene',inplace=True)

atlas_labels.rename(columns = {"Zhang2017_Primay_oocyte_3": "Zhang2017_Primary_oocyte_3"},
          inplace = True)

# pgclc_atlas.drop(['iPSC_1','iPSC_2'],axis=1,inplace=True)

samples = ['iPSC_1','iPSC_2','iMeLC_1','iMeLC_2','d6PGCLC_1','d6PGCLC_2','d6PGCLC_3','hPGCw7f_1','hPGCw7f_2','hPGCw9f_1','hPGCw9f_2','ag7_1','ag7_2','ag21_1','ag21_2','ag35_1','ag35_2','ag77','ag120','FetOva','Primordial_oocyte_2','Primordial_oocyte_3','Primary_oocyte_1','Primary_oocyte_2','Primary_oocyte_3','Secondary_oocyte_1','Secondary_oocyte_2','Secondary_oocyte_3','Antral_oocyte_1','Antral_oocyte_2']

columns = list(atlas_labels.columns)
valid_columns = []

for sample in samples:
    for col in columns:
        if sample in col:
            print(col)
            valid_columns.append(col)

valid_columns.extend(['D3-DDX4', 'DNR3'])
valid_columns.remove('iPSC_13D7')


new_df = atlas_labels[valid_columns].join(pgclc_atlas[pgclc_atlas.columns[~pgclc_atlas.columns.isin(valid_columns)]])


new_df.to_csv('atlas/merged_atlas_samples.csv')

pd.read_csv('atlas/merged_atlas_samples.csv')
.set_index('Gene',inplace=True)


sample_names = list(new_df.columns)
condition = ['Yamashiro2018_IPSC', 'Yamashiro2018_IPSC', 'Yamashiro2018_iMeLC',
       'Yamashiro2018_iMeLC', 'Yamashiro2018_d6PGCLC',
       'Yamashiro2018_d6PGCLC', 'Yamashiro2018_d6PGCLC',
       'Tang2015_hPGCw7f', 'Tang2015_hPGCw7f', 'Tang2015_hPGCw9f',
       'Tang2015_hPGCw9f', 'Yamashiro2018_ag7', 'Yamashiro2018_ag7',
       'Yamashiro2018_ag21', 'Yamashiro2018_ag21', 'Yamashiro2018_ag35',
       'Yamashiro2018_ag35', 'Yamashiro2018_ag77_1390G3_AG+VT-',
       'Yamashiro2018_ag77', 'Yamashiro2018_ag77',
       'Yamashiro2018_ag77_1390G3_VT',
       'Yamashiro2018_ag120_1390G3',
       'Yamashiro2018_ag120_1390G3_VT',
       'Yamashiro2018_ag120_1390G3_VT',
       'Yamashiro2018_ag120_1390G3_VT',
       'Yamashiro2018_ag120_1390G3', 'Yatsenko2019_FetOva',
       'Yatsenko2019_FetOva', 'Yatsenko2019_FetOva',
       'Zhang2017_Primordial_oocyte', 'Zhang2017_Primordial_oocyte',
       'Zhang2017_Primary_oocyte', 'Zhang2017_Primary_oocyte',
       'Zhang2017_Primary_oocyte', 'Zhang2017_Secondary_oocyte',
       'Zhang2017_Secondary_oocyte', 'Zhang2017_Secondary_oocyte',
       'Zhang2017_Antral_oocyte', 'Zhang2017_Antral_oocyte', 'D3-DDX4',
       'DNR3', 'Control-N3V', 'DLX5-N3V', 'HHEX-N3V', 'FIGLA-N3V',
       'Control-EB', 'Control', 'Control', 'F2_LTC', 'F3_LTC']
sampl_dict = {'SampleName': sample_names, 'Condition': condition}
sample_sheet = pd.DataFrame(sampl_dict)
sample_sheet.set_index('SampleName',inplace=True)
sample_sheet.to_csv("atlas/samplesheet_merged.csv")

for col in columns:
    if 'oocyte' in col.lower():
        print(col)



import pandas as pd
import os
import sys
import pyensembl

#ensembl = pyensembl.EnsemblRelease(release=75)
file_path = "results/"
csv_files = sorted(os.listdir(file_path))

log2fc_df = pd.DataFrame()

ensembl_ids = pd.read_csv(file_path+csv_files[0]).iloc[:,0].tolist()

log2fc_df['Gene'] = ensembl_ids

for f in csv_files:
    sample_df = pd.read_csv(file_path+f)
    log2fc_df[f.split(".")[0]] = sample_df['log2FoldChange'].tolist()

log2fc_df.set_index('Gene',inplace=True)
log2fc_df.dropna(inplace = True)



log2fc_df.to_csv('merged_atlas_log2fc.csv', index = True)
