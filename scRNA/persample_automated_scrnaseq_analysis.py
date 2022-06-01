import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import os
import scipy.io as sio
import sys
import scvelo as scv
import stream as st
import progeny
import matplotlib.pyplot as plt
import phate
import scanpy.external as sce

filepath = "./" + sys.argv[1] + "/"
last_folder = sys.argv[1].split('/')[-1]
sc.settings.verbosity = 3
sc.settings.autoshow = False
figpath = filepath + "target_figures_" + last_folder + "/"
sc.settings.set_figure_params(dpi=300, transparent = True)
sc.settings.figdir= figpath

#### load in target variable csv
if sys.argv[1]=='100kfinal' or sys.argv[1]=='oogonia_day25':
    df = pd.read_csv('100k_cell_kit.csv')
    df = df.fillna('')
    arr = df.values
elif sys.argv[1]=='1milfinal':
    df = pd.read_csv('1milscr.csv')
    df = df.fillna('')
    arr = df.values

sample_targ_dict = {}
for row in arr:
    cells_and_genes = ''

    for i in range(3,len(row)):
        if row[i]!= '':
            if cells_and_genes == '':
                cells_and_genes = cells_and_genes + row[i]

            cells_and_genes = cells_and_genes +', '+ row[i]
            cells_and_genes = cells_and_genes

    tmp = cells_and_genes.split(', ')

    sample_targ_dict[row[1]] = tmp

#--------------------------------------------------------------------------------------------------------------
# Basic Parse Filtering
print("Loading DGE")
adata = sc.read_mtx(filepath + 'DGE.mtx')
# reading in gene and cell data
gene_data = pd.read_csv(filepath + 'all_genes.csv')
cell_meta = pd.read_csv(filepath + 'cell_metadata.csv')

genes_df = pd.read_csv(filepath+"all_genes.csv")
gene_names = genes_df['gene_name'].tolist()
adata.var_names = gene_names

# find genes with nan values
notNa = gene_data[gene_data.gene_name.notnull()].index
notNa = notNa.to_list()

# remove genes with nan values
adata = adata[:,notNa]
adata.var_names_make_unique()

# add cell meta data to anndata object
print("add cell meta data to anndata object")
adata.obs = cell_meta
adata.obs.index = cell_meta.bc_wells
adata.obs_names_make_unique()

samp_dct = {}
batches = set(adata.obs['sample'].tolist())
print("Per batch loop")
for i in batches:
    try:
        ###sample specific targets
        target_gene_cell_list = sample_targ_dict[i]
        print(target_gene_cell_list)
        sample_adata = adata[adata.obs['sample'].isin([i])]
        figpath_sample = filepath + "target_figures_" + last_folder + "/"+i+"/"
        if not os.path.exists(figpath_sample):
            os.makedirs(figpath_sample)
        matrix_path_sample = filepath + last_folder + "/matrices/"+i+"/"
        if not os.path.exists(matrix_path_sample):
            os.makedirs(matrix_path_sample)

        sc.pp.filter_cells(sample_adata, min_genes=1000)
        #sc.pp.filter_genes(sample_adata, min_cells=5)

        sample_adata.var['mt'] = sample_adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(sample_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        #
        # Scanpy will prepend the string in the save argument with "violin"
        # and save it to our figure directory defined in the first step.
        # sc.pl.violin(sample_adata, ['n_genes_by_counts'], save='_n_genes', jitter=0.4)
        # sc.pl.violin(sample_adata, ['total_counts'], save='_total_counts', jitter=0.4)
        # sc.pl.violin(sample_adata, ['pct_counts_mt'], save='_mito_pct', jitter=0.4)

        # Filter the data
        sample_adata = sample_adata[sample_adata.obs.n_genes_by_counts < 8000,:]
        sample_adata = sample_adata[sample_adata.obs.total_counts < 50000,:]
        sample_adata = sample_adata[sample_adata.obs.pct_counts_mt < 15,:]

        sc.pl.scatter(sample_adata, x='total_counts', y='n_genes_by_counts', save=i+'_gene_vs_transcript_counts')
        print('median transcript count per cell: ' + str(sample_adata.obs['tscp_count'].median(0)))
        print('median gene count per cell: ' + str(sample_adata.obs['gene_count'].median(0)))

        sc.pp.normalize_total(sample_adata, target_sum=1e4)
        sc.pp.log1p(sample_adata)

        sc.pp.highly_variable_genes(sample_adata, min_mean=0.0125, max_mean=3, min_disp=0.25)
        sc.pl.highly_variable_genes(sample_adata, save='') # scanpy generates the filename automatically
        #sample_adata = sample_adata[:, sample_adata.var.highly_variable]
        #sc.pp.regress_out(sample_adata, ['total_counts'])
        sc.pp.scale(sample_adata, max_value=10)

        sc.tl.pca(sample_adata, svd_solver='arpack')
        # sc.pl.pca_variance_ratio(sample_adata, log=True, n_pcs=50, save='') # scanpy generates the filename automatically
        #--------------------------------------------------------------------------------------------------------------
        # Generate UMAP
        print('umapping')
        sc.pp.neighbors(sample_adata, n_neighbors=10, n_pcs=30)
        sc.tl.umap(sample_adata)
        sc.tl.leiden(sample_adata, resolution=0.5)
        sc.pl.umap(sample_adata, color=['leiden'], legend_fontsize=8, save=i+'_leiden')

        #--------------------------------------------------------------------------------------------------------------

        # Generate Cell-Type Score Matrix

        score_df = pd.DataFrame()
        score_df['Cell'] = (sample_adata.obs.bc_wells).tolist()
        score_df['Number of Genes Used'] = (sample_adata.obs.n_genes).tolist()
        score_df['Leiden Cluster'] = (sample_adata.obs.leiden).tolist()

        marker_df = pd.read_csv('marker_lists.csv')
        for cell in (marker_df.columns.tolist()):
            marker_list_for_cell = marker_df[cell]
            marker_list_for_scores = []
            for marker_gene in marker_list_for_cell:
                if marker_gene in sample_adata.var_names:
                    marker_list_for_scores.append(marker_gene)
            if(len(marker_list_for_scores) > 0):
                sc.tl.score_genes(sample_adata, gene_list=marker_list_for_scores)
                score_df[cell] = (sample_adata.obs.score).tolist()
                # print("Finished with " + cell)
        score_df.to_csv(matrix_path_sample+i+'_cellscorematrix.csv', index = False)

        #--------------------------------------------------------------------------------------------------------------
        # Generate Filtered Data Matrix

        sample_adata_genenames = []
        for i in (sample_adata.var_names).tolist():
            sample_adata_genenames.append(i)
        filtered_data = pd.DataFrame(data=sample_adata.X.transpose(), index=sample_adata.var_names, columns=sample_adata.obs_names)
        filtered_data.insert(0, 'Gene', sample_adata_genenames)
        filtered_data.to_csv(matrix_path_sample+i+'_filtereddata.csv', index = False)
        #--------------------------------------------------------------------------------------------------------------
        # Generate UMAPs for Specified Cell Types
        # score_df = pd.read_csv(figpath+'cellscorematrix.csv')
        score_matrix_df = score_df
        print('loaded and running cell specific umaps')
        if target_gene_cell_list != ['']:
            for j in range(len(target_gene_cell_list)):
                cell_type = target_gene_cell_list[j]
                if(cell_type not in score_matrix_df.columns):
                    try:
                        sc.pl.umap(sample_adata, color=[cell_type], legend_fontsize=8, save=i+'_' + cell_type)
                    except:
                        print('Not found: ' +cell_type)
                        continue
                else:
                    sample_adata.obs[cell_type] = score_matrix_df[cell_type].tolist()
                    sc.pl.umap(sample_adata, color=[cell_type], legend_fontsize=8, save=i+'_' + cell_type)

        #--------------------------------------------------------------------------------------------------------------
        # PHATE for leiden and cell types
        print('PHATE')
        sample_adata.obs['leiden cluster'] = score_matrix_df['Leiden Cluster'].tolist()
        sce.tl.phate(sample_adata, k=5, a=20, t=150)
        sce.pl.phate(
            sample_adata,
            color='leiden',
            color_map='tab20',
            save=i+'_leiden'
        )

        if target_gene_cell_list != ['']:
            for j in range(len(target_gene_cell_list)):
                cell_type = target_gene_cell_list[j]
                try:
                    sce.pl.phate(
                        sample_adata,
                        color=cell_type,
                        color_map='tab20',
                        save=i+'_'+cell_type
                    )
                except:
                    continue

        #--------------------------------------------------------------------------------------------------------------
        # Generate trajectories using stream
        print('Stream trajectories')

        st.select_variable_genes(sample_adata, fig_path=figpath) # stream select variable genes
        try:
            st.dimension_reduction(sample_adata,method='umap',feature='var_genes',n_components=10,nb_pct=0.025) # calculate dim reduce

            sample_adata_low = st.switch_to_low_dimension(sample_adata,n_components=2) # switch to 2d for low dim

            st.seed_elastic_principal_graph(sample_adata_low)

            init_nodes_pos,init_edges = st.infer_initial_structure(sample_adata_low)

            st.seed_elastic_principal_graph(sample_adata,init_nodes_pos=init_nodes_pos,init_edges=init_edges)

            st.elastic_principal_graph(sample_adata,incr_n_nodes=10,fig_path=figpath)
            st.extend_elastic_principal_graph(sample_adata)

            sample_adata.obs['leiden cluster'] = score_matrix_df['Leiden Cluster'].tolist() # leiden cluster to sample_adata, for coloring plots
            print('Stream plotting per cell')
            ### generate STREAM plot for cell types in args
            if target_gene_cell_list != ['']:
                for j in range(len(target_gene_cell_list)):
                    cell = target_gene_cell_list[j]
                    # if(cell in score_matrix_df.columns):
                    #     print(cell + " is not present in the current cell atlas.")
                    # else:
                    try:
                        st.plot_visualization_2D(sample_adata,color=[cell],save_fig=True,fig_path=figpath,fig_ncol=4)
                        st.plot_dimension_reduction(sample_adata,color=[cell],n_components=3,show_graph=False,show_text=True,save_fig=False,fig_path = figpath, fig_name=cell+'_dimreduction.pdf')
                        st.plot_stream_sc(sample_adata,root='S0',color=[cell],dist_scale=0.5,save_fig=True,fig_path=figpath_sample)
                        st.plot_flat_tree(sample_adata,color=[cell],dist_scale=0.3,save_fig=True,fig_path=figpath,fig_name=i+'_'+cell+'_flattree.pdf',show_graph=True,show_text=True)
                    except:
                        print('Not found: ' +cell)
                        continue

            try:
                st.plot_visualization_2D(sample_adata,color=['leiden cluster'],save_fig=True,fig_path=figpath,fig_ncol=4)
                st.plot_dimension_reduction(sample_adata,color=['leiden cluster'],n_components=3,show_graph=False,show_text=True,save_fig=False,fig_path = figpath, fig_name='leiden_dimreduction.pdf')
                st.plot_stream_sc(sample_adata,root='S0',color=['leiden cluster'],dist_scale=0.5,save_fig=True,fig_path=figpath)
                # st.plot_flat_tree(sample_adata, color = ['leiden cluster'], dist_scale=0.3,show_graph=False,show_text=False,save_fig=True,fig_path = figpath,fig_name='leiden_flattree.pdf')
            except:
                print('Error with leiden STREAM plot')
        except:
                print('Not enouch nearest')


        #--------------------------------------------------------------------------------------------------------------
        print('Progeny plotting')

        # Generate plots using progeny
        sample_adata.raw = sample_adata ## rename so it matches the pipeline

        model = progeny.load_model(
            organism='Human', # If working with mouse, set to Mouse
            top=1000          # For sc we recommend ~1k target genes since there are dropouts
        )

        progeny.run(sample_adata,        # Data to use
                    model,        # PROGENy network
                    center=True,  # Center gene expression by mean per cell
                    num_perm=100, # Simulate m random activities
                    norm=True,    # Normalize by number of edges to correct for large regulons
                    scale=True,   # Scale values per feature so that values can be compared across cells
                    use_raw=True, # Use raw sample_adata, where we have the lognorm gene expression
                    min_size=5    # Pathways with less than 5 targets will be ignored
                   )

        #Let's plot the obtained activities on the previously annotated clusters
        pw_sample_adata = progeny.extract(sample_adata)
        sc.set_figure_params(figsize=(4,4))
        sc.pl.umap(pw_sample_adata, color=pw_sample_adata.var.index, vmin=-4, vmax=5, cmap='coolwarm', ncols=3, save=i+'_activities.pdf')

        # Here is a summary per cell type:
        fig, axes = plt.subplots(5,3, figsize=(3*3,5*3), tight_layout=True)
        axes = axes.flatten()
        for pw,ax in zip(pw_sample_adata.var.index, axes):
            sc.pl.violin(pw_sample_adata, keys=pw, groupby='leiden', stripplot=False, ax=ax, show=False)
            ax.set_xlabel('')
            ax.tick_params(axis='x', rotation=90, labelsize=10)
        axes[-1].set_visible(False)
        fig.tight_layout()
        plt.savefig(figpath+i+'violin.pdf',format='pdf')
        plt.clf()

        # Comparison of pathway activities between groups
        # Once pathway activities are predicted, we can perform comparison tests to check if there are differences between groups.
        # To do, we have implemented a wrapper for the Wilcoxon rank-sum test. Let's find what pathways are activaten in cTEC cells:
        # df = progeny.rank_pws_groups(sample_adata, groupby='leiden', group='5')
        # df.to_csv(figpath+'wiloxon_ranksum.csv')

        # Let's plot the mean activity per cell type:
        pws = pw_sample_adata.var.index
        sc.pl.matrixplot(pw_sample_adata, pws, 'leiden', save=i+'_activity_per_cell.pdf', dendrogram=True, cmap='coolwarm', vmin=-1, vmax=1)
    except:
        print('Not completed all: '+i)
