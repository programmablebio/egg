library(SummarizedExperiment)
library('anndata')
library('spatstat')
library('SingleCellExperiment')
library(symphony)
library(tidyverse)
library(data.table)
library(Matrix)
library(plyr)
library(dplyr)
library(ggrastr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(DescTools)
library(harmony)

####### fetal and cell atlas integration ########################################################


ref = readRDS('PanglaoDBplus_noCancer_QCed_counts.rds')
ref_meta = readRDS('PanglaoDBplus_noCancer_metadata.rds')
fetal_atlas <- read.table(file="fetal_ovary_dataset_antilog.txt", na.strings = "-", header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
fetal_atlas = Matrix(as.matrix(fetal_atlas),sparse=TRUE)

tissue_type_SL <- rep('Fetal Ovary Cells', 1488)

cell_subset = grep(paste(paste('^',key_value,sep=''),'_',sep=''), query_exp$barcodes, value=TRUE)

source <- rep('fetal', 1488)
ref_meta$source = 'ref'

sample = colnames(fetal_atlas)

sample_list = as.character(df$sample)
sample_list[grepl("F.",df$sample)] = 'Fetal ovary cells'
sample_list[grep("M.",df$sample)] = 'Fetal spermatogonial cells'
tissue_type_SL = unlist(sample_list)

df <- data.frame(tissue_type_SL, sample, source)
str(df$tissue_type_SL)
common_cols <- intersect(colnames(ref_meta), colnames(df))
full_meta = rbind(
  subset(ref_meta, select = common_cols),
  subset(df, select = common_cols)
)

common_rows <- intersect(rownames(ref), rownames(fetal_atlas))

fetal_atlas_s = fetal_atlas[rownames(fetal_atlas) %in% common_rows]

ref_s = ref[common_rows ,]
fetal_atlas_s = fetal_atlas[common_rows ,]

library(Seurat)
###ref_s = NormalizeData(ref_s)
###fetal_atlas_s = NormalizeData(fetal_atlas_s)

full_ref = cbind(ref_s,fetal_atlas_s)

# subset to genes in our query data:
query_exp = readRDS('/newvolume/analysis/m1.rds')
genes =  read.csv('/newvolume/analysis/1milfinal/all_genes.csv')
cells =  read.csv('/newvolume/analysis/1milfinal/cell_metadata.csv')
query_metadata =  cells
mt.genes <- rownames(query_exp)[grep("^MT-",rownames(query_exp))]
C<-GetAssayData(object = query_exp, slot = "counts")
percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
query_exp <- AddMetaData(query_exp, percent.mito, col.name = "percent.mt")
query_exp <- subset(query_exp, subset = percent.mt < 15)
query_exp <- subset(query_exp, subset = nFeature_RNA < 8000 & nCount_RNA < 50000)
########################################################################################
common_rows <- intersect(rownames(full_ref), rownames(query_exp))
full_ref = NormalizeData(full_ref)
full_ref = full_ref[common_rows,]

common_rows <- intersect(rownames(full_ref), rownames(query_exp))

####var_genes = vargenes_vst(full_ref, topn = 2000)
var_genes = vargenes_vst(full_ref, groups = as.character(full_meta[['source']]), topn = 2000)
ref_exp = full_ref[var_genes, ]

vargenes_means_sds = tibble(symbol = var_genes, mean = Matrix::rowMeans(ref_exp))
vargenes_means_sds$stddev = singlecellmethods::rowSDs(ref_exp, vargenes_means_sds$mean)

ref_exp_scaled = singlecellmethods::scaleDataWithStats(ref_exp, vargenes_means_sds$mean, vargenes_means_sds$stddev, 1)

library("irlba")
set.seed(0)

s = irlba(ref_exp_scaled, nv = 20)
Z_pca_ref = diag(s$d) %*% t(s$v) # [pcs by cells]
loadings = s$u

ref_harmObj = HarmonyMatrix(
        data_mat=t(Z_pca_ref),  ## PCA embedding matrix of cells
        meta_data=full_meta, ## dataframe with cell labels
        vars_use=c('source'),
        theta = 0,
        nclust = 100,             ## number of clusters in Harmony model
        return_object = TRUE,     ## return the full Harmony model object
        do_pca = FALSE            ## don't recompute PCs
)

reference = symphony::buildReferenceFromHarmonyObj(
                           ref_harmObj,            # output object from HarmonyMatrix()
                           full_meta,           # reference cell metadata
                           vargenes_means_sds,     # gene names, means, and std devs for scaling
                           loadings,               # genes x PCs matrix
                           verbose = TRUE,         # verbose output
                           do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
                           save_uwot_path = './full_reference_uwot3')

saveRDS(reference, './full_reference_norm_fix.rds')


####### symphony map query ########################################################


plotBasic = function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                        title = 'Query',         # Plot title
                        color.by = 'cell_type',  # metadata column name for coloring
                        facet.by = NULL,         # (optional) metadata column name for faceting
                        color.mapping = NULL,    # custom color mapping
                        legend.position = 'bottom') {  # Show cell type legend

    p = umap_labels %>%
            dplyr::sample_frac(1L) %>% # permute rows randomly
            ggplot(aes(x = UMAP1, y = UMAP2)) +
            geom_point_rast(aes(col = get(color.by)), size = 0.3, stroke = 0.2, shape = 16)
        if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }

    # Default formatting
    p = p + theme_bw() +
            labs(title = title, color = color.by) +
            theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position=legend.position) +
            theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) +
            guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')

    if(!is.null(facet.by)) {
        p = p + facet_wrap(~get(facet.by)) +
                theme(strip.text.x = element_text(size = 12)) }
    return(p)
}

adata = read_h5ad('/newvolume/analysis/scvelo-objs/mega_scvelo_adata.h5ad')
adata = read_h5ad('/newvolume/analysis/scvelo-objs/M1_adata.h5ad')
as.data.frame.matrix(adata)
adata = readH5AD('/newvolume/analysis/scvelo-objs/mega_scvelo_adata.h5ad')

adata = LoadH5Seurat('/newvolume/analysis/scvelo-objs/mega_scvelo_adata.h5ad')

library('Seurat')
library(scater)
library(sceasy)
seurat <- Convert(adata, to = "seurat")

sceasy::convertFormat('/newvolume/analysis/scvelo-objs/mega_scvelo_adata.h5ad', from="anndata", to="seurat",
                       outFile='/newvolume/analysis/m1.rds')
##############

reference = readRDS('./full_reference_norm_fix.rds')
#######################

query_exp = readRDS('/newvolume/analysis/m1.rds')
genes =  read.csv('/newvolume/analysis/1milfinal/all_genes.csv')
cells =  read.csv('/newvolume/analysis/1milfinal/cell_metadata.csv')
query_metadata =  cells

mt.genes <- rownames(query_exp)[grep("^MT-",rownames(query_exp))]
C<-GetAssayData(object = query_exp, slot = "counts")

###  filter

percent.mito <- colSums(C[mt.genes,])/Matrix::colSums(C)*100
query_exp <- AddMetaData(query_exp, percent.mito, col.name = "percent.mt")
query_exp <- subset(query_exp, subset = percent.mt < 15)
query_exp <- subset(query_exp, subset = nFeature_RNA < 8000 & nCount_RNA < 50000)

x = list(query_exp$barcodes)

pdf(file="/newvolume/analysis/atlas_integration/reference_plot.pdf")
plotReference(reference,
                   as.density = TRUE,      # plot density or individual cells
                   bins = 10,              # if density, nbins parameter for stat_density_2d
                   bandwidth = 1.5,        # if density, bandwidth parameter for stat_density_2d
                   title = "Reference Atlas",    # Plot title
                   color.by = 'tissue_type_SL', # metadata column name for cell type labels
                   show.legend = TRUE,     # Show cell type legend
                   show.labels = TRUE,     # Show cell type labels
                   show.centroids = FALSE) # Plot soft cluster centroid locations)
dev.off()
key_value = 96
my_range = 1:96

### run all samples

percent_list_v2 = rep(-2,96)
percent_list = list()
for (indexing in 1:96) {
    key_value = indexing
    if (key_value < 10){key_value = paste0('0',key_value)}
    cell_subset = grep(paste(paste('^',key_value,sep=''),'_',sep=''), query_exp$barcodes, value=TRUE)
    print(str(cell_subset))
    sample_query  <- subset(query_exp, cells = cell_subset)

    query = mapQuery(sample_query[['RNA']]@counts,             # query gene expression (genes x cells)
                     sample_query@meta.data,        # query metadata (cells x attributes)
                     reference,             # Symphony reference object
                     do_normalize = TRUE,
                     do_umap = TRUE)

    query = knnPredict(query, reference, reference$meta_data$tissue_type_SL, k = 10, confidence = FALSE)
    tabulated = table(query$meta_data$cell_type_pred_knn)
    percent = 100.0*tabulated["Fetal ovary cells"]/length(query$meta_data$cell_type_pred_knn)
    percent_list <- append(percent_list, percent)
    print(percent)
    percent_list_v2[indexing] = percent

    # Sync the column names for both data frames
    reference$meta_data$cell_type_pred_knn = NA
    reference$meta_data$cell_type_pred_knn_prob = NA

    query$meta_data$ref_query = key_value
    query$meta_data$ref_query = key_value

    reference$meta_data$ref_query = 'reference'
    query$meta_data$ref_query = 'query'
    query$meta_data$sample = query$meta_data$barcodes
    query$meta_data$tissue_type_SL = query$meta_data$cell_type_pred_knn

    common_cols <- intersect(colnames(query$meta_data), colnames(reference$meta_data))
    meta_data_combined = rbind(
      subset(query$meta_data, select = common_cols),
      subset(reference$meta_data, select = common_cols))

    # Add the UMAP coordinates to the metadata

    umap_combined = rbind(query$umap, reference$umap$embedding)
    umap_combined_labels = cbind(meta_data_combined, umap_combined)


    }

######### plots


pdf(file=paste(paste("/newvolume/analysis/atlas_integration/",key_value,sep=''),'_test.pdf',sep=''))
plotBasic(umap_combined_labels, title = paste(round(percent, digits = 2),'% Fetal Ovary Cluster',sep=''),
        color.by = 'cell_type_pred_knn', legend.position='bottom')
dev.off()


pdf(file=paste(paste("/newvolume/analysis/atlas_integration/",key_value,sep=''),'_tissuetype.pdf',sep=''))
plotBasic(umap_combined_labels, title = paste(round(percent, digits = 2),'% Fetal Ovary Cluster',sep=''),
        color.by = 'tissue_type_SL', facet.by = 'ref_query', legend.position='bottom')
dev.off()

pdf(file=paste(paste("/newvolume/analysis/atlas_integration/",key_value,sep=''),'_knn.pdf',sep=''))
plotBasic(umap_combined_labels, title = paste(round(percent, digits = 2),'% Fetal Ovary Cluster',sep=''), color.by = 'cell_type_pred_knn', facet.by = 'ref_query', legend.position='bottom')
dev.off()
