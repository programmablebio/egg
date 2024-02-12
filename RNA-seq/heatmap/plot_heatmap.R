library(tidyverse)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(stringr)
library(devtools)
library(ggfortify)
library(factoextra)
library(cluster)
library(RColorBrewer)
library(biomaRt)
library(gridtext)
library(tidyHeatmap)
library(viridis)

# read counts
logCounts <- read.table('log2(norm_counts+1)_2024-02-10.csv', header=TRUE, sep=',', check.names=FALSE)#TRUE)
logCounts <- logCounts[, !grepl("^X$", colnames(logCounts))]
logCounts$Geneid <- sub("\\..*", "", logCounts$Geneid)

# read genes
new_genes <- read.table('heatmap_genes_list_2024-02-10.csv', header=TRUE, sep=",", check.names =TRUE)

# convert Ensembl gene ID to gene symbol
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- logCounts$Geneid
gene_names <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = ensembl)

logCounts <- merge(logCounts, gene_names, by.x = "Geneid", by.y = "ensembl_gene_id", all.x = TRUE)
logCounts$Geneid <- logCounts$external_gene_name
logCounts <- logCounts[, !(names(logCounts) == "external_gene_name")]
colnames(logCounts)[colnames(logCounts) == "Geneid"] <- "Genes"

# remove other samples
logCounts <- na.omit(logCounts)
logCounts <- logCounts[!is.na(logCounts$Genes) & logCounts$Genes != "", ]
new_genes <- na.omit(new_genes)

# subset based on gene symbols
new_genes_logCounts <- subset(logCounts, Genes %in% new_genes$Genes)
new_genes <- subset(new_genes, Genes %in% unique(new_genes_logCounts$Genes))
new_genes_logCounts <- new_genes_logCounts[match(new_genes$Genes, new_genes_logCounts$Genes), ]
new_genes_logCounts[, "Genes"] <- new_genes_logCounts$Genes
new_genes_logCounts <- distinct(new_genes_logCounts)
rownames(new_genes_logCounts) <- new_genes_logCounts$Genes
new_genes_logCounts <- subset(new_genes_logCounts,select=-c(Genes))
new_genes_logCounts_mat <- data.matrix(new_genes_logCounts)

#save the filtered counts
write.csv(new_genes_logCounts_mat,"2024-02-10_new_genes_logCounts_mat.csv")

# NEW GENES HEATMAP
mat <- new_genes_logCounts_mat
new_genes_column <- mat[, 1]
new_genes <- new_genes[new_genes$Genes %in% rownames(mat), ]

#Set the color scale using Magma from the viridis package
breaks <- c(0,1,2,3,4,5,6,7,8,9,12,13.5,15)
colors_for_map <- c(magma(256)[1],magma(256)[26],magma(256)[51],magma(256)[76],magma(256)[101],magma(256)[126],magma(256)[151],magma(256)[176],magma(256)[201],magma(256)[226],magma(256)[241],magma(256)[256],"#FDFEEF")
colmap <- colorRamp2(breaks, colors_for_map)

colors_stage <- c("hiPSC" = "#000000",
                  "hPGCLC" = "grey18",
                  "RA-FGC" = "grey30",
                  "Oogonia" = "snow4",
                  "Maturing_Oocyte" = "snow3",
                  "Meiotic_Oocyte" = "snow2"
)

font_colors <- c("hiPSC" = "white",
                 "hPGCLC" = 'white',
                 "RA-FGC" = "white",
                 "Oogonia" = "white",
                 "Maturing_Oocyte" = "black",
                 "Meiotic_Oocyte" = "black"
)


new_genes_heat <- ComplexHeatmap::Heatmap(mat,
                           name = 'Log2FC',
                           cluster_rows = F, 
                           cluster_columns = F,
                           col = colmap,
                           row_names_gp = gpar(fontsize = 7.8),
                           show_row_names=TRUE,
                           column_names_rot = 90,
                           #row_names_gp = NULL,
                           column_names_gp = gpar(fontsize = 9),
                           #row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 10)),
                           column_names_max_height = max_text_height(colnames(mat), gp = gpar(fontsize = 4)),
                           heatmap_width = unit(9, "in"),
                           heatmap_height = unit(7, "in"),
                           heatmap_legend_param = list(title = gt_render("<span style='color:black'>**log2(Normalized Counts + 1)**</span>"),
                                                       title_gp = gpar(fontsize=16),
                                                       labels_gp = gpar(fontsize=16)),
                           row_split=factor(c(new_genes$Partition),levels=c("TF_Control","Pluripotency","PGCLC","Oogonia","Oocyte","Meiosis")),
                           row_title_rot = 0,
                           column_split= NULL,
                           row_gap = unit(3, "mm"),
                           column_gap = unit(15, "mm"),
                           row_title_gp = gpar(
                             fill = colors_stage,
                             col = font_colors,
                             border = c(rep("white",6))),
                           #row_gap_color = "black"
                          
)
new_genes_heat
png(file="log2(norm_counts+1)_heat_2024-02-10_magma.png", width=13, height=12, units='in', res=300)
save_pdf(new_genes_heat,"Figure6_log2(norm_counts+1)_heat_2024-02-10.pdf")
dev.off()

