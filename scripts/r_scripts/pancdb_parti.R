## PARTI -- PANCDB DATA

# To make sure R uses the correct conda enviroment you can run this when you start R:
reticulate::use_condaenv("ret-pcha-3.8", conda = "auto",
                         required = TRUE) # set TRUE to force R to use reticulate_PCHA

# load packages
library(tidyverse)
library(ParetoTI)
library(Seurat)

cellz = schard::h5ad2seurat('pairedpanc_ct0.10.h5ad')

delta_cells <- subset(cellz, subset = cell_type_0.10 == 'Delta Cells')

delta_cells <- FindVariableFeatures(delta_cells, nfeatures = 2000)
delta_cells <- ScaleData(delta_cells)
var_genes <- VariableFeatures(delta_cells)
seurat_df <- GetAssayData(delta_cells)[var_genes,]
hello <- as.data.frame(seurat_df)


has_diabetes <- delta_cells@meta.data$disease_pheno
names(has_diabetes) <- rownames(delta_cells@meta.data$cell_type_0.45)

                                                      
harmony_pca <- t(Embeddings(delta_cells, reduction = "Xpca_"))

find_arcs(harmony_pca)

# find archetypes
arc_4 <- fit_pch(harmony_pca)

arc_ks <- k_fit_pch(harmony_pca, ks = 2:6, check_installed = T,
                    verbose=T,
                   bootstrap = T, bootstrap_N = 200, maxiter = 1000,
                   bootstrap_type = "m", seed = 2543, 
                   volume_ratio = "t_ratio", # set to "none" if too slow
                   delta=0, conv_crit = 1e-04, order_type = "align",
                   method="louvain",
                   sample_prop = 0.75)

# Show variance explained by a polytope with each k (cumulative)
plot_arc_var(arc_ks, type = "varexpl", point_size = 2, line_size = 1.5) + theme_bw()
