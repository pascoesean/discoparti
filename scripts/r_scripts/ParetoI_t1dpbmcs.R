# load packages
library(Seurat)
library(patchwork)
library(ParetoTI)
library(tidyverse)
library(SingleCellExperiment)

# Standard normalizing and preprocessing of the data
preprocess <- function(seuratObject) {
  seuratObject <- NormalizeData(seuratObject)
  seuratObject <- FindVariableFeatures(seuratObject)
  seuratObject <- ScaleData(seuratObject)
  seuratObject <- RunPCA(seuratObject)
  print("Done preprocessing.")
  return(seuratObject)
}

# Fit multiple polytopes to pcas and plot variance as a function of increasing dimensions
find_arcs <- function(pca) {
  
  # fit pcas to different k-dimension polytopes using PCHA algorithm
  arc_ks <- k_fit_pch(pca, ks = 2:8, check_installed = T,
                     bootstrap = T, bootstrap_N = 100, maxiter = 100,
                     bootstrap_type = "m", seed = 2543, 
                     volume_ratio = "t_ratio", # set to "none" if too slow
                     delta=0, conv_crit = 1e-04, order_type = "align",
                     sample_prop = 0.75)

  # Show variance explained by a polytope with each k (cumulative)
  var1 <- plot_arc_var(arc_ks, type = "varexpl", point_size = 2, 
                       line_size = 1.5) + theme_bw()
  
  # Show variance explained by k-vertex model on top of k-1 model (each k separately)
  var2 <- plot_arc_var(arc_ks, type = "res_varexpl", point_size = 2, 
                       line_size = 1.5) + theme_bw()

  # Show variance in position of vertices obtained using bootstraping 
  # - use this to find largest k that has low variance
  var3 <- plot_arc_var(arc_ks, type = "total_var", point_size = 2, 
                       line_size = 1.5) + theme_bw() + 
                       ylab("Mean variance in position of vertices")

  plots <- var1 / var2 / var3 + ggtitle("Variance as a function of k archetypes")
  
  print(plots)
  
  print("Done calculating vertices")
}

# Find the smallest polytope that fits most of the data
find_best_arc <- function(pca, vertices, labs, tit) {
  # Fit a polytope with bootstraping of cells to see stability of positions
  arc <- fit_pch_bootstrap(pca, n = 200, sample_prop = 0.75, seed = 235,
                          noc = vertices, delta = 0, conv_crit = 1e-04, type = "m")

  # plot the data by disease
  p_pca <- plot_arc(arc_data = arc, data = pca, 
                   which_dimensions = 1:3, line_size = 1.5,
                   data_lab = labs, data_alpha = 0.5, 
                   text_size = 60, data_size = 6) 
  arc_plot <- plotly::layout(p_pca, title = tit)
  print(arc_plot)
  print("Done finding best arc.")
  return(arc)
}

# Determine the fit of the polytope to data significance by bootstrapping 
# data points and comparing t-ratio (volume ratios) finding the
# likelihood of random polytopes having smaller t-ratios than observed polytope
fit_p <- function(pca, arc) {
  print("In p function")
  pch_rand = randomise_fit_pch(pca, arc_data = arc,
                               n_rand = 500,
                               replace = FALSE, bootstrap_N = NA,
                               volume_ratio = "t_ratio",
                               maxiter = 100, delta = 0, conv_crit = 1e-4,
                               type = "m", clust_options = list(cores = 3))
  print(pch_rand)
  print("Done finding p.")
}

# Map genes to archetypes, map genes/archetypes to GO terms
go_analyze <- function(gene_exp, vertices, pca) {
  arc1 <- fit_pch(pca, volume_ratio = "t_ratio", maxiter = 500,
                         noc = vertices, delta = 0,
                         conv_crit = 1e-04)
  
  activ = measure_activity(gene_exp, # row names are assumed to be gene identifiers
                           which = "BP", return_as_matrix = F,
                           taxonomy_id = 10090, keytype = "ALIAS",
                           lower = 20, upper = 1000,
                           aucell_options = list(aucMaxRank = nrow(gene_exp) * 0.1,
                                                 binary = F, nCores = 3,
                                                 plotStats = FALSE))
  
  # Merge distances, gene expression and gene set activity into one matrix
  data_attr = merge_arch_dist(arc_data = arc1, data = pca, 
                              feature_data = as.matrix(gene_exp),
                              colData = activ,
                              dist_metric = c("euclidean", "arch_weights")[1],
                              colData_id = "cells", rank = F) 
  
  # Use Wilcox test to find genes maximally expressed in 10% closest to each vertex
  enriched_genes = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                          features = data_attr$features_col,
                                          bin_prop = 0.1, method = "BioQC")
  
  enriched_sets = find_decreasing_wilcox(data_attr$data, data_attr$arc_col,
                                         features = data_attr$colData_col,
                                         bin_prop = 0.1, method = "BioQC")
  
  # Take a look at top genes and functions for each archetype
  labs = get_top_decreasing(summary_genes = enriched_genes, summary_sets = enriched_sets,
                            cutoff_genes = 0.01, cutoff_sets = 0.05, 
                            cutoff_metric = "wilcoxon_p_val", 
                            p.adjust.method = "fdr",
                            order_by = "mean_diff", order_decreasing = T,
                            min_max_diff_cutoff_g = 0.4, min_max_diff_cutoff_f = 0.03)
  print("Done finding gene sets.")
}

# # Load in the data
# pbmc <- readRDS("T1D_Seurat_Object_Final.rds")
# 
# # Check unique cell type labels
# unique(pbmc@meta.data$Cluster_Annotation_All)
# 
# # Extract cells of interest
# pbmc_ss<- subset(pbmc, subset = Cluster_Annotation_All %in% c("CD8_CM"))
# rm(pbmc)
# 
# # preprocess data
# pbmc_ss <- preprocess(pbmc_ss)
# 
# # Extract metadata labels by condition for plotting
# labs <- pbmc_ss@meta.data$COND
# names(labs) <- rownames(pbmc_ss@meta.data)
# 
# # Get expression data for go analysis
# gene_mat <- pbmc_ss@assays$RNA@data
# 
# # Extract pcas for ParTI and transpose so pcs are rows, samples are columns
# pcs <- t(Embeddings(pbmc_ss, reduction="pca"))
# 
# # find arcs and plot polytope
# find_arcs(pcs)
# best_arc <- find_best_arc(pcs, 3, labs)
# 
# # analyze go
# go_analyze(gene_mat, 3, pcs)
# 
# # p test
# fit_p(pcs, best_arc)

run_everything <- function(fileName, cell_type) {
  seuratObject = schard::h5ad2seurat(fileName)
  ss <- subset(seuratObject, subset = cell_type_0.10 %in% c(cell_type))
  rm(seuratObject)
  #ss <- preprocess(ss)
  labs <- ss@meta.data$disease_pheno
  names(labs) <- rownames(ss@meta.data)
  gene_mat <- ss@assays$RNA@data
  
  pcs <- t(Embeddings(ss, reduction="Xpca_"))
  
  find_arcs(pcs)
  best_arc <- find_best_arc(pcs, 4, labs, cell_type)
  
  go_analyze(gene_mat, 4, pcs)
  
  fit_p(pcs, best_arc)
}

run_everything('pairedpanc_ct0.10.h5ad', "Acinar Cells")
