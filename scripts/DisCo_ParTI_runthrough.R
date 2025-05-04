# Pareto Task Inference on T1D Data ----

## load packages, data, functions ----

library(Seurat)
library(patchwork)
library(ParetoTI)
library(tidyverse)
library(SingleCellExperiment)

source("scripts/functions/ct_props.R")
source("scripts/functions/Preprocessing.R")
source("scripts/functions/ParTI_funcs.R")



pbmc <- readRDS("/Users/duyle/Library/CloudStorage/OneDrive-SharedLibraries-MassachusettsInstituteofTechnology/Sean Pascoe - 20.440_final/discoparti/data/raw/T1D_Seurat_Object_Final.rds")


## initial processing -----


# Subset the data because too big LOL
cells_to_keep <- sample(colnames(pbmcs), 5000)
pbmc_subset <- subset(pbmcs, cells = cells_to_keep)
rm(pbmcs)


pbmc_subset <- preprocess(pbmc_subset)
pbmc_subset <- cluster(pbmc_subset)


## Cell Type proportions: these can be run with data on github ----

### PANCREAS CELL TYPE COUNTS ----

metadata_wcts <- read_csv("github_data/processed/metadata_w_0.10cts.csv") 
plot_pancreas_celltypes(metadata_wcts)

### PBMC CELL TYPE COUNTS ----

pbmcmetadata_wcts <- read_csv("github_data/processed/pbmc_metadata_w_0.60cts_noribo.csv") 
plot_pbmc_celltypes(pbmcmetadata_wcts)


## Run ParTI + generate figs ----
# these need larger files, request from spascoe@mit.edu or dale2@mit.edu

# need filepath for celltype labeled pancreas cells, and the cell type you want to look at as strings
run_everything('pairedpanc_ct0.10.h5ad', "Acinar Cells")


## Get + plot PCA loadings ----

pca_loadings <- get_pca_loadings("/data/pbmc_ct0.60_noribo.h5ad", "CD16+ monocytes")

plot_pca_loadings(pca_loadings)

