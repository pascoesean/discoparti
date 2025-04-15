# load packages
library(Seurat)
library(dplyr)
library(patchwork)

# Load in the data
pbmc <- readRDS("/Users/duyle/Library/CloudStorage/OneDrive-SharedLibraries-MassachusettsInstituteofTechnology/Sean Pascoe - 20.440_final/T1D_Seurat_Object_Final.rds")
pbmc

# Subset the data because too big LOL
cells_to_keep <- sample(colnames(pbmcs), 5000)
pbmc_subset <- subset(pbmcs, cells = cells_to_keep)
rm(pbmcs)

# Split data into healthy and t1d subsets
# split_by_cond takes a seurat object with COND in metadata labeling healthy
# controls and t1d patients. 
# split_by_cond returns a list object with (healthy, t1d) seurat objects
split_by_cond <- function(seuratObject) {
  healthy <- subset(seuratObject, subset = COND == "H")
  t1d <- subset(seuratObject, subset = COND == "T1D")
  return (list(healthy=healthy, t1d=t1d))
}

# Normalize and preprocess the data
preprocess <- function(seuratObject) {
  seuratObject <- NormalizeData(seuratObject)
  seuratObject <- FindVariableFeatures(seuratObject)
  seuratObject <- ScaleData(seuratObject)
  seuratObject <- RunPCA(seuratObject)
  seuratObject <- FindNeighbors(seuratObject, dims=1:10)
  
  return(seuratObject)
}

# Apply Leiden Cluster and view clusters
cluster <- function(seuratObject) {
  seuratObject <- FindClusters(seuratObject, resolution = 0.8, algorithm = 4)
  print(table(seuratObject$seurat_clusters))
  seuratObject <- RunUMAP(seuratObject, dims = 1:10)
  print(DimPlot(seuratObject, group.by = "Sample_ID", label = TRUE))
  print(DimPlot(seuratObject, group.by = "Cluster_Annotation_All", label = TRUE))
  print(DimPlot(seuratObject, group.by = "seurat_clusters", label = TRUE))
  print( print(DimPlot(seuratObject, group.by = "COND", label = TRUE)))
  
  return(seuratObject)
}

# splits <- split_by_cond(pbmc_ss)
# pbmc_healthy <- splits$healthy
# pbmc_t1d <- splits$t1d
# pbmc_healthy <- preprocess(pbmc_healthy)
# pbmc_healthy <- cluster(pbmc_healthy)
# pbmc_t1d <- preprocess(pbmc_t1d)
# pbmc_t1d <- cluster(pbmc_t1d)


# find markers for each cluster
find_markers <- function(seuratObject) { 
  markers <- FindAllMarkers(seuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # View top markers per cluster
  top_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  print(top_markers)
}

pbmc_subset <- preprocess(pbmc_subset)
pbmc_subset <- cluster(pbmc_subset)
