# srun -N 1 --mem=40G --time=02:00:00 --exclude=c17 --pty bash
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/conda.sh" source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
# mamba activate pancre-yas-env

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import os
import glob
import sys
import scanpy.external as sce 
import harmonypy as hm 

from discoparti.scripts.functions.figure_functions import getSetup, genFigure


pbmc_data = sc.read_h5ad("discoparti/data/processed/pbmc_qc.h5ad")

print(pbmc_data)

print(pbmc_data.obs)

## NORMALIZATION + FEATURE SELECTION

# Saving count data
pbmc_data.layers["counts"] = pbmc_data.X.copy()

# Normalizing to 10,000 counts 
sc.pp.normalize_total(pbmc_data, target_sum=10000)
# Logarithmize the data
sc.pp.log1p(pbmc_data)

print("X MATRIX\n",pbmc_data.X)


## DO SOME MORE FILTERING. NOT SUPER ALLOWED; DONT TELL ANYONE

# remove ribosomal genes
pbmc_data = pbmc_data[:, (~pbmc_data.var.ribo)].copy()


sc.pp.highly_variable_genes(pbmc_data, n_top_genes=1500, layer = "counts", flavor='seurat_v3')

## DIMENSIONALITY REDUCTION + HARMONY
sc.pp.pca(pbmc_data, mask_var = 'highly_variable')

#print("RUNNING HARMONY\n\n")
#sce.pp.harmony_integrate(pbmc_data, 'donor', max_iter_harmony=30)


#sc.pp.neighbors(pbmc_data, use_rep="X_pca_harmony")
sc.pp.neighbors(pbmc_data)

sc.tl.umap(pbmc_data)

pbmc_data.write_h5ad("pbmc_t1d/data/processed/pbmc_notharmonized_noribo.h5ad")

# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets

print("RUNNING LEIDEN\n\n")

sc.tl.leiden(pbmc_data, flavor="igraph", resolution=0.25)


# save out some figs
sc.pl.highly_variable_genes(pbmc_data, save="_pbmc_hvg.png")
sc.pl.pca_variance_ratio(pbmc_data, n_pcs=50, log=False, save="_pbmc_pvr_.png")

def makeFigure():
    ax, f = getSetup((12,28), (5,2))

    sc.pl.pca(pbmc_data, color="Sample_ID", ax=ax[0])
    sc.pl.pca(pbmc_data, color='leiden', ax=ax[1])
    sc.pl.umap(pbmc_data, color="Sample_ID", size=2, ax=ax[2])
    sc.pl.umap(pbmc_data, color="leiden", ax=ax[3])
    sc.pl.umap(pbmc_data, color='COND', ax=ax[4])
    sc.pl.umap(pbmc_data, color='Cluster_Annotation_Merged', ax=ax[5])
    sc.pl.umap(pbmc_data, color='total_counts', ax=ax[6])
    sc.pl.umap(pbmc_data, color='pct_counts_ribo', ax=ax[7])
    sc.pl.umap(pbmc_data, color='pct_counts_mt', ax=ax[8])
    sc.pl.umap(pbmc_data, color='pct_counts_hb', ax=ax[9])

    return f


ff = makeFigure()

genFigure(ff, "inital_clustering_noribo")
