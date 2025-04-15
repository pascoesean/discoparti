# srun -N 1 --mem=25G --time=02:00:00 --exclude=c17 --pty bash
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/conda.sh" source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
# mamba activate pancre-yas-env

import scanpy as sc
import anndata
import pandas as pd
import os
import glob
import sys

from discoparti.scripts.functions.figure_functions import getSetup, genFigure

resolutions_to_test = [0.6, 0.75, 1]

pbmc_data = sc.read_h5ad("discoparti/data/processed/pbmc_notharmonized_noribo.h5ad")

marker_genes = {
    "T cells": ["CD3E","CD3D", "CD4", "CD40", "IL2RA", "GZMB"], 
    "B Cells": ["CD19", "MS4A1", "CD79A"],
    "CD14+ Mono": ["CD14", "CD68", "FCGR1A", "S100A9"],
    "CD16+ Mono": ["FCGR3A","CX3CR1", "C1QA"],
    "myeloid DCs": ["CD1C", "CD83"],
    "plasmacytoid DCs": ["CD48", "IL10"],
    "NK Cells": ["NKG7", "GNLY", "FCER1G", "CCL5"],
    "Plasma Cells": ["CD38", "MZB1"],
    "Platelet contam": ["PPBP", "PF4", "ITGA2B"],
    "IFN Signature": ["IFIT2", "ISG15", "OASL"]
}

# mostly from celltypist

def makeFigure():
    ax, f = getSetup((25,30), (5,2))

    for index, value in enumerate(resolutions_to_test):
        sc.tl.leiden(pbmc_data, flavor='igraph', resolution=value, key_added=f"leiden_res_{value:4.2f}")
        sc.pl.umap(pbmc_data, color=f"leiden_res_{value:4.2f}", ax=ax[index*2])
        # get dotplot for clusters
        #sc.pl.dotplot(pbmc_data, marker_genes, groupby=f"leiden_res_{value:4.2f}", standard_scale='var', ax=ax[index*2 + 1])
        sc.tl.rank_genes_groups(pbmc_data, groupby=f"leiden_res_{value:4.2f}", method="wilcoxon")
        sc.pl.rank_genes_groups_dotplot(
            pbmc_data, groupby=f"leiden_res_{value:4.2f}", standard_scale="var", n_genes=10, ax=ax[index*2 + 1]
        )
        sc.get.rank_genes_groups_df(pbmc_data,group=None).to_csv(f'discoparti/data/processed/cluster_titration_genesets_res_{value:4.2f}_noribo.csv')

    
    return f
 

genFigure(makeFigure(), "cluster_titration_dotplot_noribolist")
