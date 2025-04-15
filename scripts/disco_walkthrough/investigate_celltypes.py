# srun -N 1 --mem=25G --time=02:00:00 --exclude=c17 --pty bash
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

from discoparti.scripts.functions.figure_functions import getSetup, genFigure


pbmc_data = sc.read_h5ad("discoparti/data/processed/pbmc_notharmonized_noribo.h5ad")


marker_genes = {
    "T cells": ["CD3E","CD3D", "IL7R", "CCR7", "TRAT1", "GZMB"], 
    "B Cells": ["CD19", "MS4A1", "CD79A", "IGJ"],
    "CD14+ Mono": ["CD14", "CD68", "FCGR1A", "S100A9"],
    "CD16+ Mono": ["FCGR3A","CX3CR1", "C1QA"],
    "DCs": ["CD74", "NPC2"],
    "NK Cells": ["NKG7", "GNLY", "FCER1G", "CCL5"],
    "Plasma Cells": ["CD38", "MZB1"],
    "Platelet contam": ["PPBP", "PF4", "ITGA2B"],
    "IFN Signature": ["IFIT2", "ISG15", "OASL"]
}




sc.tl.leiden(pbmc_data, flavor="igraph", resolution=0.60, key_added=f"leiden_res_0.60")

# assign cell types
pbmc_data.obs["cell_type_0.60"] = pbmc_data.obs["leiden_res_0.60"].map(
    {
        "0": "DCs",
        "1": "B Cells",
        "2": "T cells",
        "3": "ISGhi T Cells",
        "4": "NKT Cells",
        "5": "T cells",
        "6": "CD14+ monocytes",
        "7": "CD16+ monocytes",
        "8": "NK cells",
        "9": "Platelet Contamination"
    }
)


pbmc_data.obs.to_csv('discoparti/data/processed/pbmc_metadata_w_0.60cts_noribo.csv')

pbmc_data.write_h5ad("discoparti/data/processed/pbmc_ct0.60_noribo.h5ad")

# list partially from this paper (https://www.nature.com/articles/s41467-024-49724-w), partially from AI google
# partially from cellmarker 2.0. but some of those aren't helpful

def make_blank(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    return None

def makeFigure():
    ax, f = getSetup((8,15), (3,1))

    #make_blank(ax[1])
    #make_blank(ax[4])
    #make_blank(ax[7])

    #sc.pl.umap(pbmc_data, color=f"leiden_res_0.60", ax=ax[0])
    sc.pl.dotplot(pbmc_data, marker_genes, groupby=f"leiden_res_0.60", standard_scale='var', ax=ax[0])
    #sc.pl.umap(pbmc_data, color="cell_type_0.60", ax=ax[1])
    sc.pl.dotplot(pbmc_data, marker_genes, groupby=f"cell_type_0.60", standard_scale='var', ax=ax[1])
    #sc.pl.umap(pbmc_data, color="Cluster_Annotation_Merged", ax=ax[2])
    sc.pl.dotplot(pbmc_data, marker_genes, groupby=f"Cluster_Annotation_Merged", standard_scale='var', ax=ax[2])

    return f 

genFigure(makeFigure(), "investigate_celltypes_noribo")
