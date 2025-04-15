# srun -N 1 --mem=25G --time=02:00:00 --exclude=c17 --pty bash
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/conda.sh" source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
# mamba activate pancre-yas-env

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation
from scipy.sparse import csr_matrix
import seaborn as sns
import os
import glob
import sys

from discoparti.scripts.functions.figure_functions import getSetup, genFigure


pbmc_data_df = sc.read_csv("discoparti/data/raw/counts_matrix.csv")

counts = csr_matrix(pbmc_data_df.X.T)

print("THE COUNTS", counts)


pbmc_data = anndata.AnnData(counts)
pbmc_data.var_names = pbmc_data_df.obs_names
pbmc_data.obs_names = pbmc_data_df.var_names


metadata = pd.read_csv("discoparti/data/raw/pbmc_metadata.csv")

pbmc_data.obs = metadata
#pbmc_data.obs = pbmc_data.obs.join(metadata, how='left', validate='one_to_one')


#### PREPROCESSING

# quality control -- remove mitochondiral + hemoglobin + ribosomes


# mitochondrial genes, "MT-" for human, "Mt-" for mouse
pbmc_data.var["mt"] = pbmc_data.var_names.str.startswith("MT-")
# ribosomal genes
pbmc_data.var["ribo"] = pbmc_data.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
pbmc_data.var["hb"] = pbmc_data.var_names.str.contains("^HB[^(P)]")

print("IM CALCULATING QC")

sc.pp.calculate_qc_metrics(
    pbmc_data, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True, percent_top=[20]
)


#sc.pp.filter_cells(pbmc_data, min_genes=100)



def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

pbmc_data.obs["outlier"] = (
    is_outlier(pbmc_data, "log1p_total_counts", 5)
    | is_outlier(pbmc_data, "log1p_n_genes_by_counts", 5)
    | is_outlier(pbmc_data, "pct_counts_in_top_20_genes", 5)
)

print("I HAVE THIS MANY OUTLIERS:", pbmc_data.obs.outlier.value_counts())

pbmc_data.obs["mt_outlier"] = is_outlier(pbmc_data, "pct_counts_mt", 3)

print("MY MITOCHONDRIAL UPPER THRESHOLD IS:", np.median(pbmc_data.obs['pct_counts_mt'] + 3*median_abs_deviation(pbmc_data.obs['pct_counts_mt'])))


print("I HAVE THIS MANY MITOCHONDRIAL OUTLIERS:", pbmc_data.obs.mt_outlier.value_counts())


print(f"Total number of cells: {pbmc_data.n_obs}")

pbmc_data_qc = pbmc_data[(~pbmc_data.obs.outlier) & (~pbmc_data.obs.mt_outlier)].copy()

print(f"Number of cells after filtering of low quality cells: {pbmc_data_qc.n_obs}")



# save out some figs
#sc.pl.highly_variable_genes(pbmc_data_nodoubs, save="_hvg.png")
#sc.pl.pca_variance_ratio(pbmc_data_nodoubs, n_pcs=50, log=False, save="pvr_.png")

print("RUNNING SCRUBLET")
sc.pp.scrublet(pbmc_data_qc)

print("\nNUMBER OF DOUBLETS:", pbmc_data_qc.obs['predicted_doublet'].sum())

pbmc_data_qc.write_h5ad("discoparti/data/processed/pbmc_qc.h5ad")

def makeFigure():
    ax, f = getSetup((15,20), (3,2))

    sns.histplot(pbmc_data.obs["total_counts"], bins=100, kde=False, ax=ax[0])
    sns.histplot(pbmc_data.obs["pct_counts_ribo"], bins=100, kde=False, ax=ax[1])
    sc.pl.violin(pbmc_data, "pct_counts_mt", ax=ax[2])
    sc.pl.scatter(pbmc_data, "total_counts", "n_genes_by_counts", color="pct_counts_mt", ax=ax[3])
    sc.pl.scatter(pbmc_data_qc, "total_counts", "n_genes_by_counts", color="pct_counts_mt", ax=ax[4])

    return f


ff = makeFigure()

genFigure(ff, "quality_control")
