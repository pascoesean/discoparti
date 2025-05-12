# srun -N 1 --mem=25G --time=02:00:00 --exclude=c17 --pty bash
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/conda.sh" source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
# mamba activate pancre-yas-env


# script for exporting things to be run with parti
import scanpy as sc
import pandas as pd
import numpy as np
import os

pancreas_data = sc.read_h5ad("pancre-yas/data/processed/pairedpanc_ct0.10.h5ad")


pd.DataFrame(pancreas_data.X.toarray()).to_csv("pancre-yas/data/processed/whole_array.csv")
pd.Series(pancreas_data.var_names).to_csv("pancre-yas/data/processed/whole_genes.csv")
pancreas_data.obs['disease_pheno'].to_csv("pancre-yas/data/processed/whole_cell_labels.csv")

assert False
def save_subset(anndata, subset_name: str, column_name: str, output_dir: str):
    print("\nSAVING SUBSET:", subset_name, "\n")
    subset = anndata[anndata.obs[column_name] == subset_name]
    try:
        sc.pp.highly_variable_genes(subset, n_top_genes=1500, layer = "counts", flavor='seurat_v3', subset = True)
    except ValueError:
        sc.pp.highly_variable_genes(subset, n_top_genes=1500, subset = True)
    #subset_highlyvariablegenes = subset[:, subset.var["highly_variable"]]
    scaled_counts_subset = pd.DataFrame(subset.X.toarray())
    gene_names = pd.Series(subset.var_names)

    # save out
    try:
        os.mkdir(output_dir + subset_name)
    except FileExistsError:
        print("Directory '" + output_dir + subset_name + "' already exists. Overwriting <3")
    
    scaled_counts_subset.to_csv(output_dir + subset_name + "/matrix.csv", header=False, index=False)
    subset.obs['disease_pheno'].to_csv(output_dir + subset_name + "/cell_labels.csv")
    gene_names.to_csv(output_dir + subset_name + "/gene_names.csv")


celltype_list = pancreas_data.obs['cell_type_0.10'].unique().tolist()

for celltype in celltype_list:
    save_subset(pancreas_data, celltype, 'cell_type_0.10', "pancre-yas/data/processed/parti_inputs/")


