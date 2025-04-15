# srun -N 1 --mem=25G --time=01:00:00 --exclude=c17 --pty bash
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/conda.sh" 
# source "/net/bmc-lab8/data/lab/dallgglab/users/spascoe/conda/etc/profile.d/mamba.sh"
# mamba activate pancre-yas-env

import scanpy as sc
import anndata
import os
import glob


# really just need the paired ones i want to use
donor_list = ["055", "064", "123", "130", "135", "137", "138", "141", "156", "167", "172",
              "131", "155", "082", "140", "024", "029", "037", "049", "146", "036", "080"]


'''
IF I WANT TO LOOK AT ALL DONORS:
donor_list = ["022", "024", "026", "029", "035", "036", "037", "040", "045", "049", "050", "052", "053",
 "054", "055", "056", "059", "063", "064", "072", "074", "074", "077", "080", "082", "092", "093", "095", 
 "097", "099", "101", "103", "105", "110", "114", "117", "118", "119", "122", "123", "128", "129", "130",
 "131", "135", "136", "137", "138", "139", "140", "141", "146", "155", "156", "157", "160", "167", "169",
 "171", "172", "173", "174", "175"]
 '''

adatas = {}

for donor_id in donor_list:
    string_to_search = "hpapdata/HPAP-"+ donor_id +"/STAR*Solo.out/Gene/filtered/"
    pathz = glob.glob(string_to_search)

    for index, string in enumerate(pathz):
        print("STRING TO SEARCH:", string)
        anndata_obj = sc.read_10x_mtx(string, var_names='gene_symbols', make_unique=True, cache=True)
        donor_run = "HPAP-" + donor_id + "-" + str(index)
        print("\n\nI JUST MADE AN ANNDATA OBJECT FOR :", donor_run)
        print(anndata_obj)
        anndata_obj.obs['donor'] = "HPAP-" + donor_id
        anndata_obj.obs['donor-run'] = donor_run
        sc.pp.scrublet(anndata_obj)
        adatas[donor_run] = anndata_obj



adata = anndata.concat(adatas, index_unique="_")
print(adata)
print(adata.obs)
print(adata.obs["donor-run"].value_counts())

adata.write_h5ad("discoparti/data/raw/pancdbcounts_matchedpairs.h5ad")

# next to do:
# 1. append metadata
# 2. preprocessing: remove dead (mito counts) + doublets (scrublet) https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
# 3. leiden cluster (big!)
# 4. plot umap and stuff
