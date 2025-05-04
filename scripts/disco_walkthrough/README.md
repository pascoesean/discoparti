## Walkthrough for running DisCo Analysis

Generate celltypes from an anndata object.

Run through scripts in this order, paying attention to how different preprocessing may 
be necessary for your data:

1. `preprocessing_and_metadata.py`: appends metadata, filters cells
2. `normalization_and_integration.py`: normalizes, (optionally) integrates, and calculates neighborhoods
3. `titrate_clusters_dotplot.py`: creates dotplots of different clustering resolutions to assess
which is best for discrete chunking
4. `investigate_celltypes.py`: Look at celltypes that have been assigned manually, compare to orginally published (if applicable)