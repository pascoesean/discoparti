## Scripts relevant for reproducing pipeline

Directory of scripts needed to reproduce pipeline. 

Celltyping analysis is first performed in python, with scripts included to download pancreas data (`celltyping\`).

A full runthrough of cell type proportion analysis and ParTI is provided in `DisCo_ParTI_runthrough.R`, with some example data uploaded in `github_data/processed/`

Parti polytope fitting is done in MATLAB using the [ParTI package](https://www.weizmann.ac.il/mcb/alon/download/pareto-task-inference-parti-method) from the Uri Alon lab. Scripts for our fitting and subsequent enrichment analysis are available in (`parti/`)

Both parts reference common functions that can be found in the `functions` subdirectory.