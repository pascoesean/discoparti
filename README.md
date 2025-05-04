# :mirror_ball: :tada: DisCo ParTI 

### Repository for reproducing Discrete and Continuous Analysis with Pareto Task Inference (DisCo ParTI) for single cell RNA sequencing data.

Targeted single cell analyses that combine methods for studying discrete and continuous variance have the power to reveal specific molecular changes in disease. Here, we propose ***Dis***crete and ***Co***ntinuous analysis through ***P***areto ***T***ask ***I***nference (***DisCo ParTI***), a computational analysis method that pairs cell type proportion exploration with Pareto optimization to reveal coordinated perturbations associated with different cell populations. This is based on the [initial ParTI method](https://www.weizmann.ac.il/mcb/alon/download/pareto-task-inference-parti-method), pioneered for biology by the Uri Alon lab ([Korem et al. 2015](http://www.weizmann.ac.il/mcb/UriAlon/sites/mcb.UriAlon/files/korem_2015_-_geometry_of_the_gene_expression_space.pdf)).

This repository contains the code itself, as well as an example runthrough from two datasets from Type I Diabetes patients. DisCo ParTI recapitulates the known ablation of β cells in T1D patients and diseased cells are shown to preferentially resolve near known T1D associated gene signatures, such as HLA-DQ1. Additionally, we find a unique subset of islet-proximal acinar cells to separate by disease phenotype and are thus able to resolve finetuned molecular shifts along with broad correlates of disease.

All custom functions are located in `scripts/functions/`, and can be used à la carte or in the example pipeline presented (`scripts/DisCo_ParTI_runthrough.R`).

Files in `\data\` are not uploaded to github because of size, but can be sent upon request. Two are included for ease of testing cell type proportion plotters. Scripts are able to completely reproduce pancreas figures by downloading fastq.gz files from the HPAP consortium. Seurat object for the PBMC data from [Honardoost et al., 2024](https://doi.org/10.1186/s13073-024-01300-z) is available under synapse accession code [syn53641849](https://www.synapse.org/#!Synapse:syn53641849).

Information about python packages + versions used in `scripts/disco_walkthrough` is included in `py_package_info.txt`. R packages required include `Seurat`, [`ParetoTI`](https://github.com/vitkl/ParetoTI), `patchwork`, and the [`tidyverse`](https://www.tidyverse.org/) suite. All R packages were up to date as of 14 April 2025.

