# DisCo ParTI :mirror_ball::tada: 

Repository for reproducing Discrete and Continuous Analysis with Pareto Task Inference (DisCo ParTI) for single cell RNA sequencing data.

Targeted single cell analyses that combine methods for studying discrete and continuous variance have the power to reveal specific molecular changes in disease. Here, we propose *Dis*crete and *Co*ntinuous analysis through *P*areto *T*ask *I*nference (***DisCo ParTI***), a computational analysis method that pairs cell type proportion exploration with Pareto optimization to reveal coordinated perturbations associated with different cell populations. 


This repository contains the code itself, as well as an example runthrough from two datasets from Type I Diabetes patients. DisCo ParTI recapitulates the known ablation of β cells in T1D patients and diseased cells are shown to preferentially resolve near known T1D associated gene signatures, such as HLA-DQ1. Additionally, we find a unique subset of islet-proximal acinar cells to separate by disease phenotype and are thus able to resolve finetuned molecular shifts along with broad correlates of disease.

All custom functions are located in `scripts/functions/`, and can be used à la carte or in the example pipeline presented.

Files in `\data\` are not uploaded to github because of size, but can be sent upon request. Scripts are able to completely reproduce pancreas figures by downloading fastq.gz files from the HPAP consortium.