## Parti Polytope Fitting and Analysis

Polytope fitting can be done by feeding input data (normalized counts for each discrete unit of interest) into `polytope_fitting.m` and running in MATLAB. 

Phenotype enrichment analysis is performed as follows using `run_pea` (which references enrichment functions in the `functions` subdirectory).

1. Rank all cells from shortest to longest euclidean distance to that archetype

2. For each cell, add $s_{hit}$ to the running score if the cell is in the phenotype of interest. Otherwise, subtract $s_{miss}$, where $s_{hit}$ is the absolute value of the mean-centered distance to that archetype scaled by the total score of observations in that phenotype, and $s_{miss}$ is one over the number of observations not in that phenotype, i.e.

$$ s_{hit} = \frac{|d_{i}|}{D_p} \text{, where } D_p = \sum_{i=1}^{p}{|d_i|}$$
$$ s_{miss} = \frac{1}{N_{\neg p}}$$
Here, $D_p$ is the summed absolute value of the mean centered distances for all $p$ values in the phenotype of interest, and $N_\neg p$ is the number of observations that are not in the phenotype of interest.


3. Report maximum deviation from zero of the running score as the normalized enrichment score.

Polytope plotting can be done using the `plot_parti` R function located in `functions/parti_output_wrangling.R`. This requires inputs from the two scripts in this folder.