# cell type proportion analysis functions

library(tidyverse)

plot_pancreas_celltypes <- function(metadata_wcts){
	cells_per_donor <- metadata_wcts |>
  		group_by(donor) |>
  		summarize(total_cells_per_donor = n())

	ct_percentages <- metadata_wcts |>
  		group_by(donor, disease_pheno, SampleSex, cell_type_0.10) |>
  		summarize(ct_counts = n()) |>
  		left_join(cells_per_donor, by=join_by("donor")) |>
  		mutate(ct_percent = ct_counts/total_cells_per_donor)


# ordered from top are diabetes bottom are healthy
	ct_percentages |>
  		mutate(donor = factor(donor, 
                        levels = c("HPAP-131", "HPAP-155", "HPAP-082", "HPAP-140", "HPAP-024", "HPAP-029", 
                                   "HPAP-037", "HPAP-049", "HPAP-146", "HPAP-036", "HPAP-080",
                                   "HPAP-055", "HPAP-064", "HPAP-123", "HPAP-130", "HPAP-135", "HPAP-137",
                                   "HPAP-138", "HPAP-141", "HPAP-156", "HPAP-167", "HPAP-172"
                                   ))) |>
  	ggplot(aes(fill=cell_type_0.10, y=donor, x=ct_percent)) + 
  	scale_fill_brewer(palette = 'Set1') +
  	geom_bar(position="fill", stat="identity") +
  	labs(x = "Proportion of all cells", y = "Donor", fill="Cell Type") +
  	ggpubr::theme_classic2()
}


plot_pbmc_celltypes <- function(pbmcmetadata_wcts){
	cells_per_donor <- pbmcmetadata_wcts |>
  		group_by(Sample_ID) |>
  		summarize(total_cells_per_donor = n())

	pbmcct_percentages <- pbmcmetadata_wcts |>
  		group_by(Sample_ID, COND, cell_type_0.60) |>
  		summarize(ct_counts = n()) |>
  		left_join(cells_per_donor, by=join_by("Sample_ID")) |>
 		 mutate(ct_percent = ct_counts/total_cells_per_donor)

	pbmcct_percentages |>
  		ggplot(aes(fill=cell_type_0.60, y=Sample_ID, x=ct_percent)) + 
 		 scale_fill_brewer(palette = 'Set1') +
  		geom_bar(position="fill", stat="identity") +
  		labs(x = "Proportion of PBMCs", y = '', fill = "Cell Type") +
  		ggpubr::theme_classic2()
}