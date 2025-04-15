# 

library(tidyverse)

cluster_genemarkers_pbmcs <- read_csv("data/cluster_titration_genesets_res_0.60_noribo.csv") |>
  filter(group == 4) |>
  head(n=100) |>
  write_csv("geneslol.csv")
  View()

pbmc_metadata <- read_csv("pbmc_metadata.csv")
filtering_criteria <- read_csv("data/pancreas_filtering_criteria.csv")

filtering_criteria |>
  pivot_longer(cols = pct_counts_mt:pct_counts_hb,names_to = "measure", values_to = "percent") |>
  ggplot(aes(x = percent)) +
  geom_histogram(bins=500) +
  facet_wrap(~measure, nrow=3, scales = "free_y") +
  theme_minimal()

### PANCREAS CELL TYPE COUNTS ----

metadata_wcts <- read_csv("data/metadata_w_0.10cts.csv") 
  

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




### PBMC CELL TYPE COUNTS ----

pbmcmetadata_wcts <- read_csv("data/pbmc_metadata_w_0.60cts_noribo.csv") 

cells_per_donor <- pbmcmetadata_wcts |>
  group_by(Sample_ID) |>
  summarize(total_cells_per_donor = n())

pbmcct_percentages <- pbmcmetadata_wcts |>
  group_by(Sample_ID, COND, cell_type_0.60) |>
  summarize(ct_counts = n()) |>
  left_join(cells_per_donor, by=join_by("Sample_ID")) |>
  mutate(ct_percent = ct_counts/total_cells_per_donor)

# ordered from top are diabetes bottom are healthy
pbmcct_percentages |>
  ggplot(aes(fill=cell_type_0.60, y=Sample_ID, x=ct_percent)) + 
  scale_fill_brewer(palette = 'Set1') +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Proportion of PBMCs", y = '', fill = "Cell Type") +
  ggpubr::theme_classic2()


# can also look at violins
pbmcct_percentages |>
  filter(cell_type_0.60 %in% c('NK cells', 'CD14+ monocytes', 'T cells', "B Cells")) |>
  ggplot(aes(x=cell_type_0.60, color=COND, y=ct_percent)) + 
  geom_boxplot(outliers = F) +
  scale_color_brewer(palette = 'Set1') +
  ggpubr::stat_compare_means(label = 'p.signif', method='wilcox.test', bracket.size = 2) +
  geom_point(position = position_jitterdodge()) +
  labs(x = 'Cell Type', y = "Proportion of PBMCs", color = "Disease Status") +
  ggpubr::theme_classic2()


rownames <- pbmcct_percentages |>
  select(COND) |>
  unique() |>
  ungroup() |>
  column_to_rownames('Sample_ID')

pbmcct_percentages |>
  ungroup() |>
  select(cell_type_0.60)

matrix <- pbmcct_percentages |>
  ungroup() |>
  select(Sample_ID, ct_percent, cell_type_0.60) |>
  pivot_wider(values_from = 'ct_percent', names_from = 'cell_type_0.60', values_fill = 0) |>
  column_to_rownames('Sample_ID') |>
  as.matrix()

pheatmap::pheatmap(matrix, annotation_row = rownames, cluster_rows = F, cluster_cols = F)


# retired -----

data <- read_csv("Desktop/metadata_wpreds.csv")
immunedata <- read_csv("Desktop/metadata_wpreds_IMMUNE.csv")

markers <- read_csv("cluster_titration_genesets_res_0.20.csv") |>
  filter(group == 6) |>
  select(names) |>
  head(n=50) |>
  write_csv("geneslol.csv")
  print(n=30) 
  View()

  
## RETIRED -----

data |>
  pivot_longer(cols=PP:delta, names_to='celltype', values_to = 'probability') |>
  ggplot(aes(x = probability)) +
  geom_histogram() +
  facet_wrap(~celltype)

data |>
  pivot_longer(cols=PP:delta, names_to='celltype', values_to = 'probability') |>
  filter(predicted_labels == celltype) |>
  ggplot(aes(x = probability)) +
  geom_histogram() +
  facet_wrap(~celltype)

immunedata |>
  pivot_longer(cols=`B cells`:`pDC precursor`, names_to='celltype', values_to = 'probability') |>
  filter(predicted_labels == celltype) |>
  ggplot(aes(x = probability)) +
  geom_histogram() +
  facet_wrap(~celltype)


both <- data |> left_join(immunedata, by = join_by(...1, donor, DiseaseStatus, SampleSex, SampleAge))


# PROBABLY MACROPHAGES?
both |>
  filter(predicted_labels.y == "Macrophages") |>
  filter(Macrophages > 0.99) |>
  select(!predicted_labels.y) |>
  pivot_longer(cols=PP:`pDC precursor`, names_to='celltype', values_to = 'probability') |>
  ggplot(aes(x = probability)) +
  geom_histogram() +
  facet_wrap(~celltype)
  