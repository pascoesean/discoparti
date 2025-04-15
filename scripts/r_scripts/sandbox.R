library(tidyverse)
library(Seurat)


data <- read_rds("multiomics_t1d_60k.rds")


# 
t1d_metadata <- data@meta.data |>
  select(disease_state, donor_id, sex, development_stage, cell_label) |>
  group_by(disease_state, cell_label) |>
  summarize(count = n())


data@meta.data |>
  select(disease_state, donor_id, sex) |>
  unique()
