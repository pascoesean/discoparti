library(tidyverse)
library(Seurat)


#### 
# CELIAC -----
#####

celiac_metadata <- readxl::read_excel("ImmuneCellsSmallIntestineCeliac_metadata_16-05-2023.xlsx",
                                      sheet="Donor organism",range = "A4:R25") |>
  filter(!is.na(donor_organism.biomaterial_core.biomaterial_id))

celiac_count_matrix <- readxl::read_excel("Source_data_1_scRNAseq_expression_count_matrix.xlsx")  |>
  rename("gene" = "...1") |> 
  column_to_rownames("gene") |>
  as.matrix()

# to generate metadata; need to get the list of rownames; then append schtuff <3

celiac_cell_metadata <- colnames(celiac_count_matrix) |>
  as_tibble() |>
  rename("cellID" = "value") |>
  tidyr::separate(cellID, c("donor", "tissue", "cell", "number"), sep= "\\.", remove=FALSE) |>
  left_join(y = celiac_metadata, by= join_by("donor" == "donor_organism.biomaterial_core.biomaterial_id")) |>
  select(cellID, donor, donor_organism.biomaterial_core.biomaterial_name, 
         donor_organism.biomaterial_core.biomaterial_description, donor_organism.sex,
         donor_organism.development_stage.text, donor_organism.diseases.text,
         donor_organism.medical_history.test_results, donor_organism.medical_history.treatment) |>
  mutate(donor_organism.medical_history.treatment = case_when(
    is.na(donor_organism.medical_history.treatment) ~ "Not gluten free",
    TRUE ~ donor_organism.medical_history.treatment
  )) |>
  column_to_rownames("cellID") |>
  as.data.frame()


celiac_seurat <- CreateSeuratObject(counts = celiac_count_matrix, meta.data = celiac_cell_metadata,
                                    project = "celiac")

SaveSeuratRds(celiac_seurat, file = "celiac_seuratobj.rds")

