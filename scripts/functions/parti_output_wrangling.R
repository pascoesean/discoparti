## plotting/wrangling parti outputs

# load packages ----
library(tidyverse)
library(plotly)

# load data -----

# FUNCTION DEFINITION ---- 
plot_parti <- function(pc_filepath, labels_filepath, enrichment_filepath,
                       gene_names_filepath, arch_orig_filepath, arch_filepath, pea_filepath,
                       title) 
  {

beta_cells_scores <- read_csv(pc_filepath, col_names = str_c("PC_", 1:1500))
cell_labels <- read_csv(labels_filepath) |>
  rownames_to_column("index")

enrichment <- read_csv(enrichment_filepath) |>
  mutate(archetype = str_c("archetype_", `archetype #`)) |>
  group_by(archetype) |>
  slice_min(`P value (Mann-Whitney)`, n = 5) |>
  select(archetype, `Feature Name`) |>
  group_by(archetype) |>
  mutate(top_terms = str_c(`Feature Name`, collapse = ', <br>')) |>
  select(archetype, top_terms) |>
  unique()


gene_names <- read_csv(gene_names_filepath) |>
  select(`0`) |>
  as.vector()

archOrig <- read_csv(arch_orig_filepath, col_names = gene_names$`0`) |>
  rownames_to_column('index') |>
  mutate(archetype = paste0('archetype_', index)) |>
  select(!index) |>
  pivot_longer(cols = 1:1500, names_to = 'gene', values_to = 'expression_level')

mean_expression_bygene <- archOrig |>
  group_by(gene) |>
  summarize(mean_expression = mean(expression_level))

top_genes <- archOrig |>
  left_join(y = mean_expression_bygene) |>
  mutate(expression_above_mean = expression_level - mean_expression) |>
  group_by(archetype) |>
  slice_max(expression_above_mean, n = 15) |>
  select(archetype, gene) |>
  #write_csv(file = 'gabagook.csv')
  group_by(archetype) |>
  mutate(top_genes = str_c(gene, collapse = ', ')) |>
  select(archetype, top_genes) |>
  unique() |>
  mutate(top_genes = str_c(strwrap(top_genes, width = 40), collapse = "<br>")) 

t1dm_association <- read_csv(pea_filepath) |>
  dplyr::filter(labels == "T1DM")

beta_cells_archs <- read_csv(arch_filepath, col_names = str_c("PC_", 1:6)) |>
  select(PC_1:PC_3) |>
  rownames_to_column('index') |>
  mutate(is_datapoint = paste0('archetype_', index)) |>
  select(!index) |>
  # manually adding scores
  #add_column(`T1D score` = c(0.361, 0.759, 1.238, -0.109, -0.712, -0.899, -0.88),
  #           `T1D S level`= c(5,5,5,2,5,5,5)) |>
  left_join(y = top_genes, by = join_by('is_datapoint' == 'archetype')) |>
  left_join(y = enrichment, by = join_by('is_datapoint' == 'archetype')) |>
  left_join(y = t1dm_association, by = join_by('is_datapoint' == 'archetype'))


# really only need first 3 PCs for plotting
scores <- beta_cells_scores |>
  select(PC_1:PC_3) |>
  rownames_to_column("index") |>
  left_join(cell_labels) |>
  column_to_rownames('...1') |>
  select(!index) |>
  mutate(is_datapoint = 'yes') |>
  mutate(disease_pheno = as.character(factor(disease_pheno)))


custom_colors <- c("T1DM" = "blue", "Healthy" = "red")

plot_ly() |>
  add_trace(type='mesh3d', data=beta_cells_archs, 
            facecolor = colorRamp2::rand_color(n=8, luminosity = 'light'),
            x=~PC_1, y=~PC_2, z=~PC_3, alphahull = 0.02, opacity=0.1, 
            showlegend=T, legendgroup='1', name='show polytope', visible='legendonly') |>
  add_trace(data=beta_cells_archs, type='scatter3d',
            x=~PC_1, y=~PC_2, z=~PC_3, color=~is_datapoint, marker = list(size= 10,
                                                                          color='#8ecae6',
                                                                          line=list(
                                                                            color='#3c096c',
                                                                            width=2)
                                                                          ), 
            text = ~paste(is_datapoint, "<br>",
                          "T1D Enrichment: ", round(nes_scores, 3), 
                          '<br>Enrichment Null Z Score:', round(zscores, 3),
                          '<br> Top 15 Genes:', top_genes,
                          '<br> Top 5 Enriched Terms: <br>', top_terms),
            hoverinfo = 'text',
            legendgroup='2') |>
  add_trace(x=~PC_1, y=~PC_2, z=~PC_3, type="scatter3d", mode = 'markers',
            data = scores,
            color = ~disease_pheno,
            colors= custom_colors,
            marker = list(size = 3), 
            opacity = 0.5) |>
  layout(title= list(text=title, font = list(size=20), pad = list(t=5), y = 0.95),
         legend = list(bordercolor='#d4a373', borderwidth=1, title=list(text='Elements')))

}



